#include <assert.h>
#include <hdf5.h>
#include <libpressio.h>
#include <libpressio_ext/io/posix.h>
#include <mpi.h>
#include <papi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define MAX_POWERCAP_EVENTS 64

// Function to get current time
double get_time() {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return ts.tv_sec + ts.tv_nsec / 1e9;
}

// Function to perform compression
void compress_data(struct pressio_compressor *compressor,
                   struct pressio_data *input_data,
                   struct pressio_data *compressed_data) {
  if (pressio_compressor_compress(compressor, input_data, compressed_data)) {
    fprintf(stderr, "Compression failed: %s\n",
            pressio_compressor_error_msg(compressor));
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
}

// Function to write compressed data to HDF5 file
void write_to_hdf5(const char *filename, void *data, size_t data_size) {
  hid_t file_id, dataset_id, dataspace_id;
  hsize_t dims[1] = {data_size};

  file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  dataspace_id = H5Screate_simple(1, dims, NULL);
  dataset_id = H5Dcreate2(file_id, "/compressed_data", H5T_NATIVE_UCHAR,
                          dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  H5Dwrite(dataset_id, H5T_NATIVE_UCHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

  H5Dclose(dataset_id);
  H5Sclose(dataspace_id);
  H5Fclose(file_id);
}

int main(int argc, char **argv) {
  int rank, size, node_rank, nodes, ranks_per_node;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // Get node information
  MPI_Comm node_comm;
  MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, rank, MPI_INFO_NULL,
                      &node_comm);
  MPI_Comm_rank(node_comm, &node_rank);
  MPI_Comm_size(node_comm, &ranks_per_node);
  nodes = size / ranks_per_node;

  if (argc != 5) {
    if (rank == 0) {
      fprintf(
          stderr,
          "Usage: %s <compressor> <dataset_file> <error_bound> <output_dir>\n",
          argv[0]);
    }
    MPI_Finalize();
    return 1;
  }

  const char *compressor_id = argv[1];
  const char *dataset_file = argv[2];
  double error_bound = atof(argv[3]);
  const char *output_dir = argv[4];

  // Initialize PAPI
  int EventSet = PAPI_NULL;
  long long *values;
  int num_events = 0;
  char event_names[MAX_POWERCAP_EVENTS][PAPI_MAX_STR_LEN];
  int data_type[MAX_POWERCAP_EVENTS];

  if (node_rank == 0) { // Only measure energy on one rank per node
    assert(PAPI_library_init(PAPI_VER_CURRENT) == PAPI_VER_CURRENT);
    assert(PAPI_create_eventset(&EventSet) == PAPI_OK);

    // Find powercap component
    int numcmp = PAPI_num_components();
    int cid, powercap_cid = -1;
    for (cid = 0; cid < numcmp; cid++) {
      const PAPI_component_info_t *cmpinfo = PAPI_get_component_info(cid);
      if (strstr(cmpinfo->name, "powercap")) {
        powercap_cid = cid;
        break;
      }
    }
    assert(powercap_cid != -1);

    // Add powercap events
    int code = PAPI_NATIVE_MASK;
    int r = PAPI_enum_cmp_event(&code, PAPI_ENUM_FIRST, powercap_cid);
    while (r == PAPI_OK) {
      PAPI_event_info_t evinfo;
      assert(PAPI_event_code_to_name(code, event_names[num_events]) == PAPI_OK);
      assert(PAPI_get_event_info(code, &evinfo) == PAPI_OK);
      data_type[num_events] = evinfo.data_type;
      if (PAPI_add_event(EventSet, code) != PAPI_OK)
        break;
      num_events++;
      r = PAPI_enum_cmp_event(&code, PAPI_ENUM_EVENTS, powercap_cid);
    }

    values = (long long *)calloc(num_events, sizeof(long long));
  }

  // Initialize LibPressio
  struct pressio *library = pressio_instance();
  struct pressio_compressor *compressor =
      pressio_get_compressor(library, compressor_id);
  if (compressor == NULL) {
    fprintf(stderr, "Failed to get compressor %s: %s\n", compressor_id,
            pressio_error_msg(library));
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  // Read input data (only on rank 0 of each node)
  struct pressio_data *input_data = NULL;
  size_t data_size = 0;
  if (node_rank == 0) {
    struct pressio_data *metadata = pressio_data_new_empty(
        pressio_float_dtype, 3, (size_t[]){512, 512, 512});
    input_data = pressio_io_data_path_read(metadata, dataset_file);
    if (input_data == NULL) {
      fprintf(stderr, "Failed to read dataset %s\n", dataset_file);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
    data_size = pressio_data_get_bytes(input_data);
    pressio_data_free(metadata);
  }

  // Broadcast data size to all ranks
  MPI_Bcast(&data_size, 1, MPI_UNSIGNED_LONG, 0, node_comm);

  // Allocate memory for data on all ranks
  void *data_buffer = malloc(data_size);

  // Broadcast data to all ranks within the node
  if (node_rank == 0) {
    memcpy(data_buffer, pressio_data_ptr(input_data, NULL), data_size);
  }
  MPI_Bcast(data_buffer, data_size, MPI_BYTE, 0, node_comm);

  // Create pressio_data object for each rank
  struct pressio_data *rank_input_data = pressio_data_new_move(
      pressio_float_dtype, data_buffer, 3, (size_t[]){512, 512, 512},
      pressio_data_libc_free_fn, NULL);

  // Set compressor options
  struct pressio_options *options = pressio_options_new();
  pressio_options_set_double(options, "pressio:abs", error_bound);
  pressio_compressor_set_options(compressor, options);

  // Prepare for compression
  struct pressio_data *compressed_data =
      pressio_data_new_empty(pressio_byte_dtype, 0, NULL);

  // Start energy measurement
  if (node_rank == 0) {
    assert(PAPI_start(EventSet) == PAPI_OK);
  }

  double start_time = get_time();

  // Perform compression
  compress_data(compressor, rank_input_data, compressed_data);

  // Write compressed data to HDF5 file
  char output_filename[256];
  snprintf(output_filename, sizeof(output_filename),
           "%s/compressed_data_node%d_rank%d.h5", output_dir,
           rank / ranks_per_node, rank % ranks_per_node);
  write_to_hdf5(output_filename, pressio_data_ptr(compressed_data, NULL),
                pressio_data_get_bytes(compressed_data));

  double end_time = get_time();

  // Stop energy measurement
  if (node_rank == 0) {
    assert(PAPI_stop(EventSet, values) == PAPI_OK);

    // Calculate energy consumption
    double cpu_energy = 0.0, dram_energy = 0.0;
    for (int i = 0; i < num_events; i++) {
      if (strstr(event_names[i], "ENERGY_UJ")) {
        if (data_type[i] == PAPI_DATATYPE_UINT64) {
          if (strstr(event_names[i], "SUBZONE")) {
            dram_energy += values[i] / 1.0e6;
          } else {
            cpu_energy += values[i] / 1.0e6;
          }
        }
      }
    }

    // Write energy stats to file
    char stats_filename[256];
    snprintf(stats_filename, sizeof(stats_filename),
             "%s/energy_stats_node%d.txt", output_dir, rank / ranks_per_node);
    FILE *stats_file = fopen(stats_filename, "w");
    if (stats_file) {
      fprintf(stats_file, "Node: %d\n", rank / ranks_per_node);
      fprintf(stats_file, "CPU Energy: %.6f J\n", cpu_energy);
      fprintf(stats_file, "DRAM Energy: %.6f J\n", dram_energy);
      fprintf(stats_file, "Total Energy: %.6f J\n", cpu_energy + dram_energy);
      fprintf(stats_file, "Compression Time: %.6f s\n", end_time - start_time);
      fclose(stats_file);
    } else {
      fprintf(stderr, "Failed to open stats file for writing\n");
    }
  }

  // Clean up
  pressio_data_free(rank_input_data);
  pressio_data_free(compressed_data);
  pressio_options_free(options);
  pressio_compressor_release(compressor);
  pressio_release(library);

  if (node_rank == 0) {
    free(values);
    PAPI_cleanup_eventset(EventSet);
    PAPI_destroy_eventset(&EventSet);
  }

  MPI_Finalize();
  return 0;
}