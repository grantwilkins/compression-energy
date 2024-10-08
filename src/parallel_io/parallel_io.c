#include <assert.h>
#include <hdf5.h>
#include <libpressio.h>
#include <libpressio_ext/io/posix.h>
#include <mpi.h>
#include <netcdf.h>
#include <papi.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#define MAX_POWERCAP_EVENTS 64
#define NUM_ITERATIONS 10

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

// Function to write compressed data to NetCDF file
void write_to_netcdf(const char *output_file, void *data, size_t data_size) {
  int ncid, varid, dimid;

  nc_create(output_file, NC_NETCDF4 | NC_CLOBBER, &ncid);
  nc_def_dim(ncid, "size", data_size, &dimid);
  nc_def_var(ncid, "data", NC_BYTE, 1, &dimid, &varid);
  nc_enddef(ncid);

  nc_put_var_uchar(ncid, varid, data);

  nc_close(ncid);
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

  if (argc != 4) {
    if (rank == 0) {
      fprintf(stderr, "Usage: %s <compressor> <dataset_file> <error_bound>\n",
              argv[0]);
    }
    MPI_Finalize();
    return 1;
  }

  const char *compressor_id = argv[1];
  const char *dataset_file = argv[2];
  double error_bound = atof(argv[3]);
  const char *datadir = "/work2/10191/gfw/stampede3/";
  const char *output_dir = "/work2/10191/gfw/stampede3/compressed/";

  // Initialize PAPI
  int EventSet = PAPI_NULL;
  long long *values;
  int num_events = 0;
  char event_names[MAX_POWERCAP_EVENTS][PAPI_MAX_STR_LEN];
  int data_type[MAX_POWERCAP_EVENTS];

  if (node_rank == 0) { // Only initialize PAPI on one rank per node
    if (PAPI_library_init(PAPI_VER_CURRENT) != PAPI_VER_CURRENT) {
      fprintf(stderr, "PAPI library init error!\n");
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
    PAPI_CHECK(PAPI_thread_init(pthread_self),
               "Error initializing PAPI thread support");

    PAPI_CHECK(PAPI_create_eventset(&EventSet),
               "Error creating PAPI event set");

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
    if (powercap_cid == -1) {
      fprintf(stderr, "No powercap component found\n");
      MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Add powercap events
    int code = PAPI_NATIVE_MASK;
    int r = PAPI_enum_cmp_event(&code, PAPI_ENUM_FIRST, powercap_cid);
    while (r == PAPI_OK) {
      PAPI_CHECK(PAPI_add_event(EventSet, code), "Error adding PAPI event");
      r = PAPI_enum_cmp_event(&code, PAPI_ENUM_EVENTS, powercap_cid);
    }
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
  char full_path[1024];
  snprintf(full_path, sizeof(full_path), "%s%s", datadir, dataset_file);
  struct pressio_data *input_data = NULL;
  size_t data_size = 0;
  if (node_rank == 0) {
    struct pressio_data *metadata = pressio_data_new_empty(
        pressio_float_dtype, 3, (size_t[]){512, 512, 512});
    input_data = pressio_io_data_path_read(metadata, full_path);
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
  MPI_Barrier(MPI_COMM_WORLD);
  // Broadcast data to all ranks within the node
  if (node_rank == 0) {
    memcpy(data_buffer, pressio_data_ptr(input_data, NULL), data_size);
  }
  int mpi_err;
  mpi_err = MPI_Bcast(data_buffer, data_size, MPI_BYTE, 0, node_comm);
  if (mpi_err != MPI_SUCCESS) {
    fprintf(stderr, "MPI_Bcast error on rank %d\n", rank);
    MPI_Abort(MPI_COMM_WORLD, mpi_err);
  }
  printf("bcast data: %d %d\n", node_rank, rank);
  // Create pressio_data object for each rank
  struct pressio_data *rank_input_data = pressio_data_new_move(
      pressio_float_dtype, data_buffer, 3, (size_t[]){512, 512, 512},
      pressio_data_libc_free_fn, NULL);

  // Set compressor options
  struct pressio_options *options = pressio_options_new();
  pressio_options_set_double(options, "pressio:rel", error_bound);
  pressio_compressor_set_options(compressor, options);

  // Prepare for compression
  struct pressio_data *compressed_data =
      pressio_data_new_empty(pressio_byte_dtype, 0, NULL);

  // Compression phase
  double compress_start_time, compress_end_time;
  long long *compress_values = NULL;

  if (node_rank == 0) {
    num_events = PAPI_num_events(EventSet);
    compress_values = (long long *)calloc(num_events, sizeof(long long));
    assert(PAPI_start(EventSet) == PAPI_OK);
  }

  compress_start_time = get_time();

  // Perform compression
  printf("Compress before %d %d\n", node_rank, rank);
  compress_data(compressor, rank_input_data, compressed_data);
  printf("Compress after %d %d\n", node_rank, rank);

  compress_end_time = get_time();

  // Stop energy measurement for compression
  if (node_rank == 0) {
    assert(PAPI_stop(EventSet, compress_values) == PAPI_OK);
  }

  // Synchronize all processes after compression
  MPI_Barrier(MPI_COMM_WORLD);

  // I/O phase
  const char *io_methods[] = {"hdf5", "netcdf"};
  int num_methods = sizeof(io_methods) / sizeof(io_methods[0]);

  for (int method = 0; method < num_methods; method++) {
    for (int iteration = 0; iteration < NUM_ITERATIONS; iteration++) {
      double write_start_time, write_end_time;
      long long *write_values = NULL;

      // Start energy measurement for writing
      if (node_rank == 0) {
        write_values = (long long *)calloc(num_events, sizeof(long long));
        assert(PAPI_start(EventSet) == PAPI_OK);
      }

      write_start_time = get_time();

      // Write compressed data to file
      char output_filename[256];
      snprintf(output_filename, sizeof(output_filename),
               "%s/%s_%s_%g_%d_%d_%s_%d.%s", output_dir, dataset_file,
               compressor_id, error_bound, rank, node_rank, io_methods[method],
               iteration, (method == 0 ? "h5" : "nc"));

      if (method == 0) {
        write_to_hdf5(output_filename, pressio_data_ptr(compressed_data, NULL),
                      pressio_data_get_bytes(compressed_data));
      } else {
        write_to_netcdf(output_filename,
                        pressio_data_ptr(compressed_data, NULL),
                        pressio_data_get_bytes(compressed_data));
      }

      write_end_time = get_time();

      // Stop energy measurement for writing
      if (node_rank == 0) {
        assert(PAPI_stop(EventSet, write_values) == PAPI_OK);

        // Calculate energy consumption for writing
        double write_cpu_energy = 0.0, write_dram_energy = 0.0;
        for (int i = 0; i < num_events; i++) {
          if (strstr(event_names[i], "ENERGY_UJ")) {
            if (data_type[i] == PAPI_DATATYPE_UINT64) {
              if (strstr(event_names[i], "SUBZONE")) {
                write_dram_energy += write_values[i] / 1.0e6;
              } else {
                write_cpu_energy += write_values[i] / 1.0e6;
              }
            }
          }
        }

        // Write energy stats to file
        char stats_filename[256];
        snprintf(stats_filename, sizeof(stats_filename),
                 "parallel_write_split.csv");
        FILE *stats_file = fopen(stats_filename, "a");
        if (stats_file) {
          fprintf(stats_file, "%s,%s,%d,%d,%e,%s,%d,%e,%e,%e,%zu\n",
                  dataset_file, compressor_id, node_rank, rank, error_bound,
                  io_methods[method], iteration,
                  write_end_time - write_start_time, write_cpu_energy,
                  write_dram_energy, pressio_data_get_bytes(compressed_data));
          fclose(stats_file);
        } else {
          fprintf(stderr, "Failed to open stats file for writing\n");
        }
      }

      // Synchronize all processes after writing
      MPI_Barrier(MPI_COMM_WORLD);

      // Delete the file after writing
      if (rank == 0) {
        if (remove(output_filename) != 0) {
          fprintf(stderr, "Error deleting file: %s\n", output_filename);
        }
      }

      // Synchronize all processes after deleting
      MPI_Barrier(MPI_COMM_WORLD);
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
