#include <assert.h>
#include <hdf5.h>
#include <libpressio.h>
#include <libpressio_ext/io/posix.h>
#include <mpi.h>
#include <netcdf.h>
#include <papi.h>
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

void handle_error(int retval) {
  printf("PAPI error %d: %s\n", retval, PAPI_strerror(retval));
  exit(1);
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

  if (PAPI_library_init(PAPI_VER_CURRENT) != PAPI_VER_CURRENT) {
    handle_error(1);
  }
  if (PAPI_create_eventset(&EventSet) != PAPI_OK) {
    handle_error(1);
  }

  // Find powercap component and add events
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

  int code = PAPI_NATIVE_MASK;
  int r = PAPI_enum_cmp_event(&code, PAPI_ENUM_FIRST, powercap_cid);
  while (r == PAPI_OK) {
    PAPI_event_info_t evinfo;
    assert(PAPI_event_code_to_name(code, event_names[num_events]) == PAPI_OK);
    assert(PAPI_get_event_info(code, &evinfo) == PAPI_OK);
    data_type[num_events] = evinfo.data_type;
    if (PAPI_add_event(EventSet, code) != PAPI_OK) {
      break;
    }
    num_events++;
    r = PAPI_enum_cmp_event(&code, PAPI_ENUM_EVENTS, powercap_cid);
  }

  long long *compress_values =
      (long long *)calloc(num_events, sizeof(long long));
  long long *write_values = (long long *)calloc(num_events, sizeof(long long));

  // Initialize LibPressio
  struct pressio *library = pressio_instance();
  struct pressio_compressor *compressor = NULL;
  pressio_get_compressor(library, compressor_id);
  if (strcmp(compressor_id, "None") != 0) {
    compressor = pressio_get_compressor(library, compressor_id);
    if (compressor == NULL) {
      fprintf(stderr, "Failed to get compressor %s: %s\n", compressor_id,
              pressio_error_msg(library));
      pressio_release(library);
      return 1;
    }
  }

  struct pressio_data *input_data = NULL;
  size_t data_size = 0;
  size_t dims[4] = {0};
  size_t ndims = 0;
  enum pressio_dtype dtype;

  if (node_rank == 0) {
    if (strstr(dataset_file, "nyx") != NULL) {
      dims[0] = 512;
      dims[1] = 512;
      dims[2] = 512;
      ndims = 3;
      dtype = pressio_float_dtype;
    } else if (strstr(dataset_file, "hacc") != NULL) {
      dims[0] = 1073726487;
      ndims = 1;
      dtype = pressio_float_dtype;
    } else if (strstr(dataset_file, "s3d") != NULL) {
      dims[0] = 11;
      dims[1] = 500;
      dims[2] = 500;
      dims[3] = 500;
      ndims = 4;
      dtype = pressio_double_dtype;
    } else if (strstr(dataset_file, "miranda") != NULL) {
      dims[0] = 3072;
      dims[1] = 3072;
      dims[2] = 3072;
      ndims = 3;
      dtype = pressio_float_dtype;
    } else if (strstr(dataset_file, "cesm") != NULL) {
      dims[0] = 26;
      dims[1] = 1800;
      dims[2] = 3600;
      ndims = 3;
      dtype = pressio_float_dtype;
    } else {
      fprintf(stderr, "Unknown dataset %s\n", dataset_file);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }

    struct pressio_data *metadata = pressio_data_new_empty(dtype, ndims, dims);
    char full_path[1024];
    snprintf(full_path, sizeof(full_path), "%s%s", datadir, dataset_file);
    input_data = pressio_io_data_path_read(metadata, full_path);
    if (input_data == NULL) {
      fprintf(stderr, "Failed to read dataset %s\n", dataset_file);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
    data_size = pressio_data_get_bytes(input_data);
    pressio_data_free(metadata);
  }

  // Broadcast metadata to all ranks in the node
  MPI_Bcast(&data_size, 1, MPI_UNSIGNED_LONG, 0, node_comm);
  MPI_Bcast(&ndims, 1, MPI_UNSIGNED_LONG, 0, node_comm);
  MPI_Bcast(dims, 4, MPI_UNSIGNED_LONG, 0, node_comm);
  int dtype_int;
  if (node_rank == 0) {
    dtype_int = (int)dtype;
  }
  MPI_Bcast(&dtype_int, 1, MPI_INT, 0, node_comm);
  dtype = (enum pressio_dtype)dtype_int;

  printf("Rank %d: data_size = %zu, ndims = %zu, dims = {%zu, %zu, %zu}, "
         "dtype = %d\n",
         rank, data_size, ndims, dims[0], dims[1], dims[2], dtype);

  // Allocate memory for data on all ranks
  void *data_buffer = malloc(data_size);
  if (data_buffer == NULL) {
    fprintf(stderr, "Failed to allocate memory on rank %d\n", rank);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  // Broadcast data to all ranks within the node
  if (node_rank == 0) {
    memcpy(data_buffer, pressio_data_ptr(input_data, NULL), data_size);
  }
  MPI_Bcast(data_buffer, data_size, MPI_BYTE, 0, node_comm);
  MPI_Barrier(MPI_COMM_WORLD);
  // Create pressio_data object for each rank
  struct pressio_data *rank_input_data = pressio_data_new_move(
      dtype, data_buffer, ndims, dims, pressio_data_libc_free_fn, NULL);
  MPI_Barrier(MPI_COMM_WORLD);
  // Set compressor options
  struct pressio_options *options = pressio_options_new();
  pressio_options_set_double(options, "pressio:rel", error_bound);
  pressio_compressor_set_options(compressor, options);

  // Prepare for compression
  struct pressio_data *compressed_data =
      pressio_data_new_empty(pressio_byte_dtype, 0, NULL);
  MPI_Barrier(MPI_COMM_WORLD);
  // Compression phase
  double compress_start_time = get_time();

  // Perform compression
  compress_data(compressor, rank_input_data, compressed_data);

  double compress_end_time = get_time();
  double compress_time = compress_end_time - compress_start_time;

  // I/O phase
  const char *io_methods[] = {"hdf5", "netcdf"};
  int num_methods = sizeof(io_methods) / sizeof(io_methods[0]);

  for (int method = 0; method < num_methods; method++) {
    for (int iteration = 0; iteration < NUM_ITERATIONS; iteration++) {
      double write_start_time = get_time();

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

      double write_end_time = get_time();
      double write_time = write_end_time - write_start_time;

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

      // Collect and report statistics
      if (rank == 0) {
        char stats_filename[256];
        snprintf(stats_filename, sizeof(stats_filename),
                 "parallel_write_split.csv");
        FILE *stats_file = fopen(stats_filename, "a");
        if (stats_file) {
          fprintf(stats_file, "%s,%s,%d,%d,%e,%s,%d,%e,%e,%zu\n", dataset_file,
                  compressor_id, nodes, size, error_bound, io_methods[method],
                  iteration, compress_time, write_time,
                  pressio_data_get_bytes(compressed_data));
          fclose(stats_file);
        } else {
          fprintf(stderr, "Failed to open stats file for writing\n");
        }
      }
    }
  }

  // Stop PAPI counters
  if (PAPI_stop(EventSet, values) != PAPI_OK) {
    handle_error(1);
  }

  if (rank == 0) {
    printf("PAPI measured %lld instructions\n", values[0]);
  }

  // Clean up
  pressio_data_free(rank_input_data);
  pressio_data_free(compressed_data);
  pressio_options_free(options);
  pressio_compressor_release(compressor);
  pressio_release(library);

  PAPI_cleanup_eventset(EventSet);
  PAPI_destroy_eventset(&EventSet);

  MPI_Finalize();
  return 0;
}