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

// Function to write data to HDF5 file
void write_to_hdf5(const char *filename, void *data, size_t data_size) {
  hid_t file_id, dataset_id, dataspace_id;
  hsize_t dims[1] = {data_size};

  file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  dataspace_id = H5Screate_simple(1, dims, NULL);
  dataset_id = H5Dcreate2(file_id, "/data", H5T_NATIVE_UCHAR, dataspace_id,
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  H5Dwrite(dataset_id, H5T_NATIVE_UCHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

  H5Dclose(dataset_id);
  H5Sclose(dataspace_id);
  H5Fclose(file_id);
}

// Function to write data to NetCDF file
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
  double relative_error_bound = atof(argv[3]);
  const char *datadir = "/work2/10191/gfw/stampede3/";
  const char *output_dir = "/work2/10191/gfw/stampede3/compressed/";
  int i;

  // Initialize PAPI
  int EventSet = PAPI_NULL;
  int num_events = 0;
  char event_names[MAX_POWERCAP_EVENTS][PAPI_MAX_STR_LEN];
  int data_type[MAX_POWERCAP_EVENTS];

  if (node_rank == 0) {
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
  }

  long long *transmit_values = NULL;
  long long *write_values = NULL;
  if (node_rank == 0) {
    transmit_values = (long long *)calloc(num_events, sizeof(long long));
    write_values = (long long *)calloc(num_events, sizeof(long long));
  }

  // Initialize LibPressio
  struct pressio *library = pressio_instance();
  struct pressio_data *input_data = NULL;
  size_t data_size = 0;
  size_t dims[4] = {0};
  size_t ndims = 0;
  enum pressio_dtype dtype;

  // Read data on one rank per node
  if (rank == 0) {
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
      fprintf(stderr, "Rank %d failed to read dataset %s\n", rank,
              dataset_file);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
    data_size = pressio_data_get_bytes(input_data);
    pressio_data_free(metadata);
  }

  MPI_Bcast(&data_size, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
  MPI_Bcast(&ndims, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
  MPI_Bcast(dims, 4, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
  int dtype_int;
  if (rank == 0) {
    dtype_int = (int)dtype;
  }
  MPI_Bcast(&dtype_int, 1, MPI_INT, 0, MPI_COMM_WORLD);
  dtype = (enum pressio_dtype)dtype_int;

  // Allocate memory for data on all ranks
  void *data_buffer = malloc(data_size);
  if (data_buffer == NULL) {
    fprintf(stderr, "Failed to allocate memory on rank %d\n", rank);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  // Transmission phase
  double transmit_start_time = get_time();

  if (node_rank == 0) {
    if (PAPI_start(EventSet) != PAPI_OK) {
      handle_error(1);
    }
  }

  // Broadcast data to all ranks
  if (rank == 0) {
    memcpy(data_buffer, pressio_data_ptr(input_data, NULL), data_size);
  }
  MPI_Bcast(data_buffer, data_size, MPI_BYTE, 0, MPI_COMM_WORLD);

  MPI_Barrier(MPI_COMM_WORLD);
  if (node_rank == 0) {
    if (PAPI_stop(EventSet, transmit_values) != PAPI_OK) {
      handle_error(1);
    }
  }
  double transmit_end_time = get_time();
  double transmit_time = transmit_end_time - transmit_start_time;

  MPI_Barrier(MPI_COMM_WORLD);

  double cpu_energy_transmission = 0.0;
  if (node_rank == 0) {
    for (i = 0; i < num_events; i++) {
      if (strstr(event_names[i], "ENERGY_UJ")) {
        if (data_type[i] == PAPI_DATATYPE_UINT64) {
          if (!strstr(event_names[i], "SUBZONE")) {
            cpu_energy_transmission += transmit_values[i] / 1.0e6;
          }
        }
      }
    }
  }

  printf("Rank %d received %zu data\n", rank, data_size);

  double transmit_time_max = 0.0;
  MPI_Reduce(&transmit_time, &transmit_time_max, 1, MPI_DOUBLE, MPI_MAX, 0,
             node_comm);
  
  
  // I/O phase
  const char *io_methods[] = {"hdf5"};
  int num_methods = sizeof(io_methods) / sizeof(io_methods[0]);

  for (int method = 0; method < num_methods; method++) {
    for (int iteration = 0; iteration < NUM_ITERATIONS; iteration++) {
      double write_start_time = get_time();
      // Write data to file
      char output_filename[256];
      snprintf(output_filename, sizeof(output_filename),
               "%s/%s_%s_%g_%d_%d_%s_%d.%s", output_dir, dataset_file,
               compressor_id, relative_error_bound, rank, node_rank,
               io_methods[method], iteration, (method == 0 ? "h5" : "nc"));
      if (node_rank == 0) {
        if (PAPI_start(EventSet) != PAPI_OK) {
          handle_error(1);
        }
      }
      if (method == 0) {
        write_to_hdf5(output_filename, data_buffer, data_size);
      } else {
        write_to_netcdf(output_filename, data_buffer, data_size);
      }

      printf("Rank %d here\n", rank);

      MPI_Barrier(MPI_COMM_WORLD);
      if (node_rank == 0) {
        if (PAPI_stop(EventSet, write_values) != PAPI_OK) {
          handle_error(1);
        }
      }

      double write_end_time = get_time();
      double write_time = write_end_time - write_start_time;

      double cpu_energy_writing = 0.0;
      if (node_rank == 0) {
        for (i = 0; i < num_events; i++) {
          if (strstr(event_names[i], "ENERGY_UJ")) {
            if (data_type[i] == PAPI_DATATYPE_UINT64) {
              if (!strstr(event_names[i], "SUBZONE")) {
                cpu_energy_writing += write_values[i] / 1.0e6;
              }
            }
          }
        }
      }

      double write_time_max = 0.0;
      MPI_Reduce(&write_time, &write_time_max, 1, MPI_DOUBLE, MPI_MAX, 0,
                 node_comm);

      // Delete the file after writing
      if (remove(output_filename) != 0) {
        fprintf(stderr, "Error deleting file: %s\n", output_filename);
      }

      // Synchronize all processes after deleting
      MPI_Barrier(MPI_COMM_WORLD);

      // Collect and report statistics
      if (node_rank == 0) {
        char stats_filename[256];
        snprintf(stats_filename, sizeof(stats_filename),
                 "parallel_write_experiment.csv");
        FILE *stats_file = fopen(stats_filename, "a");
        if (stats_file) {
          int node_num = rank / ranks_per_node;
          fprintf(stats_file, "%s,%s,%d,%d,%d,%e,%s,%d,%e,%e,%e,%e,%zu\n",
                  dataset_file, compressor_id, node_num, nodes, size,
                  relative_error_bound, io_methods[method], iteration,
                  transmit_time_max, write_time_max, cpu_energy_transmission,
                  cpu_energy_writing, data_size);
          fclose(stats_file);
        } else {
          fprintf(stderr, "Failed to open stats file for writing\n");
        }
      }
    }
  }

  // Clean up
  free(data_buffer);
  if (rank == 0) {
    pressio_data_free(input_data);
  }
  pressio_release(library);
  if (node_rank == 0) {
    free(transmit_values);
    free(write_values);
    PAPI_cleanup_eventset(EventSet);
    PAPI_destroy_eventset(&EventSet);
  }

  MPI_Finalize();
  return 0;
}
