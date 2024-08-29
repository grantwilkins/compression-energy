#include <assert.h>
#include <hdf5.h>
#include <libpressio.h>
#include <libpressio_ext/io/posix.h>
#include <math.h>
#include <mpi.h>
#include <netcdf.h>
#include <papi.h>
#include <pnetcdf.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAX_ITERATIONS 10
#define MAX_POWERCAP_EVENTS 64

double get_time() {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return ts.tv_sec + ts.tv_nsec / 1e9;
}

void calculate_energy(long long *values, int num_events,
                      char event_names[][PAPI_MAX_STR_LEN], int *data_type,
                      double *cpu_energy, double *dram_energy) {
  *cpu_energy = 0.0;
  *dram_energy = 0.0;
  for (int i = 0; i < num_events; i++) {
    if (strstr(event_names[i], "ENERGY_UJ")) {
      if (data_type[i] == PAPI_DATATYPE_UINT64) {
        if (strstr(event_names[i], "PACKAGE_ENERGY") != NULL) {
          *cpu_energy += values[i] / 1.0e6; // Convert from uJ to J
        } else if (strstr(event_names[i], "DRAM_ENERGY") != NULL) {
          *dram_energy += values[i] / 1.0e6;
        }
      }
    }
  }
}

void perform_io(const char *method, const void *data, size_t data_size,
                const char *output_file, int mpi_rank, int mpi_size,
                int EventSet, int num_events,
                char event_names[][PAPI_MAX_STR_LEN], int *data_type) {
  double start_time, end_time, io_time;
  long long values[MAX_POWERCAP_EVENTS];
  double cpu_energy, dram_energy;

  // Start PAPI
  assert(PAPI_start(EventSet) == PAPI_OK);
  start_time = get_time();

  if (strcmp(method, "hdf5") == 0) {
    hid_t file_id, dataset_id, dataspace_id;
    hsize_t dims[1] = {data_size};

    file_id = H5Fcreate(output_file, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    dataspace_id = H5Screate_simple(1, dims, NULL);
    dataset_id = H5Dcreate2(file_id, "/data", H5T_NATIVE_UCHAR, dataspace_id,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    H5Dwrite(dataset_id, H5T_NATIVE_UCHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);
    H5Fclose(file_id);
  } else if (strcmp(method, "phdf5") == 0) {
    hid_t plist_id, file_id, dataset_id, dataspace_id;
    hsize_t dims[1] = {data_size};
    hsize_t chunk_dims[1] = {data_size / mpi_size};

    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
    file_id = H5Fcreate(output_file, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    H5Pclose(plist_id);

    dataspace_id = H5Screate_simple(1, dims, NULL);
    plist_id = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(plist_id, 1, chunk_dims);

    dataset_id = H5Dcreate2(file_id, "/data", H5T_NATIVE_UCHAR, dataspace_id,
                            H5P_DEFAULT, plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);

    hsize_t offset[1] = {mpi_rank * chunk_dims[0]};
    H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, chunk_dims,
                        NULL);

    hid_t memspace_id = H5Screate_simple(1, chunk_dims, NULL);

    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    H5Dwrite(dataset_id, H5T_NATIVE_UCHAR, memspace_id, dataspace_id, plist_id,
             (const unsigned char *)data + offset[0]);

    H5Pclose(plist_id);
    H5Sclose(memspace_id);
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);
    H5Fclose(file_id);
  } else if (strcmp(method, "netcdf") == 0) {
    int ncid, varid, dimid;

    nc_create(output_file, NC_CLOBBER, &ncid);
    nc_def_dim(ncid, "size", data_size, &dimid);
    nc_def_var(ncid, "data", NC_BYTE, 1, &dimid, &varid);
    nc_enddef(ncid);

    nc_put_var_uchar(ncid, varid, data);

    nc_close(ncid);
  } else if (strcmp(method, "pnetcdf") == 0) {
    int ncid, varid, dimid;
    MPI_Offset start, count;

    ncmpi_create(MPI_COMM_WORLD, output_file, NC_CLOBBER, MPI_INFO_NULL, &ncid);
    ncmpi_def_dim(ncid, "size", data_size, &dimid);
    ncmpi_def_var(ncid, "data", NC_BYTE, 1, &dimid, &varid);
    ncmpi_enddef(ncid);

    start = mpi_rank * (data_size / mpi_size);
    count = data_size / mpi_size;
    ncmpi_put_vara_uchar_all(ncid, varid, &start, &count,
                             (const unsigned char *)data + start);

    ncmpi_close(ncid);
  }

  end_time = get_time();
  assert(PAPI_stop(EventSet, values) == PAPI_OK);

  io_time = end_time - start_time;
  calculate_energy(values, num_events, event_names, data_type, &cpu_energy,
                   &dram_energy);

  if (mpi_rank == 0) {
    FILE *csv_file = fopen("io_results.csv", "a");
    if (csv_file) {
      fprintf(csv_file, "%s,%s,%f,%f,%f\n", method, output_file, io_time,
              cpu_energy, dram_energy);
      fclose(csv_file);
    } else {
      fprintf(stderr, "Error opening CSV file for writing results\n");
    }
  }
}

int main(int argc, char *argv[]) {
  if (argc != 4) {
    fprintf(stderr, "Usage: %s <compressor> <dataset_file> <error_bound>\n",
            argv[0]);
    return 1;
  }

  const char *compressor_id = argv[1];
  const char *dataset_file = argv[2];
  double relative_error_bound = atof(argv[3]);
  const char *datadir = "/work2/10191/gfw/stampede3/";
  const char *output_dir = "/work2/10191/gfw/stampede3/compressed/";

  // PAPI initialization
  int EventSet = PAPI_NULL;
  long long *values;
  int num_events = 0;
  char event_names[MAX_POWERCAP_EVENTS][PAPI_MAX_STR_LEN];
  int data_type[MAX_POWERCAP_EVENTS];

  assert(PAPI_library_init(PAPI_VER_CURRENT) == PAPI_VER_CURRENT);
  assert(PAPI_create_eventset(&EventSet) == PAPI_OK);

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
    if (PAPI_add_event(EventSet, code) != PAPI_OK)
      break;
    num_events++;
    r = PAPI_enum_cmp_event(&code, PAPI_ENUM_EVENTS, powercap_cid);
  }

  values = (long long *)calloc(num_events, sizeof(long long));

  // LibPressio initialization
  struct pressio *library = pressio_instance();
  struct pressio_compressor *compressor =
      pressio_get_compressor(library, compressor_id);
  if (compressor == NULL) {
    fprintf(stderr, "Failed to get compressor %s: %s\n", compressor_id,
            pressio_error_msg(library));
    pressio_release(library);
    return 1;
  }

  // Read input data
  char full_path[1024];
  snprintf(full_path, sizeof(full_path), "%s%s", datadir, dataset_file);
  struct pressio_data *input_data = pressio_io_data_path_read(NULL, full_path);
  if (input_data == NULL) {
    fprintf(stderr, "Failed to read dataset %s\n", full_path);
    pressio_compressor_release(compressor);
    pressio_release(library);
    return 1;
  }

  // Calculate absolute error bound
  size_t num_elements = pressio_data_num_elements(input_data);
  void *data_ptr = pressio_data_ptr(input_data, NULL);
  double data_min, data_max, data_range;

  if (pressio_data_dtype(input_data) == pressio_float_dtype) {
    float *float_data = (float *)data_ptr;
    data_min = data_max = float_data[0];
    for (size_t i = 1; i < num_elements; i++) {
      if (float_data[i] < data_min)
        data_min = float_data[i];
      if (float_data[i] > data_max)
        data_max = float_data[i];
    }
  } else {
    double *double_data = (double *)data_ptr;
    data_min = data_max = double_data[0];
    for (size_t i = 1; i < num_elements; i++) {
      if (double_data[i] < data_min)
        data_min = double_data[i];
      if (double_data[i] > data_max)
        data_max = double_data[i];
    }
  }
  data_range = data_max - data_min;
  double absolute_error_bound = relative_error_bound * data_range;

  // Configure compressor
  struct pressio_options *options = pressio_options_new();
  pressio_options_set_double(options, "pressio:abs", absolute_error_bound);
  if (pressio_compressor_set_options(compressor, options)) {
    fprintf(stderr, "Failed to set compressor options: %s\n",
            pressio_compressor_error_msg(compressor));
    return 1;
  }
  pressio_options_free(options);

  // Compress data
  struct pressio_data *compressed_data =
      pressio_data_new_empty(pressio_byte_dtype, 0, NULL);
  if (pressio_compressor_compress(compressor, input_data, compressed_data)) {
    fprintf(stderr, "Compression failed: %s\n",
            pressio_compressor_error_msg(compressor));
    return 1;
  }

  // Get compressed data pointer and size
  void *compressed_ptr = pressio_data_ptr(compressed_data, NULL);
  size_t compressed_size = pressio_data_get_bytes(compressed_data);

  // Perform I/O operations
  const char *methods[] = {"hdf5", "phdf5", "netcdf", "pnetcdf"};
  int num_methods = sizeof(methods) / sizeof(methods[0]);

  // Initialize MPI
  int mpi_rank, mpi_size;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  char output_file[256];
  for (int i = 0; i < num_methods; i++) {
    snprintf(output_file, sizeof(output_file), "%s%s_%s_%g.%s", output_dir,
             dataset_file, compressor_id, relative_error_bound,
             (strstr(methods[i], "hdf5") ? "h5" : "nc"));
    for (int iter = 0; iter < MAX_ITERATIONS; iter++) {
      void perform_io(const char *method, const void *data, size_t data_size,
                      const char *output_file, int mpi_rank, int mpi_size,
                      int EventSet, int num_events,
                      char event_names[][PAPI_MAX_STR_LEN], int *data_type)
          MPI_Barrier(MPI_COMM_WORLD);
    }
  }

  // Clean up
  pressio_data_free(input_data);
  pressio_data_free(compressed_data);
  pressio_compressor_release(compressor);
  pressio_release(library);
  free(values);
  MPI_Finalize();
  PAPI_shutdown();

  return 0;
}
