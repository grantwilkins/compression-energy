#include <assert.h>
#include <fcntl.h>
#include <hdf5.h>
#include <mpi.h>
#include <netcdf.h>
#include <netcdf_par.h>
#include <papi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>

#define MAX_ITERATIONS 10
#define MAX_POWERCAP_EVENTS 64

// PAPI related variables
int EventSet = PAPI_NULL;
int num_events = 0;
char event_names[MAX_POWERCAP_EVENTS][PAPI_MAX_STR_LEN];
int data_type[MAX_POWERCAP_EVENTS];

// Function to get current time in seconds
double get_time() {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return ts.tv_sec + ts.tv_nsec / 1e9;
}

// Function to calculate energy consumption from PAPI counters
void calculate_energy(long long *values, double *cpu_energy,
                      double *dram_energy) {
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

// Function to read raw binary file
void *read_raw_file(const char *filename, size_t *size) {
  FILE *file = fopen(filename, "rb");
  if (!file) {
    fprintf(stderr, "Error opening file: %s\n", filename);
    return NULL;
  }

  fseek(file, 0, SEEK_END);
  *size = ftell(file);
  fseek(file, 0, SEEK_SET);

  void *buffer = malloc(*size);
  if (!buffer) {
    fprintf(stderr, "Error allocating memory for file: %s\n", filename);
    fclose(file);
    return NULL;
  }

  size_t bytes_read = fread(buffer, 1, *size, file);
  if (bytes_read != *size) {
    fprintf(stderr, "Error reading file: %s\n", filename);
    free(buffer);
    fclose(file);
    return NULL;
  }

  fclose(file);
  return buffer;
}

// Function to perform I/O operation and measure performance
void perform_io(const char *method, const void *data, size_t data_size,
                const char *output_file, int mpi_rank, int mpi_size) {
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
  calculate_energy(values, &cpu_energy, &dram_energy);

  if (mpi_rank == 0) {
    printf("%s,%f,%f,%f\n", method, io_time, cpu_energy, dram_energy);
  }
}

int main(int argc, char *argv[]) {
  if (argc != 3) {
    fprintf(stderr, "Usage: %s <input_file> <output_dir>\n", argv[0]);
    return 1;
  }

  const char *input_file = argv[1];
  const char *output_dir = argv[2];

  // Initialize MPI
  int mpi_rank, mpi_size;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  // Initialize PAPI
  assert(PAPI_library_init(PAPI_VER_CURRENT) == PAPI_VER_CURRENT);
  assert(PAPI_create_eventset(&EventSet) == PAPI_OK);

  // Find powercap component and add events
  int numcmp = PAPI_num_components();
  int cid, powercap_cid = -1;
  const PAPI_component_info_t *cmpinfo = NULL;
  for (cid = 0; cid < numcmp; cid++) {
    cmpinfo = PAPI_get_component_info(cid);
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

  // Read input file
  size_t data_size;
  void *data = NULL;
  if (mpi_rank == 0) {
    data = read_raw_file(input_file, &data_size);
    if (!data) {
      MPI_Abort(MPI_COMM_WORLD, 1);
      return 1;
    }
  }
  MPI_Bcast(&data_size, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);

  // Allocate memory for all ranks
  if (mpi_rank != 0) {
    data = malloc(data_size);
  }
  MPI_Bcast(data, data_size, MPI_BYTE, 0, MPI_COMM_WORLD);

  // Perform I/O operations
  const char *methods[] = {"hdf5", "phdf5", "netcdf", "pnetcdf"};
  int num_methods = sizeof(methods) / sizeof(methods[0]);

  if (mpi_rank == 0) {
    printf("Method,Time(s),CPU_Energy(J),DRAM_Energy(J)\n");
  }

  char output_file[256];
  for (int i = 0; i < num_methods; i++) {
    snprintf(output_file, sizeof(output_file), "%s/%s_output.%s", output_dir,
             methods[i], (strstr(methods[i], "hdf5") ? "h5" : "nc"));
    for (int iter = 0; iter < MAX_ITERATIONS; iter++) {
      perform_io(methods[i], data, data_size, output_file, mpi_rank, mpi_size);
      MPI_Barrier(MPI_COMM_WORLD);
    }
  }

  // Clean up
  free(data);
  MPI_Finalize();
  return 0;
}