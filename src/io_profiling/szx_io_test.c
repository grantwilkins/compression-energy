#include "szx.h"
#include "szx_rw.h"
#include <assert.h>
#include <hdf5.h>
#include <math.h>
#include <netcdf.h>
#include <papi.h>
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

void perform_io(const char *method, const void *data, size_t data_size,
                const char *output_file) {
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
  } else if (strcmp(method, "netcdf") == 0) {
    int ncid, varid, dimid;

    nc_create(output_file, NC_NETCDF4 | NC_CLOBBER, &ncid);
    nc_def_dim(ncid, "size", data_size, &dimid);
    nc_def_var(ncid, "data", NC_BYTE, 1, &dimid, &varid);
    nc_enddef(ncid);

    nc_put_var_uchar(ncid, varid, data);

    nc_close(ncid);
  }
}

void delete_file(const char *filename) {
  if (remove(filename) != 0) {
    fprintf(stderr, "Failed to delete file %s\n", filename);
  }
}

int main(int argc, char *argv[]) {
  if (argc != 3) {
    fprintf(stderr, "Usage: %s <dataset_file> <relative_error_bound>\n",
            argv[0]);
    return 1;
  }

  const char *dataset_file = argv[1];
  double relative_error_bound = atof(argv[2]);
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

  // Read dataset
  int status = 0;
  size_t nbEle;
  void *data = NULL;
  size_t dims[5] = {0};
  int data_type_szx;

  char full_path[1024];
  snprintf(full_path, sizeof(full_path), "%s%s", datadir, dataset_file);

  if (strstr(dataset_file, "nyx") != NULL) {
    dims[0] = 512;
    dims[1] = 512;
    dims[2] = 512;
    data = SZx_readFloatData(full_path, &nbEle, &status);
    data_type_szx = SZx_FLOAT;
  } else if (strstr(dataset_file, "hacc") != NULL) {
    dims[0] = 1073726487;
    data = SZx_readFloatData(full_path, &nbEle, &status);
    data_type_szx = SZx_FLOAT;
  } else if (strstr(dataset_file, "s3d") != NULL) {
    dims[0] = 11;
    dims[1] = 500;
    dims[2] = 500;
    dims[3] = 500;
    data = SZx_readDoubleData(full_path, &nbEle, &status);
    data_type_szx = SZx_DOUBLE;
  } else if (strstr(dataset_file, "miranda") != NULL) {
    dims[0] = 3072;
    dims[1] = 3072;
    dims[2] = 3072;
    data = SZx_readFloatData(full_path, &nbEle, &status);
    data_type_szx = SZx_FLOAT;
  } else if (strstr(dataset_file, "cesm") != NULL) {
    dims[0] = 26;
    dims[1] = 1800;
    dims[2] = 3600;
    data = SZx_readFloatData(full_path, &nbEle, &status);
    data_type_szx = SZx_FLOAT;
  } else {
    fprintf(stderr, "Unknown dataset %s\n", dataset_file);
    return 1;
  }

  if (status != SZx_SCES) {
    printf("Error: data file %s cannot be read!\n", dataset_file);
    return 1;
  }

  // Calculate data range and absolute error bound
  double data_min, data_max, data_range;
  if (data_type_szx == SZx_FLOAT) {
    float *float_data = (float *)data;
    data_min = data_max = float_data[0];
    for (size_t i = 1; i < nbEle; i++) {
      if (float_data[i] < data_min)
        data_min = float_data[i];
      if (float_data[i] > data_max)
        data_max = float_data[i];
    }
  } else { // SZx_DOUBLE
    double *double_data = (double *)data;
    data_min = data_max = double_data[0];
    for (size_t i = 1; i < nbEle; i++) {
      if (double_data[i] < data_min)
        data_min = double_data[i];
      if (double_data[i] > data_max)
        data_max = double_data[i];
    }
  }
  data_range = data_max - data_min;
  double absolute_error_bound = relative_error_bound * data_range;

  // Compress data
  size_t compressed_size;
  unsigned char *compressed_data = SZx_fast_compress_args(
      SZx_NO_BLOCK_FAST_CMPR, data_type_szx, data, &compressed_size, REL,
      absolute_error_bound, relative_error_bound, 0, 0, dims[4], dims[3],
      dims[2], dims[1], dims[0]);

  // Perform I/O operations
  const char *methods[] = {"hdf5", "netcdf"};
  int num_methods = sizeof(methods) / sizeof(methods[0]);

  char output_file[256];
  double start_time, end_time;
  double cpu_energy, dram_energy;
  FILE *fp = fopen("io_results_szx.csv", "a");
  for (int i = 0; i < num_methods; i++) {
    snprintf(output_file, sizeof(output_file), "%s%s_SZx_%g.%s", output_dir,
             dataset_file, relative_error_bound,
             (strstr(methods[i], "hdf5") ? "h5" : "nc"));
    for (int iter = 0; iter < MAX_ITERATIONS; iter++) {
      start_time = get_time();
      assert(PAPI_start(EventSet) == PAPI_OK);
      perform_io(methods[i], compressed_data, compressed_size, output_file);
      assert(PAPI_stop(EventSet, values) == PAPI_OK);
      end_time = get_time();
      delete_file(output_file);

      cpu_energy = 0.0;
      dram_energy = 0.0;
      for (int d = 0; d < num_events; d++) {
        if (strstr(event_names[d], "ENERGY_UJ")) {
          if (data_type[d] == PAPI_DATATYPE_UINT64) {
            if (strstr(event_names[d], "SUBZONE")) {
              dram_energy += values[d] / 1.0e6;
            } else {
              cpu_energy += values[d] / 1.0e6;
            }
          }
        }
      }
      if (fp) {
        fprintf(fp, "%s,%s,SZx,%f,%d,%e,%e,%e\n", methods[i], dataset_file,
                relative_error_bound, iter, end_time - start_time, cpu_energy,
                dram_energy);
      } else {
        fprintf(stderr, "Failed to open file io_results_szx.csv\n");
      }
    }
  }

  // Clean up
  free(data);
  free(compressed_data);
  free(values);
  PAPI_shutdown();

  return 0;
}
