#include "mgard/compress_x.hpp"
#include <algorithm>
#include <assert.h>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <hdf5.h>
#include <iostream>
#include <limits>
#include <netcdf.h>
#include <vector>

extern "C" {
#include <papi.h>
}

#define MAX_ITERATIONS 10
#define CONFIDENCE_LEVEL 1.96
#define MAX_POWERCAP_EVENTS 64

double get_time() {
  auto now = std::chrono::high_resolution_clock::now();
  auto duration = now.time_since_epoch();
  return std::chrono::duration_cast<std::chrono::duration<double>>(duration)
      .count();
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

    nc_create(output_file, NC_CLOBBER, &ncid);
    nc_def_dim(ncid, "size", data_size, &dimid);
    nc_def_var(ncid, "data", NC_BYTE, 1, &dimid, &varid);
    nc_enddef(ncid);

    nc_put_var_uchar(ncid, varid, static_cast<const unsigned char *>(data));

    nc_close(ncid);
  }
}

void delete_file(const char *filename) {
  if (remove(filename) != 0) {
    std::cerr << "Failed to delete file " << filename << std::endl;
  }
}

int main(int argc, char *argv[]) {
  if (argc != 3) {
    std::cerr << "Usage: " << argv[0]
              << " <dataset_file> <relative_error_bound>" << std::endl;
    return 1;
  }

  const char *dataset_file = argv[1];
  double relative_error_bound = std::atof(argv[2]);
  const char *datadir = "/work2/10191/gfw/stampede3/";
  const char *output_dir = "/work2/10191/gfw/stampede3/compressed/";

  // PAPI initialization
  int EventSet = PAPI_NULL;
  long long values[MAX_POWERCAP_EVENTS];
  int num_events = 0;
  char event_names[MAX_POWERCAP_EVENTS][PAPI_MAX_STR_LEN];
  int data_type[MAX_POWERCAP_EVENTS];

  assert(PAPI_library_init(PAPI_VER_CURRENT) == PAPI_VER_CURRENT);
  assert(PAPI_create_eventset(&EventSet) == PAPI_OK);

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

  std::vector<mgard_x::SIZE> shape;
  mgard_x::data_type dtype;
  void *data = nullptr;
  size_t num_elements = 0;
  std::string dataset_name;

  char full_path[1024];
  snprintf(full_path, sizeof(full_path), "%s%s", datadir, dataset_file);

  FILE *file = fopen(full_path, "rb");
  if (!file) {
    std::cerr << "Error opening file: " << full_path << std::endl;
    return 1;
  }

  if (strstr(dataset_file, "nyx") != NULL) {
    dataset_name = "NYX";
    shape = {512, 512, 512};
    dtype = mgard_x::data_type::Float;
    num_elements = 512ULL * 512ULL * 512ULL;
  } else if (strstr(dataset_file, "hacc") != NULL) {
    dataset_name = "HACC";
    shape = {1073726487};
    dtype = mgard_x::data_type::Float;
    num_elements = 1073726487;
  } else if (strstr(dataset_file, "s3d") != NULL) {
    dataset_name = "S3D";
    shape = {11, 500, 500, 500};
    dtype = mgard_x::data_type::Double;
    num_elements = 11ULL * 500ULL * 500ULL * 500ULL;
  } else if (strstr(dataset_file, "miranda") != NULL) {
    dataset_name = "Miranda";
    shape = {3072, 3072, 3072};
    dtype = mgard_x::data_type::Float;
    num_elements = 3072ULL * 3072ULL * 3072ULL;
  } else if (strstr(dataset_file, "cesm") != NULL) {
    dataset_name = "CESM";
    shape = {26, 1800, 3600};
    dtype = mgard_x::data_type::Float;
    num_elements = 26ULL * 1800ULL * 3600ULL;
  } else {
    std::cerr << "Unknown dataset: " << dataset_file << std::endl;
    fclose(file);
    return 1;
  }
  size_t data_size =
      (dtype == mgard_x::data_type::Float) ? sizeof(float) : sizeof(double);
  data = malloc(num_elements * data_size);
  size_t read_elements = fread(data, data_size, num_elements, file);
  fclose(file);

  if (read_elements != num_elements) {
    std::cerr << "Error reading data: expected " << num_elements
              << " elements, got " << read_elements << std::endl;
    free(data);
    return 1;
  }

  mgard_x::Config config;
  void *compressed_data = nullptr;
  size_t compressed_size = 0;

  // Compression
  assert(PAPI_start(EventSet) == PAPI_OK);
  double start_time = get_time();

  mgard_x::compress_status_type compress_status = mgard_x::compress(
      D, dtype, shape, relative_error_bound, 0, mgard_x::error_bound_type::REL,
      original_data, compressed_data, compressed_size, config);

  double end_time = get_time();
  assert(PAPI_stop(EventSet, values) == PAPI_OK);

  if (compress_status != mgard_x::compress_status_type::Success) {
    std::cerr << "Compression failed with status: "
              << static_cast<int>(compress_status) << std::endl;
    free(original_data);
    return 1;
  }

  double cpu_energy_compression = 0.0;
  for (int i = 0; i < num_events; i++) {
    if (strstr(event_names[i], "ENERGY_UJ") &&
        data_type[i] == PAPI_DATATYPE_UINT64) {
      cpu_energy_compression += values[i] / 1.0e6;
    }
  }

  // Perform I/O operations
  const char *methods[] = {"hdf5", "netcdf"};
  int num_methods = sizeof(methods) / sizeof(methods[0]);

  char output_file[256];
  FILE *fp = fopen("io_results_mgardx.csv", "a");
  for (int i = 0; i < num_methods; i++) {
    snprintf(output_file, sizeof(output_file), "%s%s_MGARDX_%g.%s", output_dir,
             dataset_file, relative_error_bound,
             (strstr(methods[i], "hdf5") ? "h5" : "nc"));
    for (int iter = 0; iter < MAX_ITERATIONS; iter++) {
      start_time = get_time();
      assert(PAPI_start(EventSet) == PAPI_OK);
      perform_io(methods[i], compressed_data, compressed_size, output_file);
      assert(PAPI_stop(EventSet, values) == PAPI_OK);
      end_time = get_time();
      delete_file(output_file);

      double cpu_energy = 0.0, dram_energy = 0.0;
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
        fprintf(fp, "%s,%s,MGARDX,%f,%d,%e,%e,%e\n", methods[i], dataset_file,
                relative_error_bound, iter, end_time - start_time, cpu_energy,
                dram_energy);
      } else {
        std::cerr << "Failed to open file io_results_mgardx.csv" << std::endl;
      }
    }
  }

  if (fp) {
    fclose(fp);
  }

  // Decompression (if needed for verification)
  void *decompressed_data = nullptr;
  mgard_x::compress_status_type decompress_status = mgard_x::decompress(
      compressed_data, compressed_size, decompressed_data, config);

  if (decompress_status != mgard_x::compress_status_type::Success) {
    std::cerr << "Decompression failed with status: "
              << static_cast<int>(decompress_status) << std::endl;
  }

  // Clean up
  free(original_data);
  free(compressed_data);
  free(decompressed_data);
  PAPI_cleanup_eventset(EventSet);
  PAPI_destroy_eventset(&EventSet);

  return 0;
}