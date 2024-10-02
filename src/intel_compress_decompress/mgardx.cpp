#include "/home1/10191/gfw/MGARD/include/compress_x.hpp"
#include <algorithm>
#include <assert.h>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>

#include <cstring> // Added for strstr and strncpy
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

double calculate_mean(const std::vector<double> &data) {
  double sum = 0.0;
  for (double value : data) {
    sum += value;
  }
  return sum / data.size();
}

double calculate_std_dev(const std::vector<double> &data, double mean) {
  double sum_squared_diff = 0.0;
  for (double value : data) {
    double diff = value - mean;
    sum_squared_diff += diff * diff;
  }
  return std::sqrt(sum_squared_diff / (data.size() - 1));
}

bool within_confidence_interval(const std::vector<double> &data) {
  if (data.size() < 2)
    return false;
  double mean = calculate_mean(data);
  double std_dev = calculate_std_dev(data, mean);
  double margin_of_error =
      CONFIDENCE_LEVEL * (std_dev / std::sqrt(data.size()));
  double lower_bound = mean - margin_of_error;
  double upper_bound = mean + margin_of_error;

  for (double value : data) {
    if (value < lower_bound || value > upper_bound) {
      return false;
    }
  }
  return true;
}

int main(int argc, char *argv[]) {
  if (argc != 3) {
    std::cerr << "Usage: " << argv[0]
              << " <dataset_file> <relative_error_bound>" << std::endl;
    return 1;
  }

  const char *dataset_file = argv[1];
  double relative_error_bound = std::atof(argv[2]);
  const char *datadir = "/work2/10191/gfw/stampede3/"; // Update this path

  // PAPI initialization
  int EventSet = PAPI_NULL;
  long long values[MAX_POWERCAP_EVENTS];
  int num_events = 0;
  char event_names[MAX_POWERCAP_EVENTS][PAPI_MAX_STR_LEN];
  int data_type[MAX_POWERCAP_EVENTS];
  char event_descrs[MAX_POWERCAP_EVENTS][PAPI_MAX_STR_LEN];
  char units[MAX_POWERCAP_EVENTS][PAPI_MIN_STR_LEN];

  assert(PAPI_library_init(PAPI_VER_CURRENT) == PAPI_VER_CURRENT);
  assert(PAPI_create_eventset(&EventSet) == PAPI_OK);
  std::cout << "HERE" << std::endl;
  int numcmp = PAPI_num_components();
  int cid, powercap_cid = -1;
  const PAPI_component_info_t *cmpinfo;
  for (cid = 0; cid < numcmp; cid++) {
    cmpinfo = PAPI_get_component_info(cid);
    assert(cmpinfo != NULL);
    if (strstr(cmpinfo->name, "powercap")) {
      powercap_cid = cid;
      break;
    }
  }

  // Find and add powercap events
  int code = PAPI_NATIVE_MASK;
  PAPI_event_info_t evinfo;
  int r = PAPI_enum_cmp_event(&code, PAPI_ENUM_FIRST, powercap_cid);
  while (r == PAPI_OK && num_events < MAX_POWERCAP_EVENTS) {
    assert(PAPI_event_code_to_name(code, event_names[num_events]) == PAPI_OK);
    assert(PAPI_get_event_info(code, &evinfo) == PAPI_OK);
    strncpy(event_descrs[num_events], evinfo.long_descr,
            sizeof(event_descrs[0]) - 1);
    strncpy(units[num_events], evinfo.units, sizeof(units[0]) - 1);
    units[num_events][sizeof(units[0]) - 1] = '\0';
    data_type[num_events] = evinfo.data_type;

    if (PAPI_add_event(EventSet, code) != PAPI_OK) {
      break; // We've hit an event limit
    }
    num_events++;
    r = PAPI_enum_cmp_event(&code, PAPI_ENUM_EVENTS, powercap_cid);
  }
  std::cout << "HERE" << std::endl;

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

  // Prepare for compression
  mgard_x::Config config;
  config.dev_type =
      mgard_x::device_type::SERIAL; // Use CUDA for GPU compression

  std::vector<double> compression_times;
  std::vector<double> decompression_times;
  std::vector<double> compression_energy;
  std::vector<double> decompression_energy;

  int iteration = 0;
  bool confidence_interval_reached = false;

  std::cout << "Dataset: " << dataset_name << std::endl;
  std::cout << "Shape: ";
  for (auto dim : shape)
    std::cout << dim << " ";
  std::cout << std::endl;
  std::cout << "Number of elements: " << num_elements << std::endl;
  std::cout << "Data type: "
            << (dtype == mgard_x::data_type::Double ? "Double" : "Float")
            << std::endl;

  while (iteration < MAX_ITERATIONS && !confidence_interval_reached) {
    size_t max_compressed_size =
        num_elements *
        (dtype == mgard_x::data_type::Double ? sizeof(double) : sizeof(float));
    void *compressed_data = malloc(max_compressed_size);
    size_t compressed_size = max_compressed_size;
    // size_t compressed_size = 0;

    // Compression
    assert(PAPI_start(EventSet) == PAPI_OK);
    double start_time = get_time();

    mgard_x::compress_status_type compress_status =
        mgard_x::compress(shape.size(), // DIM
                          dtype,
                          shape, // Pass shape as std::vector<SIZE>
                          relative_error_bound,
                          0, // s = 0 for L2 norm
                          mgard_x::error_bound_type::REL, data, compressed_data,
                          compressed_size, config,
                          true // output_pre_allocated
        );

    double end_time = get_time();
    assert(PAPI_stop(EventSet, values) == PAPI_OK);

    if (compress_status != mgard_x::compress_status_type::Success) {
      std::cerr << "Compression failed" << std::endl;
      free(data);
      return 1;
    }

    double cpu_energy_compression = 0.0;
    for (int i = 0; i < num_events; i++) {
      if (strstr(event_names[i], "ENERGY_UJ") &&
          data_type[i] == PAPI_DATATYPE_UINT64) {
        cpu_energy_compression += values[i] / 1.0e6;
      }
    }

    compression_energy.push_back(cpu_energy_compression);
    compression_times.push_back(end_time - start_time);

    // Decompression
    void *decompressed_data = nullptr;

    assert(PAPI_start(EventSet) == PAPI_OK);
    start_time = get_time();

    mgard_x::compress_status_type decompress_status = mgard_x::decompress(
        compressed_data, compressed_size, decompressed_data, config,
        false // output_pre_allocated
    );
    end_time = get_time();
    assert(PAPI_stop(EventSet, values) == PAPI_OK);

    if (decompress_status != mgard_x::compress_status_type::Success) {
      std::cerr << "Decompression failed" << std::endl;
      free(data);
      free(compressed_data);
      return 1;
    }

    double cpu_energy_decompression = 0.0;
    for (int i = 0; i < num_events; i++) {
      if (strstr(event_names[i], "ENERGY_UJ") &&
          data_type[i] == PAPI_DATATYPE_UINT64) {
        cpu_energy_decompression += values[i] / 1.0e6;
      }
    }

    decompression_energy.push_back(cpu_energy_decompression);
    decompression_times.push_back(end_time - start_time);

    // Calculate metrics
    double max_error = 0.0, min_error = std::numeric_limits<double>::max();
    double max_rel_error = 0.0,
           min_rel_error = std::numeric_limits<double>::max();
    double max_pw_rel_error = 0.0,
           min_pw_rel_error = std::numeric_limits<double>::max();
    double mse = 0.0, sum_diff = 0.0;
    double data_min = std::numeric_limits<double>::max();
    double data_max = std::numeric_limits<double>::lowest();
    double sum = 0.0, sum_squared = 0.0;

    for (size_t i = 0; i < num_elements; i++) {
      double orig, decomp;
      if (dtype == mgard_x::data_type::Double) {
        orig = ((double *)data)[i];
        decomp = ((double *)decompressed_data)[i];
      } else {
        orig = ((float *)data)[i];
        decomp = ((float *)decompressed_data)[i];
      }
      double diff = decomp - orig;
      double abs_diff = std::abs(diff);
      double rel_diff = abs_diff / (std::abs(orig) + 1e-6);

      max_error = std::max(max_error, abs_diff);
      min_error = std::min(min_error, abs_diff);
      max_rel_error = std::max(max_rel_error, rel_diff);
      min_rel_error = std::min(min_rel_error, rel_diff);
      max_pw_rel_error =
          std::max(max_pw_rel_error, abs_diff / (std::abs(orig) + 1e-6));
      min_pw_rel_error =
          std::min(min_pw_rel_error, abs_diff / (std::abs(orig) + 1e-6));

      mse += diff * diff;
      sum_diff += diff;
      data_min = std::min(data_min, orig);
      data_max = std::max(data_max, orig);
      sum += orig;
      sum_squared += orig * orig;
    }

    double avg_difference = sum_diff / num_elements;
    double avg_error = std::sqrt(mse / num_elements);
    mse /= num_elements;
    double psnr = 20 * std::log10((data_max - data_min) / std::sqrt(mse));
    double nrmse = std::sqrt(mse) / (data_max - data_min);
    double compression_ratio =
        (double)(num_elements * (dtype == mgard_x::data_type::Double
                                     ? sizeof(double)
                                     : sizeof(float))) /
        compressed_size;
    double bit_rate = (double)compressed_size * 8 / num_elements;

    double value_range = data_max - data_min;
    double value_mean = sum / num_elements;
    double value_std =
        std::sqrt(sum_squared / num_elements - value_mean * value_mean);

    // Write metrics to CSV file
    FILE *csv_file = fopen("compression_metrics_mgard_x.csv", "a");
    if (csv_file) {
      fprintf(csv_file,
              "MGARD-X,%s,%e,%d,%f,%f,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%zu,%e,%"
              "e,%e,"
              "%e,%e,%e,%e,%e,%zu,%f,%zu,%f,%f,%f,%f\n",
              dataset_file, relative_error_bound, iteration,
              num_elements / compression_times.back(),   // compression_rate
              num_elements / decompression_times.back(), // decompression_rate
              avg_difference, avg_error,
              max_error - min_error,         // diff_range
              max_rel_error - min_rel_error, // error_range
              max_error, max_pw_rel_error, max_rel_error, min_error,
              min_pw_rel_error, min_rel_error, mse, num_elements, psnr, nrmse,
              data_max, value_mean, data_min, value_range, value_std, bit_rate,
              compressed_size, compression_ratio,
              num_elements * (dtype == mgard_x::data_type::Double
                                  ? sizeof(double)
                                  : sizeof(float)),
              compression_times.back(), decompression_times.back(),
              compression_energy.back(), decompression_energy.back());
      fclose(csv_file);
    } else {
      std::cerr << "Error opening CSV file" << std::endl;
    }

    // Clean up
    free(compressed_data);
    free(decompressed_data);

    iteration++;

    // Check if we've reached the confidence interval
    if (iteration >= 5) {
      confidence_interval_reached =
          within_confidence_interval(compression_times) &&
          within_confidence_interval(decompression_times);
    }
  }

  // Clean up
  free(data);
  PAPI_cleanup_eventset(EventSet);
  PAPI_destroy_eventset(&EventSet);

  return 0;
}
