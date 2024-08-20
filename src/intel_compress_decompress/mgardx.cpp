#include <algorithm>
#include <assert.h>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <compress_x.hpp>
#include <papi.h>
#include <vector>

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
  const char *datadir = "/path/to/datasets/"; // Update this path

  // PAPI initialization
  int EventSet = PAPI_NULL;
  long long values[MAX_POWERCAP_EVENTS];
  int num_events = 0;
  char event_names[MAX_POWERCAP_EVENTS][PAPI_MAX_STR_LEN];
  int data_type[MAX_POWERCAP_EVENTS];

  assert(PAPI_library_init(PAPI_VER_CURRENT) == PAPI_VER_CURRENT);
  assert(PAPI_create_eventset(&EventSet) == PAPI_OK);

  // Find and add powercap events
  int code = PAPI_NATIVE_MASK;
  PAPI_event_info_t info;
  while (PAPI_enum_cmp_event(&code, PAPI_ENUM_FIRST, 0) == PAPI_OK) {
    if (PAPI_get_event_info(code, &info) == PAPI_OK) {
      if (strstr(info.symbol, "powercap")) {
        if (PAPI_add_event(EventSet, code) == PAPI_OK) {
          strncpy(event_names[num_events], info.symbol, PAPI_MAX_STR_LEN);
          data_type[num_events] = info.data_type;
          num_events++;
        }
      }
    }
  }

  if (strstr(dataset_file, "nyx") != NULL) {
    dataset_name = "NYX";
    shape = {512, 512, 512};
    dtype = mgard_x::data_type::Float;
    num_elements = 512 * 512 * 512;
  } else if (strstr(dataset_file, "hacc") != NULL) {
    dataset_name = "HACC";
    shape = {1073726487};
    dtype = mgard_x::data_type::Float;
    num_elements = 1073726487;
  } else if (strstr(dataset_file, "s3d") != NULL) {
    dataset_name = "S3D";
    shape = {11, 500, 500, 500};
    dtype = mgard_x::data_type::Double;
    num_elements = 11 * 500 * 500 * 500;
  } else if (strstr(dataset_file, "miranda") != NULL) {
    dataset_name = "Miranda";
    shape = {3072, 3072, 3072};
    dtype = mgard_x::data_type::Float;
    num_elements = 3072 * 3072 * 3072;
  } else {
    std::cerr << "Unknown dataset: " << dataset_file << std::endl;
    fclose(file);
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

  while (iteration < MAX_ITERATIONS && !confidence_interval_reached) {
    void *compressed_data = nullptr;
    size_t compressed_size = 0;

    // Compression
    assert(PAPI_start(EventSet) == PAPI_OK);
    double start_time = get_time();

    mgard_x::compress_status_type compress_status = mgard_x::compress(
        3, dtype, shape, relative_error_bound, 0, // s = 0 for L2 norm
        mgard_x::error_bound_type::REL, data, compressed_data, compressed_size,
        config);

    double end_time = get_time();
    assert(PAPI_stop(EventSet, values) == PAPI_OK);

    if (compress_status != mgard_x::compress_status_type::Success) {
      std::cerr << "Compression failed" << std::endl;
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
        compressed_data, compressed_size, decompressed_data, config);

    end_time = get_time();
    assert(PAPI_stop(EventSet, values) == PAPI_OK);

    if (decompress_status != mgard_x::compress_status_type::Success) {
      std::cerr << "Decompression failed" << std::endl;
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
    double data_min = std::numeric_limits<float>::max();
    double data_max = std::numeric_limits<float>::lowest();
    double sum = 0.0, sum_squared = 0.0;

    for (size_t i = 0; i < num_elements; i++) {
      float orig = ((float *)data)[i];
      float decomp = ((float *)decompressed_data)[i];
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
      data_min = std::min(data_min, (double)orig);
      data_max = std::max(data_max, (double)orig);
      sum += orig;
      sum_squared += orig * orig;
    }

    double avg_difference = sum_diff / num_elements;
    double avg_error = std::sqrt(mse / num_elements);
    mse /= num_elements;
    double psnr = 20 * std::log10((data_max - data_min) / std::sqrt(mse));
    double nrmse = std::sqrt(mse) / (data_max - data_min);
    double compression_ratio =
        (double)(num_elements * sizeof(float)) / compressed_size;
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
              "%e,%e,%e,%e,%e,%e,%zu,%f,%zu,%zu,%f,%f,%f,%f\n",
              dataset_file, relative_error_bound, iteration,
              num_elements / compression_times.back(),   // compression_rate
              num_elements / decompression_times.back(), // decompression_rate
              avg_difference, avg_error,
              max_error - min_error,         // diff_range
              max_rel_error - min_rel_error, // error_range
              max_error, max_pw_rel_error, max_rel_error, min_error,
              min_pw_rel_error, min_rel_error, mse, num_elements, psnr, nrmse,
              data_max, value_mean, data_min, value_range, value_std, bit_rate,
              compressed_size, compression_ratio, num_elements * sizeof(float),
              compressed_size, compression_times.back(),
              decompression_times.back(), compression_energy.back(),
              decompression_energy.back());
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
