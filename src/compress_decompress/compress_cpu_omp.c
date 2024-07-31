#include <libpressio.h>
#include <libpressio_ext/io/posix.h>
#include <math.h>
#include <omp.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAX_ITERATIONS 50
#define CONFIDENCE_LEVEL 1.96

double calculate_mean(uint32_t *data, int n) {
  double sum = 0.0;
  for (int i = 0; i < n; i++) {
    sum += (double)data[i];
  }
  return sum / n;
}

double calculate_std_dev(uint32_t *data, int n, double mean) {
  double sum_squared_diff = 0.0;
  for (int i = 0; i < n; i++) {
    double diff = (double)data[i] - mean;
    sum_squared_diff += diff * diff;
  }
  return sqrt(sum_squared_diff / (n - 1));
}

bool within_confidence_interval(uint32_t *data, int n) {
  if (n < 2)
    return false;
  double mean = calculate_mean(data, n);
  double std_dev = calculate_std_dev(data, n, mean);
  double margin_of_error = CONFIDENCE_LEVEL * (std_dev / sqrt(n));
  double lower_bound = mean - margin_of_error;
  double upper_bound = mean + margin_of_error;

  for (int i = 0; i < n; i++) {
    if ((double)data[i] < lower_bound || (double)data[i] > upper_bound) {
      return false;
    }
  }
  return true;
}

int main(int argc, char *argv[]) {
  if (argc != 3) {
    fprintf(stderr, "Usage: %s <compressor> <dataset_file>\n", argv[0]);
    return 1;
  }

  const char *compressor_id = argv[1];
  const char *dataset_file = argv[2];
  const char *datadir = "/lcrc/project/ECP-EZ/sdrbench/";

  int num_threads;
#pragma omp parallel
  {
#pragma omp master
    {
      num_threads = omp_get_num_threads();
      printf("Running with %d OpenMP threads\n", num_threads);
    }
  }

  // get the compressor
  struct pressio *library = pressio_instance();
  struct pressio_compressor *compressor =
      pressio_get_compressor(library, compressor_id);
  if (compressor == NULL) {
    fprintf(stderr, "Failed to get compressor %s: %s\n", compressor_id,
            pressio_error_msg(library));
    pressio_release(library);
    return 1;
  }

  // configure metrics for the compressor
  struct pressio_options *metrics_options = pressio_options_new();
  const char *metrics_ids[] = {"time", "size", "error_stat"};
  size_t n_metrics_ids = sizeof(metrics_ids) / sizeof(metrics_ids[0]);
  pressio_options_set_string(metrics_options, "pressio:metric", "composite");
  pressio_options_set_strings(metrics_options, "composite:plugins",
                              n_metrics_ids, metrics_ids);

  // Set OpenMP options for the compressor
  struct pressio_options *options = pressio_compressor_get_options(compressor);
  if (strstr(compressor_id, "sz3") != NULL) {
    pressio_options_set_bool(options, "sz3:openmp", true);
  } else if (strstr(compressor_id, "mgard") != NULL) {
    pressio_options_set_bool(options, "mgard:openmp_enabled", true);
    pressio_options_set_string(options, "mgard:dev_type_str", "openmp");
    pressio_options_set_uinteger(options, "mgard:nthreads", num_threads);
  } else if (strstr(compressor_id, "zfp") != NULL) {
    pressio_options_set_string(options, "zfp:execution_name", "omp");
    pressio_options_set_uinteger(options, "zfp:omp_threads", num_threads);
  }
  pressio_compressor_set_options(compressor, options);

  // Read the dataset
  size_t ndims;
  struct pressio_data *metadata, *input_data;
  if (strstr(dataset_file, "nyx") != NULL) {
    size_t dims[] = {512, 512, 512};
    ndims = sizeof(dims) / sizeof(dims[0]);
    metadata = pressio_data_new_empty(pressio_float_dtype, ndims, dims);
  } else if (strstr(dataset_file, "hacc") != NULL) {
    size_t dims[] = {1073726487};
    ndims = sizeof(dims) / sizeof(dims[0]);
    metadata = pressio_data_new_empty(pressio_float_dtype, ndims, dims);
  } else if (strstr(dataset_file, "s3d") != NULL) {
    size_t dims[] = {11, 500, 500, 500};
    ndims = sizeof(dims) / sizeof(dims[0]);
    metadata = pressio_data_new_empty(pressio_double_dtype, ndims, dims);
  } else if (strstr(dataset_file, "miranda") != NULL) {
    size_t dims[] = {3072, 3072, 3072};
    ndims = sizeof(dims) / sizeof(dims[0]);
    metadata = pressio_data_new_empty(pressio_float_dtype, ndims, dims);
  } else {
    fprintf(stderr, "Unknown dataset %s\n", dataset_file);
    pressio_compressor_release(compressor);
    pressio_release(library);
    return 1;
  }

  char full_path[1024];
  snprintf(full_path, sizeof(full_path), "%s%s", datadir, dataset_file);
  input_data = pressio_io_data_path_read(metadata, full_path);
  if (input_data == NULL) {
    fprintf(stderr, "Failed to read dataset %s\n", dataset_file);
    pressio_compressor_release(compressor);
    pressio_release(library);
    return 1;
  }

  // Create compressed and output data structures
  struct pressio_data *compressed =
      pressio_data_new_empty(pressio_byte_dtype, 0, NULL);
  struct pressio_data *output = pressio_data_new_clone(input_data);

  double bounds[] = {1e-6, 5e-6, 1e-5, 5e-5, 1e-4, 5e-4,
                     1e-3, 5e-3, 1e-2, 5e-2, 1e-1};
  size_t n_bounds = sizeof(bounds) / sizeof(bounds[0]);

  printf("Compressor: %s\n", compressor_id);
  printf("Dataset: %s\n", dataset_file);

  for (size_t i = 0; i < n_bounds; ++i) {
    // configure the compressor error bound
    struct pressio_options *bound_options = pressio_options_new();
    pressio_options_set_double(bound_options, "pressio:rel", bounds[i]);
    pressio_compressor_set_options(compressor, bound_options);
    pressio_options_free(bound_options);

    uint32_t compression_times[MAX_ITERATIONS];
    uint32_t decompression_times[MAX_ITERATIONS];
    int iteration = 0;
    bool confidence_interval_reached = false;

    while (iteration < MAX_ITERATIONS && !confidence_interval_reached) {
      // run the compression and decompression
      if (pressio_compressor_compress(compressor, input_data, compressed)) {
        fprintf(stderr, "%s\n", pressio_compressor_error_msg(compressor));
        break;
      }
      if (pressio_compressor_decompress(compressor, compressed, output)) {
        fprintf(stderr, "%s\n", pressio_compressor_error_msg(compressor));
        break;
      }

      // get metrics results and write to CSV
      struct pressio_options *metrics_results =
          pressio_compressor_get_metrics_results(compressor);

      // ... (keep the metric extraction and CSV writing code as it is)

      pressio_options_free(metrics_results);

      iteration++;

      // check if we've reached the confidence interval
      if (iteration >= 5) {
        confidence_interval_reached =
            within_confidence_interval(compression_times, iteration) &&
            within_confidence_interval(decompression_times, iteration);
      }
    }
  }

  // Clean up
  pressio_data_free(metadata);
  pressio_data_free(input_data);
  pressio_data_free(compressed);
  pressio_data_free(output);
  pressio_options_free(metrics_options);
  pressio_compressor_release(compressor);
  pressio_release(library);

  return 0;
}