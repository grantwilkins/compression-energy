#include <libpressio.h>
#include <libpressio_ext/io/posix.h>
#include <math.h>
#include <sched.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAX_ITERATIONS 10
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
  if (argc != 4) {
    fprintf(stderr, "Usage: %s <compressor> <dataset_file> <error_bound>\n",
            argv[0]);
    return 1;
  }

  const char *compressor_id = argv[1];
  const char *dataset_file = argv[2];
  double error_bound = atof(argv[3]);
  const char *datadir = "/lcrc/project/ECP-EZ/sdrbench/";

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
  } else if (strstr(dataset_file, "cesm") != NULL) {
    size_t dims[] = {26, 1800, 3600};
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

  // verify that options passed exist
  if (pressio_compressor_check_options(compressor, metrics_options)) {
    fprintf(stderr, "%s\n", pressio_compressor_error_msg(compressor));
    pressio_compressor_release(compressor);
    pressio_release(library);
    return 1;
  }

  // set the new options for the compressor
  if (pressio_compressor_set_options(compressor, metrics_options)) {
    fprintf(stderr, "%s\n", pressio_compressor_error_msg(compressor));
    pressio_compressor_release(compressor);
    pressio_release(library);
    return 1;
  }

  printf("Compressor: %s\n", compressor_id);
  printf("Dataset: %s\n", dataset_file);
  printf("Error bound: %e\n", error_bound);
  struct pressio_data *compressed =
      pressio_data_new_empty(pressio_byte_dtype, 0, NULL);
  struct pressio_data *output = pressio_data_new_clone(input_data);

  // configure the compressor error bound
  struct pressio_options *options = pressio_options_new();
  pressio_options_set_double(options, "pressio:abs", error_bound);

  // verify that options passed exist
  if (pressio_compressor_check_options(compressor, options)) {
    fprintf(stderr, "%s\n", pressio_compressor_error_msg(compressor));
    pressio_options_free(options);
    return 1;
  }

  // set the new options for the compressor
  if (pressio_compressor_set_options(compressor, options)) {
    fprintf(stderr, "%s\n", pressio_compressor_error_msg(compressor));
    pressio_options_free(options);
    return 1;
  }
  pressio_options_free(options);

  uint32_t compression_times[MAX_ITERATIONS];
  uint32_t decompression_times[MAX_ITERATIONS];
  int iteration = 0;
  bool confidence_interval_reached = false;

  while (iteration < MAX_ITERATIONS && !confidence_interval_reached) {
    // run the compression and decompression
    int cpu = sched_getcpu();
    if (pressio_compressor_compress(compressor, input_data, compressed)) {
      fprintf(stderr, "%s\n", pressio_compressor_error_msg(compressor));
      break;
    }
    if (pressio_compressor_decompress(compressor, compressed, output)) {
      fprintf(stderr, "%s\n", pressio_compressor_error_msg(compressor));
      break;
    }

    // get metrics results
    struct pressio_options *metrics_results =
        pressio_compressor_get_metrics_results(compressor);

    // Extract metrics and store them in variables
    double compression_rate, decompression_rate, avg_difference, avg_error,
        diff_range, error_range;
    double max_error, max_pw_rel_error, max_rel_error, min_error,
        min_pw_rel_error, min_rel_error;
    double mse, psnr, rmse, value_max, value_mean, value_min, value_range,
        value_std, bit_rate, nrmse;
    double compression_ratio;
    uint64_t n, compressed_size, decompressed_size, uncompressed_size;
    uint32_t compression_time, decompression_time;

    pressio_options_get_double(metrics_results, "composite:compression_rate",
                               &compression_rate);
    pressio_options_get_double(metrics_results, "composite:decompression_rate",
                               &decompression_rate);
    pressio_options_get_double(metrics_results, "error_stat:average_difference",
                               &avg_difference);
    pressio_options_get_double(metrics_results, "error_stat:average_error",
                               &avg_error);
    pressio_options_get_double(metrics_results, "error_stat:difference_range",
                               &diff_range);
    pressio_options_get_double(metrics_results, "error_stat:error_range",
                               &error_range);
    pressio_options_get_double(metrics_results, "error_stat:max_error",
                               &max_error);
    pressio_options_get_double(metrics_results, "error_stat:max_pw_rel_error",
                               &max_pw_rel_error);
    pressio_options_get_double(metrics_results, "error_stat:max_rel_error",
                               &max_rel_error);
    pressio_options_get_double(metrics_results, "error_stat:min_error",
                               &min_error);
    pressio_options_get_double(metrics_results, "error_stat:min_pw_rel_error",
                               &min_pw_rel_error);
    pressio_options_get_double(metrics_results, "error_stat:min_rel_error",
                               &min_rel_error);
    pressio_options_get_double(metrics_results, "error_stat:mse", &mse);
    pressio_options_get_uinteger64(metrics_results, "error_stat:n", &n);
    pressio_options_get_double(metrics_results, "error_stat:psnr", &psnr);
    pressio_options_get_double(metrics_results, "error_stat:rmse", &rmse);
    pressio_options_get_double(metrics_results, "error_stat:value_max",
                               &value_max);
    pressio_options_get_double(metrics_results, "error_stat:value_mean",
                               &value_mean);
    pressio_options_get_double(metrics_results, "error_stat:value_min",
                               &value_min);
    pressio_options_get_double(metrics_results, "error_stat:value_range",
                               &value_range);
    pressio_options_get_double(metrics_results, "error_stat:value_std",
                               &value_std);
    pressio_options_get_double(metrics_results, "size:bit_rate", &bit_rate);
    pressio_options_get_uinteger64(metrics_results, "size:compressed_size",
                                   &compressed_size);
    pressio_options_get_double(metrics_results, "size:compression_ratio",
                               &compression_ratio);
    pressio_options_get_uinteger64(metrics_results, "size:decompressed_size",
                                   &decompressed_size);
    pressio_options_get_uinteger64(metrics_results, "size:uncompressed_size",
                                   &uncompressed_size);
    pressio_options_get_uinteger(metrics_results, "time:compress",
                                 &compression_time);
    pressio_options_get_uinteger(metrics_results, "time:decompress",
                                 &decompression_time);

    nrmse = rmse / value_range;

    // Write metrics to CSV file
    FILE *csv_file = fopen("compression_metrics.csv", "a");
    if (csv_file == NULL) {
      fprintf(stderr, "Error opening CSV file\n");
    } else {
      fprintf(csv_file,
              "%s,%s,%e,%d,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%lu,%e,%e,%e,"
              "%e,%e,%e,%e,%e,%e,%lu,%e,%lu,%lu,%u,%u,%d\n",
              compressor_id, dataset_file, error_bound, iteration,
              compression_rate, decompression_rate, avg_difference, avg_error,
              diff_range, error_range, max_error, max_pw_rel_error,
              max_rel_error, min_error, min_pw_rel_error, min_rel_error, mse, n,
              psnr, rmse, nrmse, value_max, value_mean, value_min, value_range,
              value_std, bit_rate, compressed_size, compression_ratio,
              decompressed_size, uncompressed_size, compression_time,
              decompression_time, cpu);
      fclose(csv_file);
    }
    compression_times[iteration] = compression_time;
    decompression_times[iteration] = decompression_time;

    pressio_options_free(metrics_results);

    iteration++;

    // check if we've reached the confidence interval
    if (iteration >= 5) {
      confidence_interval_reached =
          within_confidence_interval(compression_times, iteration) &&
          within_confidence_interval(decompression_times, iteration);
    }
  }

  pressio_data_free(metadata);
  pressio_data_free(input_data);
  pressio_data_free(compressed);
  pressio_data_free(output);

  pressio_options_free(metrics_options);
  pressio_compressor_release(compressor);
  pressio_release(library);
  return 0;
}