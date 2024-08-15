#include <cuda_runtime.h>
#include <libpressio.h>
#include <libpressio_ext/io/posix.h>
#include <math.h>
#include <nvml.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zfp.h>

#define MAX_ITERATIONS 10
#define CONFIDENCE_LEVEL 1.96

#define CHECK_NVML(call)                                                       \
  do {                                                                         \
    nvmlReturn_t result = call;                                                \
    if (result != NVML_SUCCESS) {                                              \
      fprintf(stderr, "NVML Error: %s\n", nvmlErrorString(result));            \
      exit(1);                                                                 \
    }                                                                          \
  } while (0)

void lp_cuda_free(void *ptr, void *metadata) { cudaFree(ptr); }

void lp_check_cuda(cudaError_t err) {
  if (err != cudaSuccess) {
    fprintf(stderr, "%s\n", cudaGetErrorString(err));
    exit(1);
  }
}

double calculate_mean(double *data, int n) {
  double sum = 0.0;
  for (int i = 0; i < n; i++) {
    sum += data[i];
  }
  return sum / n;
}

double calculate_std_dev(double *data, int n, double mean) {
  double sum_squared_diff = 0.0;
  for (int i = 0; i < n; i++) {
    double diff = data[i] - mean;
    sum_squared_diff += diff * diff;
  }
  return sqrt(sum_squared_diff / (n - 1));
}

bool within_confidence_interval(double *data, int n) {
  if (n < 2)
    return false;
  double mean = calculate_mean(data, n);
  double std_dev = calculate_std_dev(data, n, mean);
  double margin_of_error = CONFIDENCE_LEVEL * (std_dev / sqrt(n));
  double lower_bound = mean - margin_of_error;
  double upper_bound = mean + margin_of_error;

  for (int i = 0; i < n; i++) {
    if (data[i] < lower_bound || data[i] > upper_bound) {
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
  double relative_error_bound = atof(argv[3]);
  nvmlDevice_t nvml_device;
  int cuda_device;
  unsigned long long start_gpu_energy = 0, end_gpu_energy = 0;
  double compression_rate = 0.0, decompression_rate = 0.0, avg_difference = 0.0,
         avg_error = 0.0, diff_range = 0.0, error_range = 0.0;
  double max_error = 0.0, max_pw_rel_error = 0.0, max_rel_error = 0.0,
         min_error = 0.0, min_pw_rel_error = 0.0, min_rel_error = 0.0;
  double mse = 0.0, psnr = 0.0, rmse = 0.0, value_max = 0.0, value_mean = 0.0,
         value_min = 0.0, value_range = 0.0, value_std = 0.0, bit_rate = 0.0,
         nrmse = 0.0;
  double compression_ratio = 0.0;
  uint64_t n = 0, compressed_size = 0, decompressed_size = 0,
           uncompressed_size = 0;
  double compression_time = 0.0, decompression_time = 0.0;
  int cpu = -1;

  // Initialize NVML
  CHECK_NVML(nvmlInit());

  // Get the handle for the first GPU (index 0)
  CHECK_NVML(nvmlDeviceGetHandleByIndex(0, &nvml_device));

  // Initialize the library
  struct pressio *library = pressio_instance();
  struct pressio_compressor *comp =
      pressio_get_compressor(library, compressor_id);

  if (comp == NULL) {
    fprintf(stderr, "Failed to get compressor %s\n", compressor_id);
    pressio_release(library);
    return 1;
  }

  // Read the dataset
  size_t dims[] = {500, 500, 100}; // Adjust based on your dataset
  size_t ndims = sizeof(dims) / sizeof(dims[0]);
  size_t buf_size = sizeof(float);
  for (size_t i = 0; i < ndims; ++i) {
    buf_size *= dims[i];
  }

  float *d_input;
  cudaMallocManaged((void **)&d_input, buf_size, 0);
  struct pressio_data *metadata = pressio_data_new_move(
      pressio_float_dtype, d_input, ndims, dims, lp_cuda_free, NULL);
  struct pressio_data *input_data_shared =
      pressio_io_data_path_read(metadata, dataset_file);

  float *d_output;
  cudaMallocManaged((void **)&d_output, buf_size, 0);
  struct pressio_data *output_shared = pressio_data_new_move(
      pressio_float_dtype, d_output, ndims, dims, lp_cuda_free, NULL);

  double data_min, data_max, data_range;
  size_t num_elements = 1;
  for (int t = 0; t < ndims; t++)
    num_elements *= dims[t];
  void *data_ptr = (void *)d_input;

  if (strstr(dataset_file, "s3d") == NULL) {
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

  // Allocate compressed buffer
  size_t comp_buf_size;
  if (strcmp(compressor_id, "zfp") == 0) {
    double rate = 4.0; // Adjust as needed
    zfp_stream *zfp = zfp_stream_open(NULL);
    zfp_stream_set_rate(zfp, rate, zfp_type_float, 3, 0);
    zfp_field *field =
        zfp_field_3d(NULL, zfp_type_float, dims[0], dims[1], dims[2]);
    comp_buf_size = zfp_stream_maximum_size(zfp, field);
    zfp_stream_close(zfp);
    zfp_field_free(field);
  } else {
    fprintf(stderr, "Unsupported compressor: %s\n", compressor_id);
    return 1;
  }

  float *d_compressed;
  cudaMallocManaged((void **)&d_compressed, comp_buf_size, 0);
  struct pressio_data *comp_shared = pressio_data_new_move(
      pressio_byte_dtype, d_compressed, 1, &comp_buf_size, lp_cuda_free, NULL);

  // Put the memory on the device
  cudaGetDevice(&cuda_device);
  cudaMemPrefetchAsync(d_input, buf_size, cuda_device, NULL);

  // Set compressor options
  struct pressio_options *options = pressio_options_new();
  if (strcmp(compressor_id, "zfp") == 0) {
    pressio_options_set_string(options, "zfp:execution_name", "cuda");
  } else if (strcmp(compressor_id, "mgard") == 0) {
    pressio_options_set_string(options, "mgard:execution_mode", "cuda");
  }

  // Set metrics
  pressio_options_set_string(options, "pressio:metric", "composite");
  const char *plugins[] = {"time", "size", "error_stat"};
  size_t n_plugins = sizeof(plugins) / sizeof(plugins[0]);
  pressio_options_set_strings(options, "composite:plugins", n_plugins, plugins);
  pressio_options_set_double(options, "pressio:abs", absolute_error_bound);

  pressio_compressor_set_options(comp, options);
  double compression_times[MAX_ITERATIONS];
  double decompression_times[MAX_ITERATIONS];
  int iteration = 0;
  bool confidence_interval_reached = false;

  while (iteration < MAX_ITERATIONS && !confidence_interval_reached) {
    // Compression
    CHECK_NVML(
        nvmlDeviceGetTotalEnergyConsumption(nvml_device, &start_gpu_energy));
    if (pressio_compressor_compress(comp, input_data_shared, comp_shared)) {
      fprintf(stderr, "Compression error: %s\n",
              pressio_compressor_error_msg(comp));
      exit(1);
    }
    CHECK_NVML(
        nvmlDeviceGetTotalEnergyConsumption(nvml_device, &end_gpu_energy));

    unsigned long long compression_energy_consumed =
        end_gpu_energy - start_gpu_energy;

    // Decompression
    CHECK_NVML(
        nvmlDeviceGetTotalEnergyConsumption(nvml_device, &start_gpu_energy));
    if (pressio_compressor_decompress(comp, comp_shared, output_shared)) {
      fprintf(stderr, "Decompression error: %s\n",
              pressio_compressor_error_msg(comp));
      exit(1);
    }
    CHECK_NVML(
        nvmlDeviceGetTotalEnergyConsumption(nvml_device, &end_gpu_energy));

    unsigned long long decompression_energy_consumed =
        end_gpu_energy - start_gpu_energy;

    // get metrics results and write to CSV
    struct pressio_options *metrics_results =
        pressio_compressor_get_metrics_results(comp);
    cpu = sched_getcpu();
    // Extract metrics and store them in variables

    // you should really check the return code for these to make sure you get
    // them all pressio_options_key_set is returned on success, you could use a
    // macro or a function to check if this successes.
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

    // Write metrics to CSV file
    // Write metrics to CSV file
    FILE *csv_file = fopen("compression_metrics_omp.csv", "a");
    if (csv_file == NULL) {
      fprintf(stderr, "Error opening CSV file\n");
    } else {
      fprintf(
          csv_file,
          "%s,%s,%e,%e,%d,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%lu,%e,%e,%e,"
          "%e,%e,%e,%e,%e,%e,%lu,%e,%lu,%lu,%f,%f,%d,%d\n",
          compressor_id, dataset_file, relative_error_bound,
          absolute_error_bound, iteration, compression_rate, decompression_rate,
          avg_difference, avg_error, diff_range, error_range, max_error,
          max_pw_rel_error, max_rel_error, min_error, min_pw_rel_error,
          min_rel_error, mse, n, psnr, rmse, nrmse, value_max, value_mean,
          value_min, value_range, value_std, bit_rate, compressed_size,
          compression_ratio, decompressed_size, uncompressed_size,
          compression_time, decompression_time, cpu);
      fclose(csv_file);
    }

    pressio_options_free(metrics_results);

    iteration++;

    // check if we've reached the confidence interval
    if (iteration >= 5) {
      confidence_interval_reached =
          within_confidence_interval(compression_times, iteration) &&
          within_confidence_interval(decompression_times, iteration);
    }
  }
  // Clean up
  pressio_data_free(comp_shared);
  pressio_data_free(output_shared);
  pressio_data_free(input_data_shared);
  pressio_options_free(options);
  pressio_compressor_release(comp);
  pressio_release(library);

  return 0;
}
