#include <assert.h>
#include <libpressio.h>
#include <libpressio_ext/io/posix.h>
#include <math.h>
#include <omp.h>
#include <papi.h>
#include <sched.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAX_ITERATIONS 25
#define CONFIDENCE_LEVEL 1.96
#define MAX_powercap_EVENTS 64

double get_time() {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return ts.tv_sec + ts.tv_nsec / 1e9;
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
  const char *datadir = "/ocean/projects/cis240100p/gwilkins/";
  const char *cluster_name = getenv("CLUSTER_NAME");

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

  int num_threads =
      omp_get_max_threads(); // I think this does the same as the following
#pragma omp parallel
  {
#pragma omp master
    {
      num_threads = omp_get_num_threads();
      printf("Running with %d OpenMP threads\n", num_threads);
    }
  }

  // PAPI initialization
  int EventSet = PAPI_NULL;
  long long *values;
  int num_events = 0;
  int code;
  char event_names[MAX_powercap_EVENTS][PAPI_MAX_STR_LEN];
  char event_descrs[MAX_powercap_EVENTS][PAPI_MAX_STR_LEN];
  char units[MAX_powercap_EVENTS][PAPI_MIN_STR_LEN];
  int data_type[MAX_powercap_EVENTS];
  int r, i;

  assert(PAPI_library_init(PAPI_VER_CURRENT) == PAPI_VER_CURRENT);

  // Find powercap component
  int numcmp = PAPI_num_components();
  int cid, powercap_cid = -1;
  const PAPI_component_info_t *cmpinfo = NULL;
  for (cid = 0; cid < numcmp; cid++) {
    cmpinfo = PAPI_get_component_info(cid);
    assert(cmpinfo != NULL);
    if (strstr(cmpinfo->name, "powercap")) {
      powercap_cid = cid;
      break;
    }
  }
  assert(cid != numcmp);

  // Create EventSet
  assert(PAPI_create_eventset(&EventSet) == PAPI_OK);

  // Find all events
  code = PAPI_NATIVE_MASK;
  r = PAPI_enum_cmp_event(&code, PAPI_ENUM_FIRST, powercap_cid);
  while (r == PAPI_OK) {
    PAPI_event_info_t evinfo;
    assert(PAPI_event_code_to_name(code, event_names[num_events]) == PAPI_OK);
    assert(PAPI_get_event_info(code, &evinfo) == PAPI_OK);
    strncpy(event_descrs[num_events], evinfo.long_descr,
            sizeof(event_descrs[0]) - 1);
    strncpy(units[num_events], evinfo.units, sizeof(units[0]) - 1);
    units[num_events][sizeof(units[0]) - 1] = '\0';
    data_type[num_events] = evinfo.data_type;
    if (PAPI_add_event(EventSet, code) != PAPI_OK)
      break;
    num_events++;
    r = PAPI_enum_cmp_event(&code, PAPI_ENUM_EVENTS, powercap_cid);
  }

  values = (long long *)calloc(num_events, sizeof(long long));

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

  // Create compressed and output data structures
  struct pressio_data *compressed =
      pressio_data_new_empty(pressio_byte_dtype, 0, NULL);
  struct pressio_data *output = pressio_data_new_clone(input_data);

  double data_min, data_max, data_range;
  size_t num_elements = pressio_data_num_elements(input_data);
  void *data_ptr = pressio_data_ptr(input_data, NULL);

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

  printf("Compressor: %s\n", compressor_id);
  printf("Dataset: %s\n", dataset_file);
  printf("REL Error bound: %e\n", relative_error_bound);

  // configure the compressor error bound
  struct pressio_options *bound_options = pressio_options_new();
  pressio_options_set_double(bound_options, "pressio:abs",
                             absolute_error_bound);
  pressio_options_set_uinteger(bound_options, "pressio:nthreads", num_threads);

  const char *metrics_ids[] = {"size", "error_stat", "time"};
  size_t n_metrics_ids = sizeof(metrics_ids) / sizeof(metrics_ids[0]);
  pressio_options_set_string(bound_options, "pressio:metric", "composite");
  pressio_options_set_strings(bound_options, "composite:plugins", n_metrics_ids,
                              metrics_ids);

  pressio_compressor_set_options(compressor, bound_options);
  pressio_options_free(bound_options);

  double compression_times[MAX_ITERATIONS];
  double decompression_times[MAX_ITERATIONS];
  double compression_energy[MAX_ITERATIONS];
  double decompression_energy[MAX_ITERATIONS];
  int iteration = 0;
  bool confidence_interval_reached = false;

  while (iteration < MAX_ITERATIONS && !confidence_interval_reached) {
    // run the compression and decompression
    double start_time, end_time;

    assert(PAPI_start(EventSet) == PAPI_OK);
    start_time = get_time();

    if (pressio_compressor_compress(compressor, input_data, compressed)) {
      fprintf(stderr, "%s\n", pressio_compressor_error_msg(compressor));
      break;
    }

    end_time = get_time();
    assert(PAPI_stop(EventSet, values) == PAPI_OK);

    double cpu_energy_compression = 0.0, dram_energy_compression = 0.0;
    for (i = 0; i < num_events; i++) {
      if (strstr(event_names[i], "ENERGY_UJ")) {
        if (data_type[i] == PAPI_DATATYPE_UINT64) {
          if (strstr(event_names[i], "SUBZONE")) {
            dram_energy_compression += values[i] / 1.0e6;
          } else {
            cpu_energy_compression += values[i] / 1.0e6;
          }
        }
      }
    }
    compression_energy[iteration] = cpu_energy_compression;
    compression_times[iteration] = (end_time - start_time);

    // Measure decompression energy
    assert(PAPI_start(EventSet) == PAPI_OK);
    start_time = get_time();

    if (pressio_compressor_decompress(compressor, compressed, output)) {
      fprintf(stderr, "%s\n", pressio_compressor_error_msg(compressor));
      break;
    }

    end_time = get_time();
    assert(PAPI_stop(EventSet, values) == PAPI_OK);

    double cpu_energy_decompression = 0.0, dram_energy_decompression = 0.0;
    for (i = 0; i < num_events; i++) {
      if (strstr(event_names[i], "ENERGY_UJ")) {
        if (data_type[i] == PAPI_DATATYPE_UINT64) {
          if (strstr(event_names[i], "SUBZONE")) {
            dram_energy_decompression += values[i] / 1.0e6;
          } else {
            cpu_energy_decompression += values[i] / 1.0e6;
          }
        }
      }
    }
    decompression_energy[iteration] = cpu_energy_decompression;
    decompression_times[iteration] = (end_time - start_time);

    // get metrics results and write to CSV
    struct pressio_options *metrics_results =
        pressio_compressor_get_metrics_results(compressor);
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

    compression_time = compression_times[iteration];
    decompression_time = decompression_times[iteration];

    nrmse = rmse / value_range;

#ifdef DEBUG
    // helpful debugging check to print out what libpressio got as
    // metrics_results
    char *str = pressio_options_to_string(metrics_results);
    fprintf(stderr, "%s\n", str);
    free(str);
#endif

    // Write metrics to CSV file
    FILE *csv_file = fopen("compression_metrics_qoz_omp.csv", "a");
    if (csv_file == NULL) {
      fprintf(stderr, "Error opening CSV file\n");
    } else {
      fprintf(
          csv_file,
          "%s,%s,%e,%e,%d,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%lu,%e,%e,%e,"
          "%e,%e,%e,%e,%e,%e,%lu,%e,%lu,%lu,%f,%f,%d,%d,%e,%e\n",
          compressor_id, dataset_file, relative_error_bound,
          absolute_error_bound, iteration, compression_rate, decompression_rate,
          avg_difference, avg_error, diff_range, error_range, max_error,
          max_pw_rel_error, max_rel_error, min_error, min_pw_rel_error,
          min_rel_error, mse, n, psnr, rmse, nrmse, value_max, value_mean,
          value_min, value_range, value_std, bit_rate, compressed_size,
          compression_ratio, decompressed_size, uncompressed_size,
          compression_time, decompression_time, cpu, num_threads,
          cpu_energy_compression, cpu_energy_decompression);
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
  pressio_data_free(metadata);
  pressio_data_free(input_data);
  pressio_data_free(compressed);
  pressio_data_free(output);
  pressio_compressor_release(compressor);
  pressio_release(library);

  return 0;
}
