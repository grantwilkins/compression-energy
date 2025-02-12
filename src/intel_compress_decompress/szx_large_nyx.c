#include "szx.h"
#include "szx_rw.h"
#include <assert.h>
#include <math.h>
#include <papi.h>
#include <sched.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAX_ITERATIONS 10
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
  if (argc != 3) {
    printf("Usage: %s <dataset_file> <relative_error_bound>\n", argv[0]);
    return 1;
  }

  const char *dataset_file = argv[1];
  double relative_error_bound = atof(argv[2]);
  const char *datadir = "/work2/10191/gfw/stampede3/";
  const char *cluster_name = getenv("CLUSTER_NAME");

  double compression_rate, decompression_rate, avg_difference, avg_error,
      diff_range, error_range, max_pw_rel_error, max_rel_error, min_error,
      min_pw_rel_error, min_rel_error, value_mean, value_std, bit_rate;
  uint64_t compressed_size, decompressed_size, uncompressed_size;
  double max_error, mse, psnr, nrmse;
  double value_range, value_min, value_max;
  double sum_squared_diff, sum_diff, sum, sum_squared;

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

  int status = 0;
  size_t nbEle;
  void *data = NULL;
  size_t dims[5] = {0};
  int data_type_szx;

  char full_path[1024];
  snprintf(full_path, sizeof(full_path), "%s%s", datadir, dataset_file);

  if (strstr(dataset_file, "nyx") != NULL) {
    if(strstr(dataset_file, "_2") != NULL) {
      printf("Inflated by 2\n");
      dims[0] = dims[1] = dims[2] = 1024;
    } else if (strstr(dataset_file, "_3") != NULL) {
      printf("Inflated by 3\n");
      dims[0] = dims[1] = dims[2] = 1536;
    } else if (strstr(dataset_file, "_4") != NULL) {
      printf("Inflated by 4\n");
      dims[0] = dims[1] = dims[2] = 2048;
    } else if (strstr(dataset_file, "_5") != NULL) {
      printf("Inflated by 5\n");
      dims[0] = dims[1] = dims[2] = 2560;
    } else {
      printf("Original size\n");
      dims[0] = dims[1] = dims[2] = 512;
    }
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

    float *float_data = (float *)data;
    data_min = data_max = float_data[0];
    for (size_t i = 1; i < nbEle; i++) {
        if (float_data[i] < data_min)
        data_min = float_data[i];
        if (float_data[i] > data_max)
        data_max = float_data[i];
    }
  
  data_range = data_max - data_min;
  double absolute_error_bound = relative_error_bound * data_range;

  printf("Dataset: %s\n", dataset_file);
  printf("REL Error bound: %e\n", relative_error_bound);
  printf("ABS Error bound: %e\n", absolute_error_bound);

  double compression_times[MAX_ITERATIONS];
  double decompression_times[MAX_ITERATIONS];
  double compression_energy[MAX_ITERATIONS];
  double decompression_energy[MAX_ITERATIONS];
  int iteration = 0;
  bool confidence_interval_reached = false;

  while (iteration < MAX_ITERATIONS && !confidence_interval_reached) {
    size_t outSize;
    max_error = 0.0;
    min_error = INFINITY;
    max_rel_error = 0.0;
    min_rel_error = INFINITY;
    max_pw_rel_error = 0.0;
    min_pw_rel_error = INFINITY;
    sum_squared_diff = 0.0;
    sum_diff = 0.0;
    sum = 0.0;
    sum_squared = 0.0;
    value_min = INFINITY;
    value_max = -INFINITY;

    // Compression
    assert(PAPI_start(EventSet) == PAPI_OK);
    double start_time = get_time();

    unsigned char *compressed_data = SZx_fast_compress_args(
        SZx_NO_BLOCK_FAST_CMPR, data_type_szx, data, &outSize, REL,
        absolute_error_bound, relative_error_bound, 0, 0, dims[4], dims[3],
        dims[2], dims[1], dims[0]);

    double end_time = get_time();
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
    compression_energy[iteration] =
        cpu_energy_compression + dram_energy_compression;
    compression_times[iteration] = (end_time - start_time);

    // Decompression
    assert(PAPI_start(EventSet) == PAPI_OK);
    start_time = get_time();

    void *decompressed_data = SZx_fast_decompress(
        SZx_NO_BLOCK_FAST_CMPR, data_type_szx, compressed_data, outSize,
        dims[4], dims[3], dims[2], dims[1], dims[0]);

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
    decompression_energy[iteration] =
        cpu_energy_decompression + dram_energy_decompression;
    decompression_times[iteration] = (end_time - start_time);


    float *float_data = (float *)data;
    float *float_decompressed = (float *)decompressed_data;

    for (size_t i = 0; i < nbEle; i++) {
        double orig = float_data[i];
        double decomp = float_decompressed[i];
        double diff = decomp - orig;
        double abs_diff = fabs(diff);
        double rel_diff = abs_diff / (fabs(orig) + 1e-6);
        double pw_rel_diff = abs_diff / (fabs(orig) + 1e-6);

        max_error = fmax(max_error, abs_diff);
        min_error = fmin(min_error, abs_diff);
        max_rel_error = fmax(max_rel_error, rel_diff);
        min_rel_error = fmin(min_rel_error, rel_diff);
        max_pw_rel_error = fmax(max_pw_rel_error, pw_rel_diff);
        min_pw_rel_error = fmin(min_pw_rel_error, pw_rel_diff);

        sum_squared_diff += diff * diff;
        sum_diff += diff;
        sum += orig;
        sum_squared += orig * orig;

        value_min = fmin(value_min, orig);
        value_max = fmax(value_max, orig);
    }

    value_range = value_max - value_min;
    mse = sum_squared_diff / nbEle;
    avg_difference = sum_diff / nbEle;
    avg_error = sqrt(mse);
    psnr = 20 * log10(value_range) - 10 * log10(mse);
    nrmse = sqrt(mse) / value_range;

    compression_rate = (double)nbEle / compression_times[iteration];
    decompression_rate = (double)nbEle / decompression_times[iteration];
    diff_range = max_error - min_error;
    error_range = max_rel_error - min_rel_error;

    value_mean = sum / nbEle;
    value_std = sqrt(sum_squared / nbEle - value_mean * value_mean);

    compressed_size = outSize;
    decompressed_size =
        nbEle * (data_type_szx == SZx_FLOAT ? sizeof(float) : sizeof(double));
    bit_rate = (double)compressed_size * 8 / nbEle;

    double compression_ratio = (double)decompressed_size / compressed_size;
    printf("ENERGY: %e\n", decompression_energy[iteration]);
    printf("TIME: %e\n", decompression_times[iteration]);
    // Write metrics to CSV file
    FILE *csv_file = fopen("compression_metrics_szx_serial.csv", "a");
    if (csv_file == NULL) {
      fprintf(stderr, "Error opening CSV file\n");
    } else {
      fprintf(
          csv_file,
          "SZx,%s,%e,%e,%d,%f,%f,%e,%e,%e,%e,%e,%zu,%e,%e,%e,"
          "%zu,%f,%zu,%f,%f,%f,%f\n",
          dataset_file, relative_error_bound, absolute_error_bound, iteration,
          compression_rate, decompression_rate, avg_difference, avg_error,
          diff_range, error_range,mse, nbEle, psnr, nrmse, bit_rate,
          compressed_size, compression_ratio, decompressed_size,
          compression_times[iteration],
          decompression_times[iteration], compression_energy[iteration],
          decompression_energy[iteration]);
      fclose(csv_file);
    }

    free(compressed_data);
    free(decompressed_data);

    iteration++;

    // Check if we've reached the confidence interval
    if (iteration >= 5) {
      confidence_interval_reached =
          within_confidence_interval(compression_times, iteration) &&
          within_confidence_interval(decompression_times, iteration);
    }
  }

  // Clean up
  free(data);
  PAPI_cleanup_eventset(EventSet);
  PAPI_destroy_eventset(&EventSet);
  PAPI_shutdown();

  return 0;
}
