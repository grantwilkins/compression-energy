#include "SPERR_C_API.h"
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAX_ITERATIONS 25
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

size_t read_file(const char *filename, void **dst) {
  FILE *f = fopen(filename, "r");
  if (!f)
    return 0;
  fseek(f, 0, SEEK_END);
  const size_t len = ftell(f);
  fseek(f, 0, SEEK_SET);
  uint8_t *buf = malloc(len);
  size_t rtn = fread(buf, 1, len, f);
  if (rtn != len) {
    free(buf);
    fclose(f);
    return 0;
  }
  fclose(f);
  *dst = buf;
  return len;
}

int main(int argc, char *argv[]) {
  if (argc != 3) {
    printf("Usage: %s <dataset_file> <num_threads>\n", argv[0]);
    return 1;
  }

  const char *dataset_file = argv[1];
  int num_threads = atoi(argv[2]);
  omp_set_num_threads(num_threads);

  double bounds[] = {1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1};
  size_t n_bounds = sizeof(bounds) / sizeof(bounds[0]);

  void *data = NULL;
  size_t nbEle;
  size_t dims[3] = {0};
  int is_float = 1;

  if (strstr(dataset_file, "nyx") != NULL) {
    dims[0] = dims[1] = dims[2] = 512;
  } else if (strstr(dataset_file, "hacc") != NULL) {
    dims[0] = 1073726487;
    dims[1] = dims[2] = 1;
  } else if (strstr(dataset_file, "s3d") != NULL) {
    dims[0] = dims[1] = dims[2] = 500;
    is_float = 0;
  } else if (strstr(dataset_file, "miranda") != NULL) {
    dims[0] = dims[1] = dims[2] = 3072;
  } else {
    fprintf(stderr, "Unknown dataset %s\n", dataset_file);
    return 1;
  }

  nbEle = read_file(dataset_file, &data);
  if (nbEle == 0) {
    fprintf(stderr, "Failed to read dataset %s\n", dataset_file);
    return 1;
  }

  for (size_t i = 0; i < n_bounds; i++) {
    double error_bound = bounds[i];
    uint32_t compression_times[MAX_ITERATIONS];
    uint32_t decompression_times[MAX_ITERATIONS];
    int iteration = 0;
    bool confidence_interval_reached = false;

    while (iteration < MAX_ITERATIONS && !confidence_interval_reached) {
      void *compressed_data = NULL;
      size_t compressed_size = 0;
      double start, end;

      // Compression
      start = omp_get_wtime();
      int rtn = sperr_comp_3d(data, is_float, dims[0], dims[1], dims[2], 256,
                              256, 256, 0, error_bound, num_threads,
                              &compressed_data, &compressed_size);
      end = omp_get_wtime();
      compression_times[iteration] = (uint32_t)((end - start) * 1e6);

      if (rtn != 0) {
        fprintf(stderr, "Compression error with code %d\n", rtn);
        break;
      }

      // Decompression
      void *decompressed_data = NULL;
      size_t out_dimx, out_dimy, out_dimz;
      start = omp_get_wtime();
      rtn = sperr_decomp_3d(compressed_data, compressed_size, is_float,
                            num_threads, &out_dimx, &out_dimy, &out_dimz,
                            &decompressed_data);
      end = omp_get_wtime();
      decompression_times[iteration] = (uint32_t)((end - start) * 1e6);

      if (rtn != 0) {
        fprintf(stderr, "Decompression error with code %d\n", rtn);
        free(compressed_data);
        break;
      }

      // Calculate metrics
      double max_error = 0.0, mse = 0.0, psnr, nrmse;
      double value_range, value_min, value_max;
      double sum_squared_diff = 0.0;

      if (is_float) {
        float *float_data = (float *)data;
        float *float_decompressed = (float *)decompressed_data;
        value_range = value_min = value_max = float_data[0];

#pragma omp parallel for reduction(max : max_error)                            \
    reduction(+ : sum_squared_diff) reduction(min : value_min)                 \
    reduction(max : value_max)
        for (size_t i = 0; i < nbEle; i++) {
          double diff = fabs(float_data[i] - float_decompressed[i]);
          max_error = fmax(max_error, diff);
          sum_squared_diff += diff * diff;
          value_min = fmin(value_min, float_data[i]);
          value_max = fmax(value_max, float_data[i]);
        }
      } else {
        double *double_data = (double *)data;
        double *double_decompressed = (double *)decompressed_data;
        value_range = value_min = value_max = double_data[0];

#pragma omp parallel for reduction(max : max_error)                            \
    reduction(+ : sum_squared_diff) reduction(min : value_min)                 \
    reduction(max : value_max)
        for (size_t i = 0; i < nbEle; i++) {
          double diff = fabs(double_data[i] - double_decompressed[i]);
          max_error = fmax(max_error, diff);
          sum_squared_diff += diff * diff;
          value_min = fmin(value_min, double_data[i]);
          value_max = fmax(value_max, double_data[i]);
        }
      }

      value_range = value_max - value_min;
      mse = sum_squared_diff / nbEle;
      psnr = 20 * log10(value_range) - 10 * log10(mse);
      nrmse = sqrt(mse) / value_range;

      double compression_ratio =
          (double)(nbEle * (is_float ? sizeof(float) : sizeof(double))) /
          compressed_size;

      // Write metrics to CSV file
      FILE *csv_file = fopen("compression_metrics.csv", "a");
      if (csv_file == NULL) {
        fprintf(stderr, "Error opening CSV file\n");
      } else {
        fprintf(
            csv_file,
            "SPERR-OMP,%s,%e,%d,%f,%f,%f,%f,%f,%f,%zu,%f,%zu,%zu,%u,%u,%d\n",
            dataset_file, error_bound, iteration,
            (double)nbEle / (compression_times[iteration] / 1e6),
            (double)nbEle / (decompression_times[iteration] / 1e6), max_error,
            mse, psnr, nrmse, compressed_size, compression_ratio,
            nbEle * (is_float ? sizeof(float) : sizeof(double)),
            compressed_size, compression_times[iteration],
            decompression_times[iteration], num_threads);
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
  }

  free(data);
  return 0;
}