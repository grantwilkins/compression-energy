/**
 *  @file test_compress_szx.c
 *  @author [Your Name]
 *  @date [Current Date]
 *  @brief This is an example of using SZx compression interface with similar
 * metrics to libpressio
 */

#include "szx.h"
#include "szx_rw.h"
#include <math.h>
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

int main(int argc, char *argv[]) {
  if (argc != 3) {
    printf("Usage: %s <dataset_file> <error_bound>\n", argv[0]);
    return 1;
  }

  const char *dataset_file = argv[1];
  float errorBound = atof(argv[2]);

  int status = 0;
  size_t nbEle;
  float *data = SZx_readFloatData(dataset_file, &nbEle, &status);
  if (status != SZx_SCES) {
    printf("Error: data file %s cannot be read!\n", dataset_file);
    return 1;
  }

  uint32_t compression_times[MAX_ITERATIONS];
  uint32_t decompression_times[MAX_ITERATIONS];
  int iteration = 0;
  bool confidence_interval_reached = false;

  while (iteration < MAX_ITERATIONS && !confidence_interval_reached) {
    size_t outSize;
    clock_t start, end;
    double cpu_time_used;

    // Compression
    start = clock();
    unsigned char *compressed_data =
        SZx_fast_compress_args(SZx_OPENMP_FAST_CMPR, SZx_FLOAT, data, &outSize,
                               ABS, errorBound, 0.001, 0, 0, 0, 0, 0, 0, nbEle);
    end = clock();
    cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
    compression_times[iteration] =
        (uint32_t)(cpu_time_used * 1000000); // Convert to microseconds

    // Decompression
    start = clock();
    float *decompressed_data = (float *)SZx_fast_decompress(
        SZx_OPENMP_FAST_CMPR, SZx_FLOAT, compressed_data, outSize, 0, 0, 0, 0,
        nbEle);
    end = clock();
    cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
    decompression_times[iteration] =
        (uint32_t)(cpu_time_used * 1000000); // Convert to microseconds

    // Calculate metrics
    double max_error = 0.0, mse = 0.0, psnr, nrmse;
    double value_range = data[0], value_min = data[0], value_max = data[0];
    double sum_squared_diff = 0.0;

    for (size_t i = 0; i < nbEle; i++) {
      double diff = fabs(data[i] - decompressed_data[i]);
      if (diff > max_error)
        max_error = diff;
      sum_squared_diff += diff * diff;
      if (data[i] < value_min)
        value_min = data[i];
      if (data[i] > value_max)
        value_max = data[i];
    }

    value_range = value_max - value_min;
    mse = sum_squared_diff / nbEle;
    psnr = 20 * log10(value_range) - 10 * log10(mse);
    nrmse = sqrt(mse) / value_range;

    double compression_ratio = (double)(nbEle * sizeof(float)) / outSize;

    // Write metrics to CSV file
    FILE *csv_file = fopen("compression_metrics.csv", "a");
    if (csv_file == NULL) {
      fprintf(stderr, "Error opening CSV file\n");
    } else {
      fprintf(csv_file, "SZx,%s,%e,%d,%f,%f,%f,%f,%f,%f,%zu,%f,%zu,%zu,%u,%u\n",
              dataset_file, errorBound, iteration,
              (double)nbEle / (compression_times[iteration] /
                               1000000.0), // compression_rate
              (double)nbEle / (decompression_times[iteration] /
                               1000000.0), // decompression_rate
              max_error, mse, psnr, nrmse, outSize, compression_ratio,
              nbEle * sizeof(float), outSize, compression_times[iteration],
              decompression_times[iteration]);
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

  free(data);
  return 0;
}