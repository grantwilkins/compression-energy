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
#include <stdbool.h>
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
  printf("BERGINN");
  if (argc != 2) {
    printf("Usage: ./szx_serial.out <dataset_file>\n");
    return 1;
  }
  double bounds[] = {1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1};
  size_t n_bounds = sizeof(bounds) / sizeof(bounds[0]);

  char dataset_file[640];
  sprintf(dataset_file,"%s",argv[1]);
  int ndims = 0;
  void *data;
  int status = 0;
  size_t nbEle;
  size_t dims[5] = {0};
  int data_type;

  if (strstr(dataset_file, "nyx") != NULL) {
    dims[0] = 512;
    dims[1] = 512;
    dims[2] = 512;
    data = SZx_readFloatData(dataset_file, &nbEle, &status);
    data_type = SZx_FLOAT;
  } else if (strstr(dataset_file, "hacc") != NULL) {
    dims[0] = 1073726487;
    data = SZx_readFloatData(dataset_file, &nbEle, &status);
    data_type = SZx_FLOAT;
  } else if (strstr(dataset_file, "s3d") != NULL) {
    dims[0] = 11;
    dims[1] = 500;
    dims[2] = 500;
    dims[3] = 500;
    data = SZx_readDoubleData(dataset_file, &nbEle, &status);
    data_type = SZx_DOUBLE;
  } else if (strstr(dataset_file, "miranda") != NULL) {
    dims[0] = 3072;
    dims[1] = 3072;
    dims[2] = 3072;
    data = SZx_readFloatData(dataset_file, &nbEle, &status);
    data_type = SZx_FLOAT;
  } else {
    fprintf(stderr, "Unknown dataset %s\n", dataset_file);
    return 1;
  }

  if (status != SZx_SCES) {
    printf("Error: data file %s cannot be read!\n", dataset_file);
    return 1;
  }

  uint32_t compression_times[MAX_ITERATIONS];
  uint32_t decompression_times[MAX_ITERATIONS];
  int iteration = 0;
  bool confidence_interval_reached = false;
  for (size_t i = 0; i < n_bounds; i++) {
    double error_bound = bounds[i];
    while (iteration < MAX_ITERATIONS && !confidence_interval_reached) {
      size_t outSize;
      clock_t start, end;
      double cpu_time_used;

      // Compression
      start = clock();
      unsigned char *compressed_data = SZx_fast_compress_args(
          SZx_NO_BLOCK_FAST_CMPR, data_type, data, &outSize, REL, 0,
          error_bound, 0, 0, dims[4], dims[3], dims[2], dims[1], dims[0]);
      end = clock();
      cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
      compression_times[iteration] =
          (uint32_t)(cpu_time_used * 1000000); // Convert to microseconds

      // Decompression
      start = clock();
      void *decompressed_data = SZx_fast_decompress(
          SZx_NO_BLOCK_FAST_CMPR, data_type, compressed_data, outSize, dims[4],
          dims[3], dims[2], dims[1], dims[0]);
      end = clock();
      cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
      decompression_times[iteration] =
          (uint32_t)(cpu_time_used * 1000000); // Convert to microseconds

      // Calculate metrics
      double max_error = 0.0, mse = 0.0, psnr, nrmse;
      double value_range, value_min, value_max;
      double sum_squared_diff = 0.0;

      if (data_type == SZx_FLOAT) {
        float *float_data = (float *)data;
        float *float_decompressed = (float *)decompressed_data;
        value_range = value_min = value_max = float_data[0];

        for (size_t i = 0; i < nbEle; i++) {
          double diff = fabs(float_data[i] - float_decompressed[i]);
          if (diff > max_error)
            max_error = diff;
          sum_squared_diff += diff * diff;
          if (float_data[i] < value_min)
            value_min = float_data[i];
          if (float_data[i] > value_max)
            value_max = float_data[i];
        }
      } else { // SZx_DOUBLE
        double *double_data = (double *)data;
        double *double_decompressed = (double *)decompressed_data;
        value_range = value_min = value_max = double_data[0];

        for (size_t i = 0; i < nbEle; i++) {
          double diff = fabs(double_data[i] - double_decompressed[i]);
          if (diff > max_error)
            max_error = diff;
          sum_squared_diff += diff * diff;
          if (double_data[i] < value_min)
            value_min = double_data[i];
          if (double_data[i] > value_max)
            value_max = double_data[i];
        }
      }

      value_range = value_max - value_min;
      mse = sum_squared_diff / nbEle;
      psnr = 20 * log10(value_range) - 10 * log10(mse);
      nrmse = sqrt(mse) / value_range;

      double compression_ratio =
          (double)(nbEle *
                   (data_type == SZx_FLOAT ? sizeof(float) : sizeof(double))) /
          outSize;

      // Write metrics to CSV file
      FILE *csv_file = fopen("compression_metrics.csv", "a");
      if (csv_file == NULL) {
        fprintf(stderr, "Error opening CSV file\n");
      } else {
        fprintf(
            csv_file, "SZx,%s,%e,%d,%f,%f,%f,%f,%f,%f,%zu,%f,%zu,%zu,%u,%u\n",
            dataset_file, error_bound, iteration,
            (double)nbEle /
                (compression_times[iteration] / 1000000.0), // compression_rate
            (double)nbEle / (decompression_times[iteration] /
                             1000000.0), // decompression_rate
            max_error, mse, psnr, nrmse, outSize, compression_ratio,
            nbEle * (data_type == SZx_FLOAT ? sizeof(float) : sizeof(double)),
            outSize, compression_times[iteration],
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
  }

  free(data);
  return 0;
}
