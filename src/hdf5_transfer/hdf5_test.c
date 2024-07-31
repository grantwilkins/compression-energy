#include <hdf5.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define FILE_NAME_ORIGINAL "data_transfer_original.h5"
#define FILE_NAME_COMPRESSED "data_transfer_compressed.h5"
#define DATASET_NAME "data"
#define MAX_ITERATIONS 25
#define CONFIDENCE_LEVEL 1.96

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

void perform_io_test(const char *file_name, void *data, long data_size,
                     double *write_times) {
  int iteration = 0;
  bool confidence_interval_reached = false;

  while (iteration < MAX_ITERATIONS && !confidence_interval_reached) {
    // Create a new HDF5 file (or open existing and truncate)
    hid_t file_id =
        H5Fcreate(file_name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    // Create dataspace for the dataset
    hsize_t dims[1] = {data_size};
    hid_t dataspace = H5Screate_simple(1, dims, NULL);

    // Create dataset
    hid_t dataset =
        H5Dcreate2(file_id, DATASET_NAME, H5T_NATIVE_UCHAR, dataspace,
                   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    // Write data
    double start_time = get_time();
    H5Dwrite(dataset, H5T_NATIVE_UCHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    write_times[iteration] = get_time() - start_time;

    // Close and release resources
    H5Dclose(dataset);
    H5Sclose(dataspace);
    H5Fclose(file_id);

    iteration++;

    // Check if we've reached the confidence interval
    if (iteration >= 5) {
      confidence_interval_reached =
          within_confidence_interval(write_times, iteration);
    }
  }

  printf("Number of iterations for %s: %d\n", file_name, iteration);
}

int main(int argc, char *argv[]) {
  if (argc != 3) {
    fprintf(stderr, "Usage: %s <original_file> <compressed_file>\n", argv[0]);
    return 1;
  }

  const char *original_file = argv[1];
  const char *compressed_file = argv[2];

  // Read original data
  FILE *f_original = fopen(original_file, "rb");
  if (!f_original) {
    fprintf(stderr, "Error opening original file\n");
    return 1;
  }
  fseek(f_original, 0, SEEK_END);
  long original_size = ftell(f_original);
  fseek(f_original, 0, SEEK_SET);
  void *original_data = malloc(original_size);
  fread(original_data, 1, original_size, f_original);
  fclose(f_original);

  // Read compressed data
  FILE *f_compressed = fopen(compressed_file, "rb");
  if (!f_compressed) {
    fprintf(stderr, "Error opening compressed file\n");
    free(original_data);
    return 1;
  }
  fseek(f_compressed, 0, SEEK_END);
  long compressed_size = ftell(f_compressed);
  fseek(f_compressed, 0, SEEK_SET);
  void *compressed_data = malloc(compressed_size);
  fread(compressed_data, 1, compressed_size, f_compressed);
  fclose(f_compressed);

  double original_write_times[MAX_ITERATIONS] = {0};
  double compressed_write_times[MAX_ITERATIONS] = {0};

  // Perform I/O test for original data
  perform_io_test(FILE_NAME_ORIGINAL, original_data, original_size,
                  original_write_times);

  // Perform I/O test for compressed data
  perform_io_test(FILE_NAME_COMPRESSED, compressed_data, compressed_size,
                  compressed_write_times);

  free(original_data);
  free(compressed_data);

  // Calculate and print average times
  double avg_original_write_time =
      calculate_mean(original_write_times, MAX_ITERATIONS);
  double avg_compressed_write_time =
      calculate_mean(compressed_write_times, MAX_ITERATIONS);

  printf("Average original data write time: %f seconds\n",
         avg_original_write_time);
  printf("Average compressed data write time: %f seconds\n",
         avg_compressed_write_time);

  // Write results to CSV
  FILE *csv_file = fopen("hdf5_transfer_results.csv", "a");
  if (csv_file == NULL) {
    fprintf(stderr, "Error opening CSV file\n");
  } else {
    fprintf(csv_file, "%s,%s,%f,%f\n", original_file, compressed_file,
            avg_original_write_time, avg_compressed_write_time);
    fclose(csv_file);
  }

  return 0;
}