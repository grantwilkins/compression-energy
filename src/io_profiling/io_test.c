#include <fcntl.h>
#include <libgen.h>
#include <libpressio.h>
#include <libpressio_ext/io/posix.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>

#define MAX_ITERATIONS 25
#define CONFIDENCE_LEVEL 1.96
#define BUFFER_SIZE 8192

const char *GLOBAL_SCRATCH = "/lcrc/globalscratch/ac.gwilkins";
const char *PERSISTENT_DIR = "/lcrc/project/ECP-EZ/ac.gwilkins";

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

int copy_file(const char *src, const char *dst) {
  int fd_src, fd_dst;
  char buffer[BUFFER_SIZE];
  ssize_t bytes_read, bytes_written;
  int result = 0;

  fd_src = open(src, O_RDONLY);
  if (fd_src == -1) {
    perror("Error opening source file");
    return -1;
  }

  fd_dst = open(dst, O_WRONLY | O_CREAT | O_TRUNC, 0666);
  if (fd_dst == -1) {
    perror("Error opening destination file");
    close(fd_src);
    return -1;
  }

  while ((bytes_read = read(fd_src, buffer, BUFFER_SIZE)) > 0) {
    bytes_written = write(fd_dst, buffer, bytes_read);
    if (bytes_written != bytes_read) {
      perror("Error writing to destination file");
      result = -1;
      break;
    }
  }

  if (bytes_read == -1) {
    perror("Error reading from source file");
    result = -1;
  }

  close(fd_src);
  close(fd_dst);
  return result;
}

int main(int argc, char *argv[]) {
  if (argc != 5) {
    fprintf(stderr,
            "Usage: %s <compressor> <error_bound> <dataset> <io_method>\n",
            argv[0]);
    return 1;
  }

  const char *compressor_id = argv[1];
  double error_bound = atof(argv[2]);
  const char *dataset_file = argv[3];
  const char *io_method = argv[4];

  char compressed_filename[512];
  char decompressed_filename[512];

  // Create directories if they don't exist
  mkdir(GLOBAL_SCRATCH, 0777);
  mkdir(PERSISTENT_DIR, 0777);

  // Copy dataset to global scratch
  char scratch_dataset[512];
  snprintf(scratch_dataset, sizeof(scratch_dataset), "%s/%s", GLOBAL_SCRATCH,
           basename((char *)dataset_file));
  if (copy_file(dataset_file, scratch_dataset) != 0) {
    fprintf(stderr, "Failed to copy dataset to global scratch\n");
    return 1;
  }

  // Initialize LibPressio
  struct pressio *library = pressio_instance();
  struct pressio_compressor *compressor =
      pressio_get_compressor(library, compressor_id);
  struct pressio_io *io = pressio_get_io(library, io_method);

  if (compressor == NULL) {
    fprintf(stderr, "Failed to get compressor %s: %s\n", compressor_id,
            pressio_error_msg(library));
    return 1;
  }

  if (io == NULL) {
    fprintf(stderr, "Failed to get I/O method %s: %s\n", io_method,
            pressio_error_msg(library));
    return 1;
  }

  // Configure compressor
  struct pressio_options *comp_options = pressio_options_new();
  pressio_options_set_double(comp_options, "pressio:abs", error_bound);
  if (pressio_compressor_set_options(compressor, comp_options)) {
    fprintf(stderr, "Failed to set compressor options: %s\n",
            pressio_compressor_error_msg(compressor));
    return 1;
  }

  // Configure I/O
  struct pressio_options *io_options = pressio_options_new();
  pressio_options_set_string(io_options, "io:path", scratch_dataset);
  if (strcmp(io_method, "hdf5") == 0) {
    pressio_options_set_string(io_options, "hdf5:dataset", "/data");
  } else if (strcmp(io_method, "netcdf") == 0) {
    pressio_options_set_string(io_options, "netcdf:variable", "data");
  }
  pressio_io_set_options(io, io_options);

  double read_times[MAX_ITERATIONS] = {0};
  double compression_times[MAX_ITERATIONS] = {0};
  double write_compressed_times[MAX_ITERATIONS] = {0};
  double read_compressed_times[MAX_ITERATIONS] = {0};
  double decompression_times[MAX_ITERATIONS] = {0};
  double write_decompressed_times[MAX_ITERATIONS] = {0};

  int iteration = 0;
  bool confidence_interval_reached = false;

  // Open CSV file for writing
  FILE *csv_file = fopen("compression_results.csv", "a");
  if (csv_file == NULL) {
    fprintf(stderr, "Error opening CSV file\n");
    return 1;
  }

  // Write CSV header if file is empty
  fseek(csv_file, 0, SEEK_END);
  if (ftell(csv_file) == 0) {
    fprintf(csv_file,
            "Compressor,Dataset,IO_Method,Error_Bound,Iteration,Read_Time,"
            "Compression_Time,Write_Compressed_Time,Read_Compressed_Time,"
            "Decompression_Time,Write_Decompressed_Time\n");
  }

  while (iteration < MAX_ITERATIONS && !confidence_interval_reached) {
    // Read input data
    double start_time = get_time();
    struct pressio_data *input_data = pressio_io_read(io, NULL);
    read_times[iteration] = get_time() - start_time;
    if (input_data == NULL) {
      fprintf(stderr, "Failed to read input data: %s\n",
              pressio_io_error_msg(io));
      fclose(csv_file);
      return 1;
    }

    // Compress data
    start_time = get_time();
    struct pressio_data *compressed_data =
        pressio_data_new_empty(pressio_byte_dtype, 0, NULL);
    if (pressio_compressor_compress(compressor, input_data, compressed_data)) {
      fprintf(stderr, "Compression failed: %s\n",
              pressio_compressor_error_msg(compressor));
      fclose(csv_file);
      return 1;
    }
    compression_times[iteration] = get_time() - start_time;

    // Write compressed data to persistent storage
    snprintf(compressed_filename, sizeof(compressed_filename),
             "%s/%s.compressed", PERSISTENT_DIR,
             basename((char *)dataset_file));
    pressio_options_set_string(io_options, "io:path", compressed_filename);
    pressio_io_set_options(io, io_options);
    start_time = get_time();
    if (pressio_io_write(io, compressed_data)) {
      fprintf(stderr, "Failed to write compressed data: %s\n",
              pressio_io_error_msg(io));
      fclose(csv_file);
      return 1;
    }
    write_compressed_times[iteration] = get_time() - start_time;

    // Read compressed data back from persistent storage
    start_time = get_time();
    struct pressio_data *read_compressed_data = pressio_io_read(io, NULL);
    read_compressed_times[iteration] = get_time() - start_time;
    if (read_compressed_data == NULL) {
      fprintf(stderr, "Failed to read compressed data: %s\n",
              pressio_io_error_msg(io));
      fclose(csv_file);
      return 1;
    }

    // Decompress data
    start_time = get_time();
    struct pressio_data *decompressed_data = pressio_data_new_clone(input_data);
    if (pressio_compressor_decompress(compressor, read_compressed_data,
                                      decompressed_data)) {
      fprintf(stderr, "Decompression failed: %s\n",
              pressio_compressor_error_msg(compressor));
      fclose(csv_file);
      return 1;
    }
    decompression_times[iteration] = get_time() - start_time;

    // Write decompressed data to persistent storage
    snprintf(decompressed_filename, sizeof(decompressed_filename),
             "%s/%s.decompressed", PERSISTENT_DIR,
             basename((char *)dataset_file));
    pressio_options_set_string(io_options, "io:path", decompressed_filename);
    pressio_io_set_options(io, io_options);
    start_time = get_time();
    if (pressio_io_write(io, decompressed_data)) {
      fprintf(stderr, "Failed to write decompressed data: %s\n",
              pressio_io_error_msg(io));
      fclose(csv_file);
      return 1;
    }
    write_decompressed_times[iteration] = get_time() - start_time;

    // Write results to CSV file for this iteration
    fprintf(csv_file, "%s,%s,%s,%f,%d,%f,%f,%f,%f,%f,%f\n", compressor_id,
            dataset_file, io_method, error_bound, iteration,
            read_times[iteration], compression_times[iteration],
            write_compressed_times[iteration], read_compressed_times[iteration],
            decompression_times[iteration],
            write_decompressed_times[iteration]);

    iteration++;

    // Check if we've reached the confidence interval
    if (iteration >= 5) {
      confidence_interval_reached =
          within_confidence_interval(read_times, iteration) &&
          within_confidence_interval(compression_times, iteration) &&
          within_confidence_interval(write_compressed_times, iteration) &&
          within_confidence_interval(read_compressed_times, iteration) &&
          within_confidence_interval(decompression_times, iteration) &&
          within_confidence_interval(write_decompressed_times, iteration);
    }
  }

  // Close CSV file
  fclose(csv_file);

  pressio_data_free(input_data);
  pressio_data_free(compressed_data);
  pressio_data_free(read_compressed_data);
  pressio_data_free(decompressed_data);
  pressio_options_free(comp_options);
  pressio_options_free(io_options);
  pressio_compressor_release(compressor);
  pressio_io_free(io);
  pressio_release(library);

  // Remove the copy in global scratch
  remove(scratch_dataset);

  printf("Compression and decompression completed successfully.\n");
  printf("Results written to compression_results.csv\n");
  printf("Compressed file: %s\n", compressed_filename);
  printf("Decompressed file: %s\n", decompressed_filename);

  return 0;
}