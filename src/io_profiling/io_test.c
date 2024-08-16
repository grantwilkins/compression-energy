#include <assert.h>
#include <fcntl.h>
#include <libgen.h>
#include <libpressio.h>
#include <libpressio_ext/io/posix.h>
#include <math.h>
#include <papi.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>

#define MAX_ITERATIONS 10
#define CONFIDENCE_LEVEL 1.96
#define BUFFER_SIZE 8192
#define MAX_powercap_EVENTS 64

// Get current time in seconds
double get_time() {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return ts.tv_sec + ts.tv_nsec / 1e9;
}

// Calculate mean of data
double calculate_mean(double *data, int n) {
  double sum = 0.0;
  for (int i = 0; i < n; i++) {
    sum += data[i];
  }
  return sum / n;
}

// Calculate standard deviation of data
double calculate_std_dev(double *data, int n, double mean) {
  double sum_squared_diff = 0.0;
  for (int i = 0; i < n; i++) {
    double diff = data[i] - mean;
    sum_squared_diff += diff * diff;
  }
  return sqrt(sum_squared_diff / (n - 1));
}

// Check if data is within the confidence interval
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

// Calculate total energy consumption from powercap events in Joules
void calculate_energy(long long *values, int num_events,
                      char event_names[][PAPI_MAX_STR_LEN], int *data_type,
                      double *cpu_energy, double *dram_energy) {
  *cpu_energy = 0.0;
  *dram_energy = 0.0;
  for (int i = 0; i < num_events; i++) {
    if (strstr(event_names[i], "ENERGY_UJ")) {
      if (data_type[i] == PAPI_DATATYPE_UINT64) {
        if (strstr(event_names[i], "SUBZONE") == NULL) {
          *cpu_energy += values[i] / 1.0e6; // Convert from uJ to J
        } else {
          *dram_energy += values[i] / 1.0e6;
        }
      }
    }
  }
}

int main(int argc, char *argv[]) {
  if (argc != 5) {
    fprintf(
        stderr,
        "Usage: %s <compressor> <error_bound> <dataset/field> <io_method>\n",
        argv[0]);
    return 1;
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

  const char *compressor_id = argv[1];
  double relative_error_bound = atof(argv[2]);
  const char *dataset_field = argv[3];
  const char *io_method = argv[4];

  const char *LOCAL_SCRATCH = getenv("LOCAL");
  const char *PERSISTENT_DIR = "/ocean/projects/cis240100p/gwilkins/tmp/";

  char dataset_dir[256];
  char field_name[256];
  char scratch_dataset[512];
  char compressed_filename[512];
  char decompressed_filename[512];
  char original_filename[512];

  // Split dataset/field into directory and field name
  char *slash = strchr(dataset_field, '/');
  if (slash == NULL) {
    fprintf(stderr,
            "Invalid dataset/field format. Expected format: dataset/field\n");
    return 1;
  }
  strncpy(dataset_dir, dataset_field, slash - dataset_field);
  dataset_dir[slash - dataset_field] = '\0';
  strcpy(field_name, slash + 1);

  snprintf(scratch_dataset, sizeof(scratch_dataset), "%s/%s", LOCAL_SCRATCH,
           field_name);

  // Check if directories exist
  struct stat st = {0};
  if (stat(LOCAL_SCRATCH, &st) == -1 || stat(PERSISTENT_DIR, &st) == -1) {
    fprintf(stderr, "Error: Required directories do not exist\n");
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

  // Configure I/O
  struct pressio_options *io_options = pressio_options_new();
  pressio_options_set_string(io_options, "io:path", scratch_dataset);
  if (strcmp(io_method, "hdf5") == 0) {
    pressio_options_set_string(io_options, "hdf5:dataset", "/data");
  } else if (strcmp(io_method, "netcdf") == 0) {
    pressio_options_set_string(io_options, "netcdf:variable", "data");
  }
  pressio_io_set_options(io, io_options);

  size_t ndims;
  struct pressio_data *metadata, *input_data;
  if (strstr(dataset_field, "nyx") != NULL) {
    size_t dims[] = {512, 512, 512};
    ndims = sizeof(dims) / sizeof(dims[0]);
    metadata = pressio_data_new_empty(pressio_float_dtype, ndims, dims);
  } else if (strstr(dataset_field, "hacc") != NULL) {
    size_t dims[] = {1073726487};
    ndims = sizeof(dims) / sizeof(dims[0]);
    metadata = pressio_data_new_empty(pressio_float_dtype, ndims, dims);
  } else if (strstr(dataset_field, "s3d") != NULL) {
    size_t dims[] = {11, 500, 500, 500};
    ndims = sizeof(dims) / sizeof(dims[0]);
    metadata = pressio_data_new_empty(pressio_double_dtype, ndims, dims);
  } else if (strstr(dataset_field, "miranda") != NULL) {
    size_t dims[] = {3072, 3072, 3072};
    ndims = sizeof(dims) / sizeof(dims[0]);
    metadata = pressio_data_new_empty(pressio_float_dtype, ndims, dims);
  } else if (strstr(dataset_field, "cesm") != NULL) {
    size_t dims[] = {26, 1800, 3600};
    ndims = sizeof(dims) / sizeof(dims[0]);
    metadata = pressio_data_new_empty(pressio_float_dtype, ndims, dims);
  } else {
    fprintf(stderr, "Unknown dataset %s\n", dataset_field);
    pressio_compressor_release(compressor);
    pressio_release(library);
    return 1;
  }

  input_data = pressio_io_data_path_read(metadata, scratch_dataset);

  double data_min, data_max, data_range;
  size_t num_elements = pressio_data_num_elements(input_data);
  void *data_ptr = pressio_data_ptr(input_data, NULL);

  if (strstr(dataset_field, "s3d") == NULL) {
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

  // Configure compressor
  struct pressio_options *comp_options = pressio_options_new();
  pressio_options_set_double(comp_options, "pressio:abs", absolute_error_bound);

  const char *metrics_ids[] = {"size"};
  size_t n_metrics_ids = sizeof(metrics_ids) / sizeof(metrics_ids[0]);
  pressio_options_set_string(comp_options, "pressio:metric", "composite");
  pressio_options_set_strings(comp_options, "composite:plugins", n_metrics_ids,
                              metrics_ids);
  if (pressio_compressor_set_options(compressor, comp_options)) {
    fprintf(stderr, "Failed to set compressor options: %s\n",
            pressio_compressor_error_msg(compressor));
    return 1;
  }

  pressio_options_free(comp_options);
  pressio_data_free(input_data);
  free(data_ptr);

  double read_times[MAX_ITERATIONS] = {0};
  double compression_times[MAX_ITERATIONS] = {0};
  double write_compressed_times[MAX_ITERATIONS] = {0};
  double decompression_times[MAX_ITERATIONS] = {0};
  double write_original_times[MAX_ITERATIONS] = {0};
  double read_cpu_energy[MAX_ITERATIONS] = {0};
  double read_dram_energy[MAX_ITERATIONS] = {0};
  double compression_cpu_energy[MAX_ITERATIONS] = {0};
  double compression_dram_energy[MAX_ITERATIONS] = {0};
  double write_compressed_cpu_energy[MAX_ITERATIONS] = {0};
  double write_compressed_dram_energy[MAX_ITERATIONS] = {0};
  double decompression_cpu_energy[MAX_ITERATIONS] = {0};
  double decompression_dram_energy[MAX_ITERATIONS] = {0};
  double write_original_cpu_energy[MAX_ITERATIONS] = {0};
  double write_original_dram_energy[MAX_ITERATIONS] = {0};
  double cpu_energy, dram_energy;
  int iteration = 0;
  bool confidence_interval_reached = false;
  struct pressio_data *compressed_data = NULL;
  struct pressio_data *read_compressed_data = NULL;
  struct pressio_data *decompressed_data = NULL;
  uint64_t compressed_size, decompressed_size, uncompressed_size;
  double bit_rate, compression_ratio;

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
            "Compressor,Dataset,IO Method,Error Bound,Iteration,"
            "Read Time,Compression Time,Write Compressed Time,Read Compressed "
            "Time,Decompression Time,Write Original Time,"
            "Read Energy,Compression Energy,Write Compressed Energy,Read "
            "Compressed Energy,Decompression Energy,Write Original Energy\n");
  }
  fclose(csv_file);
  double start_time = 0.0;

  while (iteration < MAX_ITERATIONS && !confidence_interval_reached) {
    // Read input data from scratch
    assert(PAPI_start(EventSet) == PAPI_OK);
    start_time = get_time();
    input_data = pressio_io_data_path_read(metadata, scratch_dataset);
    read_times[iteration] = get_time() - start_time;
    assert(PAPI_stop(EventSet, values) == PAPI_OK);

    calculate_energy(values, num_events, event_names, data_type, &cpu_energy,
                     &dram_energy);
    read_cpu_energy[iteration] = cpu_energy;
    read_dram_energy[iteration] = dram_energy;

    if (input_data == NULL) {
      fprintf(stderr, "Failed to read dataset %s\n", scratch_dataset);
      pressio_compressor_release(compressor);
      pressio_release(library);
      return 1;
    }

    // Compress data
    assert(PAPI_start(EventSet) == PAPI_OK);
    start_time = get_time();
    compressed_data = pressio_data_new_empty(pressio_byte_dtype, 0, NULL);
    if (pressio_compressor_compress(compressor, input_data, compressed_data)) {
      fprintf(stderr, "Compression failed: %s\n",
              pressio_compressor_error_msg(compressor));
      return 1;
    }
    compression_times[iteration] = get_time() - start_time;
    assert(PAPI_stop(EventSet, values) == PAPI_OK);

    calculate_energy(values, num_events, event_names, data_type, &cpu_energy,
                     &dram_energy);
    compression_cpu_energy[iteration] = cpu_energy;
    compression_dram_energy[iteration] = dram_energy;

    // Write compressed data to persistent storage
    assert(PAPI_start(EventSet) == PAPI_OK);
    start_time = get_time();
    snprintf(compressed_filename, sizeof(compressed_filename), "%s/%s.%s",
             PERSISTENT_DIR, field_name, compressor_id);
    pressio_options_set_string(io_options, "io:path", compressed_filename);
    pressio_io_set_options(io, io_options);
    if (pressio_io_write(io, compressed_data)) {
      fprintf(stderr, "Failed to write compressed data: %s\n",
              pressio_io_error_msg(io));
      return 1;
    }
    write_compressed_times[iteration] = get_time() - start_time;
    assert(PAPI_stop(EventSet, values) == PAPI_OK);

    calculate_energy(values, num_events, event_names, data_type, &cpu_energy,
                     &dram_energy);
    write_compressed_cpu_energy[iteration] = cpu_energy;
    write_compressed_dram_energy[iteration] = dram_energy;

    read_compressed_data = pressio_io_read(io, NULL);
    // Decompress data
    assert(PAPI_start(EventSet) == PAPI_OK);
    start_time = get_time();
    decompressed_data = pressio_data_new_clone(input_data);
    if (pressio_compressor_decompress(compressor, read_compressed_data,
                                      decompressed_data)) {
      fprintf(stderr, "Decompression failed: %s\n",
              pressio_compressor_error_msg(compressor));
      return 1;
    }
    decompression_times[iteration] = get_time() - start_time;
    assert(PAPI_stop(EventSet, values) == PAPI_OK);
    calculate_energy(values, num_events, event_names, data_type, &cpu_energy,
                     &dram_energy);
    decompression_cpu_energy[iteration] = cpu_energy;
    decompression_dram_energy[iteration] = dram_energy;

    // Write original data to persistent storage
    assert(PAPI_start(EventSet) == PAPI_OK);
    start_time = get_time();
    snprintf(original_filename, sizeof(original_filename), "%s/%s.original",
             PERSISTENT_DIR, field_name);
    pressio_options_set_string(io_options, "io:path", original_filename);
    pressio_io_set_options(io, io_options);
    if (pressio_io_write(io, input_data)) {
      fprintf(stderr, "Failed to write original data: %s\n",
              pressio_io_error_msg(io));
      return 1;
    }
    write_original_times[iteration] = get_time() - start_time;
    assert(PAPI_stop(EventSet, values) == PAPI_OK);

    calculate_energy(values, num_events, event_names, data_type, &cpu_energy,
                     &dram_energy);
    write_original_cpu_energy[iteration] = cpu_energy;
    write_original_dram_energy[iteration] = dram_energy;

    struct pressio_options *metrics_results =
        pressio_compressor_get_metrics_results(compressor);
    pressio_options_get_double(metrics_results, "size:bit_rate", &bit_rate);
    pressio_options_get_uinteger64(metrics_results, "size:compressed_size",
                                   &compressed_size);
    pressio_options_get_double(metrics_results, "size:compression_ratio",
                               &compression_ratio);
    pressio_options_get_uinteger64(metrics_results, "size:decompressed_size",
                                   &decompressed_size);
    pressio_options_get_uinteger64(metrics_results, "size:uncompressed_size",
                                   &uncompressed_size);
    // Write results to CSV file for this iteration
    csv_file = fopen("io_test_results.csv", "a");
    if (csv_file == NULL) {
      fprintf(stderr, "Error opening CSV file\n");
      return 1;
    }
    fprintf(csv_file,
            "%s,%s,%s,%e,%e,%d,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,"
            "%d,%e,%d,%d\n",
            compressor_id, dataset_field, io_method, relative_error_bound,
            absolute_error_bound, iteration, read_times[iteration],
            compression_times[iteration], write_compressed_times[iteration],
            decompression_times[iteration], write_original_times[iteration],
            read_cpu_energy[iteration], read_dram_energy[iteration],
            compression_cpu_energy[iteration],
            compression_dram_energy[iteration],
            write_compressed_cpu_energy[iteration],
            write_compressed_dram_energy[iteration],
            decompression_cpu_energy[iteration],
            decompression_dram_energy[iteration],
            write_original_cpu_energy[iteration],
            write_original_dram_energy[iteration], bit_rate, compressed_size,
            compression_ratio, decompressed_size, uncompressed_size);
    fclose(csv_file);

    // Clean up iteration-specific data
    pressio_data_free(input_data);
    pressio_data_free(compressed_data);
    pressio_data_free(read_compressed_data);
    pressio_data_free(decompressed_data);

    iteration++;

    // Check if we've reached the confidence interval
    if (iteration >= 5) {
      confidence_interval_reached =
          within_confidence_interval(read_times, iteration) &&
          within_confidence_interval(compression_times, iteration) &&
          within_confidence_interval(write_compressed_times, iteration) &&
          within_confidence_interval(decompression_times, iteration);
    }
  }

  // Clean up
  pressio_options_free(comp_options);
  pressio_options_free(io_options);
  pressio_compressor_release(compressor);
  pressio_io_free(io);
  pressio_release(library);

  printf("Compression and decompression completed successfully.\n");
  printf("Results written to compression_results.csv\n");
  printf("Compressed file: %s\n", compressed_filename);
  printf("Decompressed file: %s\n", decompressed_filename);

  // Cleanup: Remove files from scratch and persistent memory
  printf("Cleaning up temporary files...\n");

  // Remove the copy in global scratch
  if (remove(scratch_dataset) == 0) {
    printf("Removed scratch file: %s\n", scratch_dataset);
  } else {
    perror("Failed to remove scratch file");
  }

  // Remove the compressed file from persistent memory
  if (remove(compressed_filename) == 0) {
    printf("Removed compressed file: %s\n", compressed_filename);
  } else {
    perror("Failed to remove compressed file");
  }

  // Remove the decompressed file from persistent memory
  if (remove(decompressed_filename) == 0) {
    printf("Removed decompressed file: %s\n", decompressed_filename);
  } else {
    perror("Failed to remove decompressed file");
  }

  printf("Cleanup completed.\n");

  return 0;
}
