#include <assert.h>
#include <libpressio.h>
#include <libpressio_ext/io/posix.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

int main(int argc, char *argv[]) {
  if (argc != 5) {
    fprintf(stderr, "Usage: %s <compressor> <error_bound> <dataset> <field>\n",
            argv[0]);
    return 1;
  }

  const char *compressor_id = argv[1];
  const char *dataset_file = argv[2];
  double relative_error_bound = atof(argv[3]);
  const char *datadir = "/ocean/projects/cis240100p/gwilkins/";
  const char *cluster_name = getenv("CLUSTER_NAME");

  const char *PERSISTENT_DIR = "/ocean/projects/cis240100p/gwilkins/tmp/";
  char input_filename[512];
  char output_filename[512];

  snprintf(input_filename, sizeof(input_filename), "%s/%s/%s", PERSISTENT_DIR,
           dataset, field);
  snprintf(output_filename, sizeof(output_filename), "%s/%s_comp/%s.%g.%s",
           PERSISTENT_DIR, dataset, field, relative_error_bound, compressor_id);

  // Initialize LibPressio
  struct pressio *library = pressio_instance();

  struct pressio_compressor *compressor =
      pressio_get_compressor(library, compressor_id);
  if (compressor == NULL) {
    fprintf(stderr, "Failed to get compressor %s: %s\n", compressor_id,
            pressio_error_msg(library));
    pressio_release(library);
    return 1;
  }

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
  struct pressio_options *options = pressio_options_new();
  pressio_options_set_double(options, "pressio:abs", absolute_error_bound);

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

  // Compress data
  struct pressio_data *compressed_data =
      pressio_data_new_empty(pressio_byte_dtype, 0, NULL);
  if (pressio_compressor_compress(compressor, input_data, compressed_data)) {
    fprintf(stderr, "Compression failed: %s\n",
            pressio_compressor_error_msg(compressor));
    return 1;
  }

  // Write compressed data to persistent storage
  struct pressio_io *io = pressio_get_io(library, "posix");
  struct pressio_options *io_options = pressio_options_new();
  pressio_options_set_string(io_options, "io:path", output_filename);
  pressio_io_set_options(io, io_options);
  if (pressio_io_write(io, compressed_data)) {
    fprintf(stderr, "Failed to write compressed data: %s\n",
            pressio_io_error_msg(io));
    return 1;
  }

  printf("Compressed data written to: %s\n", output_filename);

  // Clean up
  pressio_data_free(input_data);
  pressio_data_free(compressed_data);
  pressio_data_free(metadata);
  pressio_options_free(io_options);
  pressio_compressor_release(compressor);
  pressio_io_free(io);
  pressio_release(library);

  return 0;
}