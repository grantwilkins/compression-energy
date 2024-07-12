#include <libpressio.h>
#include <libpressio_ext/io/posix.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {
  // read in the dataset
  size_t dims[] = {100, 500, 500};
  size_t ndims = sizeof(dims) / sizeof(dims[0]);
  struct pressio_data *metadata =
      pressio_data_new_empty(pressio_float_dtype, ndims, dims);
  struct pressio_data *input_data = pressio_io_data_path_read(
      metadata, "/home/gwilkins/data/nyx/temperature.f32");

  // create output locations
  struct pressio_data *compressed =
      pressio_data_new_empty(pressio_byte_dtype, 0, NULL);
  struct pressio_data *output = pressio_data_new_clone(input_data);

  // get the compressor
  struct pressio *library = pressio_instance();

  const char *compressor_ids[] = {"sz", "sz3", "mgard", "mgardx", "zfp"};
  size_t n_compressors = sizeof(compressor_ids) / sizeof(compressor_ids[0]);
  struct pressio_compressor *compressors[5];

  for (size_t i = 0; i < n_compressors; ++i) {
    compressors[i] = pressio_get_compressor(library, compressor_ids[i]);
    if (compressors[i] == NULL) {
      fprintf(stderr, "Failed to get compressor %s: %s\n", compressor_ids[i],
              pressio_error_msg(library));
      // Handle error (e.g., continue to next compressor or exit)
    }
  }

  // configure metrics for the compressors
  struct pressio_options *metrics_options = pressio_options_new();
  const char *metrics_ids[] = {"time", "size", "error_stat"};
  size_t n_metrics_ids = sizeof(metrics_ids) / sizeof(metrics_ids[0]);
  pressio_options_set_string(metrics_options, "pressio:metric", "composite");
  pressio_options_set_strings(metrics_options, "composite:plugins",
                              n_metrics_ids, metrics_ids);

  double bounds[] = {1e-6, 5e-6, 1e-5, 5e-5, 1e-4, 5e-4,
                     1e-3, 5e-3, 1e-2, 5e-2, 1e-1};
  size_t n_bounds = sizeof(bounds) / sizeof(bounds[0]);

  for (size_t c = 0; c < n_compressors; ++c) {
    struct pressio_compressor *comp = compressors[c];
    if (comp == NULL)
      continue; // Skip if compressor wasn't initialized

    // verify that options passed exist
    if (pressio_compressor_check_options(comp, metrics_options)) {
      fprintf(stderr, "%s\n", pressio_compressor_error_msg(comp));
      continue;
    }

    // set the new options for the compressor
    if (pressio_compressor_set_options(comp, metrics_options)) {
      fprintf(stderr, "%s\n", pressio_compressor_error_msg(comp));
      continue;
    }

    printf("Compressor: %s\n", compressor_ids[c]);

    for (size_t i = 0; i < n_bounds; ++i) {
      // configure the compressor error bound
      struct pressio_options *options = pressio_options_new();
      pressio_options_set_double(options, "pressio:rel", bounds[i]);

      // verify that options passed exist
      if (pressio_compressor_check_options(comp, options)) {
        fprintf(stderr, "%s\n", pressio_compressor_error_msg(comp));
        pressio_options_free(options);
        continue;
      }

      // set the new options for the compressor
      if (pressio_compressor_set_options(comp, options)) {
        fprintf(stderr, "%s\n", pressio_compressor_error_msg(comp));
        pressio_options_free(options);
        continue;
      }
      pressio_options_free(options);

      // Variables for statistical analysis
      double compression_times[50] = {0};
      double decompression_times[50] = {0};
      int run_count = 0;
      double compression_mean = 0, decompression_mean = 0;
      double compression_variance = 0, decompression_variance = 0;
      double compression_z_score = 0, decompression_z_score = 0;
      const double confidence_interval = 1.96; // 95% confidence interval

      do {
        // run the compression and decompression
        if (pressio_compressor_compress(comp, input_data, compressed)) {
          fprintf(stderr, "%s\n", pressio_compressor_error_msg(comp));
          break;
        }
        if (pressio_compressor_decompress(comp, compressed, output)) {
          fprintf(stderr, "%s\n", pressio_compressor_error_msg(comp));
          break;
        }

        // get metrics results
        struct pressio_options *metrics_results =
            pressio_compressor_get_metrics_results(comp);

        double compression_time, decompression_time;
        pressio_options_get_double(metrics_results, "time:compress",
                                   &compression_time);
        pressio_options_get_double(metrics_results, "time:decompress",
                                   &decompression_time);

        compression_times[run_count] = compression_time;
        decompression_times[run_count] = decompression_time;

        run_count++;

        // Calculate mean and variance
        double old_compression_mean = compression_mean;
        double old_decompression_mean = decompression_mean;
        compression_mean += (compression_time - compression_mean) / run_count;
        decompression_mean +=
            (decompression_time - decompression_mean) / run_count;
        compression_variance += (compression_time - old_compression_mean) *
                                (compression_time - compression_mean);
        decompression_variance +=
            (decompression_time - old_decompression_mean) *
            (decompression_time - decompression_mean);

        double compression_std_dev = 0;
        double decompression_std_dev = 0;
        if (run_count > 1) {
          compression_std_dev = sqrt(compression_variance / (run_count - 1));
          decompression_std_dev =
              sqrt(decompression_variance / (run_count - 1));

          // Add a small epsilon to avoid division by zero
          double epsilon = 1e-10;
          compression_z_score =
              (compression_time - compression_mean) /
              (compression_std_dev / sqrt(run_count) + epsilon);
          decompression_z_score =
              (decompression_time - decompression_mean) /
              (decompression_std_dev / sqrt(run_count) + epsilon);
        } else {
          // Set z-scores to a large value for the first run to ensure the loop
          // continues
          compression_z_score = decompression_z_score = 1000;
        }

        pressio_options_free(metrics_results);
        printf(
            "Run %d: compression_z_score=%1.2f decompression_z_score=%1.2f\n",
            run_count, compression_z_score, decompression_z_score);

      } while ((fabs(compression_z_score) > confidence_interval ||
                fabs(decompression_z_score) > confidence_interval) &&
               run_count < 50);

      // print out the final metrics results
      struct pressio_options *final_metrics =
          pressio_compressor_get_metrics_results(comp);
      char *metrics_results_str = pressio_options_to_string(final_metrics);
      printf("bound=%1.1e\nruns=%d\n%s\n", bounds[i], run_count,
             metrics_results_str);
      free(metrics_results_str);
      pressio_options_free(final_metrics);
    }
    printf("\n"); // Add a newline between compressors for readability
  }

  pressio_options_free(metrics_options);

  pressio_data_free(metadata);
  pressio_data_free(input_data);
  pressio_data_free(compressed);
  pressio_data_free(output);
  for (size_t i = 0; i < n_compressors; ++i) {
    if (compressors[i] != NULL) {
      pressio_compressor_release(compressors[i]);
    }
  }
  pressio_release(library);
  return 0;
}