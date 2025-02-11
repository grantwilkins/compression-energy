// inflate.c
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>

#define ORIG_DIM 512  // The original dataset is assumed to be 512x512x512

int main(int argc, char *argv[]) {
    if(argc != 4) {
        fprintf(stderr, "Usage: %s <input_file> <expansion_factor> <output_file>\n", argv[0]);
        fprintf(stderr, "  (expansion_factor should be a perfect cube, e.g., 8, 27, 64, etc.)\n");
        return EXIT_FAILURE;
    }
    
    const char *input_file = argv[1];
    long expansion_factor = atol(argv[2]);
    const char *output_file = argv[3];
    
    // Compute the per-dimension replication factor r = cbrt(expansion_factor)
    double d_r = cbrt((double) expansion_factor);
    size_t r = (size_t)(d_r + 1e-6);  // add a tiny epsilon to help rounding
    if(r * r * r != (size_t) expansion_factor) {
        fprintf(stderr, "Warning: expansion_factor %ld is not a perfect cube; using r = %zu\n",
                expansion_factor, r);
    }
    
    // Original dimensions
    size_t orig_x = ORIG_DIM;
    size_t orig_y = ORIG_DIM;
    size_t orig_z = ORIG_DIM;
    size_t orig_total = orig_x * orig_y * orig_z;
    
    // New (inflated) dimensions
    size_t new_x = orig_x * r;
    size_t new_y = orig_y * r;
    size_t new_z = orig_z * r;
    size_t new_total = new_x * new_y * new_z;
    
    printf("Original dimensions: %zux%zux%zu (%zu elements)\n", orig_x, orig_y, orig_z, orig_total);
    printf("Replication factor: %zu\n", r);
    printf("New dimensions: %zux%zux%zu (%zu elements)\n", new_x, new_y, new_z, new_total);
    
    // Allocate memory for the original dataset.
    float *orig_data = (float *)malloc(orig_total * sizeof(float));
    if (orig_data == NULL) {
        fprintf(stderr, "Error: failed to allocate memory for original dataset (%zu bytes)\n",
                orig_total * sizeof(float));
        return EXIT_FAILURE;
    }
    
    // Open and read the input binary file.
    FILE *fin = fopen(input_file, "rb");
    if (fin == NULL) {
        fprintf(stderr, "Error: failed to open input file '%s'\n", input_file);
        free(orig_data);
        return EXIT_FAILURE;
    }
    
    size_t items_read = fread(orig_data, sizeof(float), orig_total, fin);
    if (items_read != orig_total) {
        fprintf(stderr, "Error: expected %zu float elements, but read %zu\n", orig_total, items_read);
        free(orig_data);
        fclose(fin);
        return EXIT_FAILURE;
    }
    fclose(fin);
    
    // Allocate memory for the expanded dataset.
    float *expanded_data = (float *)malloc(new_total * sizeof(float));
    if (expanded_data == NULL) {
        fprintf(stderr, "Error: failed to allocate memory for expanded dataset (%zu bytes)\n",
                new_total * sizeof(float));
        free(orig_data);
        return EXIT_FAILURE;
    }
    
    // For clarity, we store the new dimensions in an array.
    size_t new_dims[3] = { new_x, new_y, new_z };
    
    // For each voxel in the original dataset, copy its value into an r x r x r block in the new array.
    // The original dataset is stored in row-major order:
    //   index = (z * orig_y + y) * orig_x + x.
    // The new dataset is also stored in row-major order with dimensions new_x, new_y, new_z:
    //   new_index = (new_z_index * new_y + new_y_index) * new_x + new_x_index.
    for (size_t z = 0; z < orig_z; z++) {
        for (size_t y = 0; y < orig_y; y++) {
            for (size_t x = 0; x < orig_x; x++) {
                size_t orig_index = (z * orig_y + y) * orig_x + x;
                float value = orig_data[orig_index];
                
                // Compute the starting index in the expanded array.
                size_t new_x0 = x * r;
                size_t new_y0 = y * r;
                size_t new_z0 = z * r;
                
                // Replicate the value into an r x r x r block.
                for (size_t dz = 0; dz < r; dz++) {
                    for (size_t dy = 0; dy < r; dy++) {
                        for (size_t dx = 0; dx < r; dx++) {
                            size_t new_x_index = new_x0 + dx;
                            size_t new_y_index = new_y0 + dy;
                            size_t new_z_index = new_z0 + dz;
                            size_t new_index = (new_z_index * new_dims[1] + new_y_index) * new_dims[0] + new_x_index;
                            expanded_data[new_index] = value;
                        }
                    }
                }
            }
        }
    }
    
    // Write the expanded dataset to the output file as a binary file of float values.
    FILE *fout = fopen(output_file, "wb");
    if (fout == NULL) {
        fprintf(stderr, "Error: failed to open output file '%s' for writing\n", output_file);
        free(orig_data);
        free(expanded_data);
        return EXIT_FAILURE;
    }
    
    size_t items_written = fwrite(expanded_data, sizeof(float), new_total, fout);
    if (items_written != new_total) {
        fprintf(stderr, "Error: expected to write %zu float elements, but wrote %zu\n", new_total, items_written);
        free(orig_data);
        free(expanded_data);
        fclose(fout);
        return EXIT_FAILURE;
    }
    fclose(fout);
    
    // Cleanup.
    free(orig_data);
    free(expanded_data);
    
    printf("Successfully wrote the expanded dataset to '%s'\n", output_file);
    return EXIT_SUCCESS;
}

