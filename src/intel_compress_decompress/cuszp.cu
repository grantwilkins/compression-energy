#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#include <cuda_runtime.h>
#include <nvml.h>
#include <papi.h>
#include <cuSZp_utility.h>
#include <cuSZp_entry_f32.h>
#include <cuSZp_entry_f64.h>
#include <cuSZp_timer.h>

#define MAX_ITERATIONS 25
#define CONFIDENCE_LEVEL 1.96
#define MAX_POWERCAP_EVENTS 64
#define CHECK_NVML(call) { nvmlReturn_t result = call; if (result != NVML_SUCCESS) { fprintf(stderr, "NVML Error: %s\n", nvmlErrorString(result)); exit(1); } }
#define CHECK_PAPI(call) { int retval = call; if (retval != PAPI_OK) { fprintf(stderr, "PAPI Error: %s\n", PAPI_strerror(retval)); exit(1); } }

typedef struct {
    double compression_time;
    double decompression_time;
    double compression_throughput;
    double decompression_throughput;
    double compression_ratio;
    double max_error;
    double avg_error;
    double mse;
    double psnr;
    double nrmse;
    unsigned long compressed_size;
    double cpu_comp_energy;
    double cpu_decomp_energy;
    unsigned long gpu_comp_energy;
    unsigned long gpu_decomp_energy;
} CompressionMetrics;

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
    if (n < 2) return false;
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

void calculate_error_metrics(void *original, void *decompressed, size_t num_elements, int data_type, CompressionMetrics *metrics) {
    double max_error = 0.0, sum_squared_error = 0.0, sum_error = 0.0;
    double min_val = INFINITY, max_val = -INFINITY;

    if (data_type == 0) { // float
        float *orig = (float *)original;
        float *decomp = (float *)decompressed;
        for (size_t i = 0; i < num_elements; i++) {
            double error = fabs(orig[i] - decomp[i]);
            max_error = fmax(max_error, error);
            sum_squared_error += error * error;
            sum_error += error;
            min_val = fmin(min_val, orig[i]);
            max_val = fmax(max_val, orig[i]);
        }
    } else { // double
        double *orig = (double *)original;
        double *decomp = (double *)decompressed;
        for (size_t i = 0; i < num_elements; i++) {
            double error = fabs(orig[i] - decomp[i]);
            max_error = fmax(max_error, error);
            sum_squared_error += error * error;
            sum_error += error;
            min_val = fmin(min_val, orig[i]);
            max_val = fmax(max_val, orig[i]);
        }
    }

    metrics->max_error = max_error;
    metrics->avg_error = sum_error / num_elements;
    metrics->mse = sum_squared_error / num_elements;
    double value_range = max_val - min_val;
    metrics->psnr = 20 * log10(value_range) - 10 * log10(metrics->mse);
    metrics->nrmse = sqrt(metrics->mse) / value_range;
}

int main(int argc, char *argv[]) {
    if (argc != 4) {
        printf("Usage: %s <dataset_file> <error_mode> <error_bound>\n", argv[0]);
        return 1;
    }

    const char *dataset_file = argv[1];
    const char *error_mode = argv[2];
    double error_bound = atof(argv[3]);

    // Initialize CUDA
    cudaSetDevice(0);

    // Initialize NVML
    nvmlInit();
    nvmlDevice_t device;
    CHECK_NVML(nvmlDeviceGetHandleByIndex(0, &device));

    // Initialize PAPI
    int EventSet = PAPI_NULL;
    long long values[MAX_POWERCAP_EVENTS];
    int num_events = 0;
    char event_names[MAX_POWERCAP_EVENTS][PAPI_MAX_STR_LEN];
    int data_type[MAX_POWERCAP_EVENTS];

    CHECK_PAPI(PAPI_library_init(PAPI_VER_CURRENT));
    CHECK_PAPI(PAPI_create_eventset(&EventSet));

    // Find and add powercap events
    int code = PAPI_NATIVE_MASK;
    PAPI_event_info_t info;
    while (PAPI_enum_cmp_event(&code, PAPI_ENUM_FIRST, 0) == PAPI_OK) {
        if (PAPI_get_event_info(code, &info) == PAPI_OK) {
            if (strstr(info.symbol, "powercap")) {
                if (PAPI_add_event(EventSet, code) == PAPI_OK) {
                    strncpy(event_names[num_events], info.symbol, PAPI_MAX_STR_LEN);
                    data_type[num_events] = info.data_type;
                    num_events++;
                }
            }
        }
    }

    // Read dataset
    size_t num_elements;
    int status = 0;
    void *data = NULL;
    int data_type;

    if (strstr(dataset_file, "nyx") != NULL || strstr(dataset_file, "hacc") != NULL || strstr(dataset_file, "miranda") != NULL) {
        data = (void *)readFloatData_Yafan(dataset_file, &num_elements, &status);
        data_type = 0; // float
    } else if (strstr(dataset_file, "s3d") != NULL) {
        data = (void *)readDoubleData_Yafan(dataset_file, &num_elements, &status);
        data_type = 1; // double
    } else {
        fprintf(stderr, "Unknown dataset %s\n", dataset_file);
        return 1;
    }

    if (status != 0) {
        fprintf(stderr, "Failed to read dataset %s\n", dataset_file);
        return 1;
    }

    // Allocate memory on GPU
    void *d_data, *d_compressed, *d_decompressed;
    size_t data_size = num_elements * (data_type == 0 ? sizeof(float) : sizeof(double));
    cudaMalloc(&d_data, data_size);
    cudaMalloc(&d_compressed, data_size);
    cudaMalloc(&d_decompressed, data_size);
    cudaMemcpy(d_data, data, data_size, cudaMemcpyHostToDevice);

    // Prepare for compression
    cudaStream_t stream;
    cudaStreamCreate(&stream);

    // Metrics arrays
    double compression_times[MAX_ITERATIONS];
    double decompression_times[MAX_ITERATIONS];
    CompressionMetrics metrics[MAX_ITERATIONS];

    int iteration = 0;
    bool confidence_interval_reached = false;

    while (iteration < MAX_ITERATIONS && !confidence_interval_reached) {
        // Compression
        size_t compressed_size;
        
        // Start CPU and GPU energy measurement for compression
        unsigned long long gpu_energy_start, gpu_energy_end;
        CHECK_PAPI(PAPI_start(EventSet));
        CHECK_NVML(nvmlDeviceGetTotalEnergyConsumption(device, &gpu_energy_start));

        double start_time = get_time();
        if (data_type == 0) {
            SZp_compress_deviceptr_f32((float *)d_data, (unsigned char *)d_compressed, num_elements, &compressed_size, error_bound, stream);
        } else {
            SZp_compress_deviceptr_f64((double *)d_data, (unsigned char *)d_compressed, num_elements, &compressed_size, error_bound, stream);
        }
        cudaStreamSynchronize(stream);
        double end_time = get_time();
        
        CHECK_NVML(nvmlDeviceGetTotalEnergyConsumption(device, &gpu_energy_end));
        CHECK_PAPI(PAPI_stop(EventSet, values));

        compression_times[iteration] = end_time - start_time;
        metrics[iteration].gpu_comp_energy = gpu_energy_end - gpu_energy_start;

        // Calculate CPU energy for compression
        metrics[iteration].cpu_comp_energy = 0.0;
        for (int i = 0; i < num_events; i++) {
            if (strstr(event_names[i], "ENERGY_UJ") && data_type[i] == PAPI_DATATYPE_UINT64) {
                metrics[iteration].cpu_comp_energy += values[i] / 1.0e6;
            }
        }

        // Decompression
        // Start CPU and GPU energy measurement for decompression
        CHECK_PAPI(PAPI_start(EventSet));
        CHECK_NVML(nvmlDeviceGetTotalEnergyConsumption(device, &gpu_energy_start));

        start_time = get_time();
        if (data_type == 0) {
            SZp_decompress_deviceptr_f32((float *)d_decompressed, (unsigned char *)d_compressed, num_elements, compressed_size, error_bound, stream);
        } else {
            SZp_decompress_deviceptr_f64((double *)d_decompressed, (unsigned char *)d_compressed, num_elements, compressed_size, error_bound, stream);
        }
        cudaStreamSynchronize(stream);
        end_time = get_time();

        CHECK_NVML(nvmlDeviceGetTotalEnergyConsumption(device, &gpu_energy_end));
        CHECK_PAPI(PAPI_stop(EventSet, values));

        decompression_times[iteration] = end_time - start_time;
        metrics[iteration].gpu_decomp_energy = gpu_energy_end - gpu_energy_start;

        // Calculate CPU energy for decompression
        metrics[iteration].cpu_decomp_energy = 0.0;
        for (int i = 0; i < num_events; i++) {
            if (strstr(event_names[i], "ENERGY_UJ") && data_type[i] == PAPI_DATATYPE_UINT64) {
                metrics[iteration].cpu_decomp_energy += values[i] / 1.0e6;
            }
        }

        // Calculate metrics
        void *h_decompressed = malloc(data_size);
        cudaMemcpy(h_decompressed, d_decompressed, data_size, cudaMemcpyDeviceToHost);

        calculate_error_metrics(data, h_decompressed, num_elements, data_type, &metrics[iteration]);

        metrics[iteration].compression_time = compression_times[iteration];
        metrics[iteration].decompression_time = decompression_times[iteration];
        metrics[iteration].compression_throughput = (data_size / 1e9) / compression_times[iteration];
        metrics[iteration].decompression_throughput = (data_size / 1e9) / decompression_times[iteration];
        metrics[iteration].compression_ratio = (double)data_size / compressed_size;
        metrics[iteration].compressed_size = compressed_size;

        free(h_decompressed);

        // Write metrics to CSV file
        FILE *csv_file = fopen("cuszip_compression_metrics.csv", "a");
        if (csv_file == NULL) {
            fprintf(stderr, "Error opening CSV file\n");
        } else {
            fprintf(csv_file, "cuSZp,%s,%s,%e,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%lu,%f,%f,%lu,%lu\n",
                    dataset_file, error_mode, error_bound, iteration,
                    metrics[iteration].compression_time,
                    metrics[iteration].decompression_time,
                    metrics[iteration].compression_throughput,
                    metrics[iteration].decompression_throughput,
                    metrics[iteration].compression_ratio,
                    metrics[iteration].max_error,
                    metrics[iteration].avg_error,
                    metrics[iteration].mse,
                    metrics[iteration].psnr,
                    metrics[iteration].nrmse,
                    metrics[iteration].compressed_size,
                    metrics[iteration].cpu_comp_energy,
                    metrics[iteration].cpu_decomp_energy,
                    metrics[iteration].gpu_comp_energy,
                    metrics[iteration].gpu_decomp_energy);
            fclose(csv_file);
        }

        iteration++;

        // Check if we've reached the confidence interval
        if (iteration >= 5) {
            confidence_interval_reached =
                within_confidence_interval(compression_times, iteration) &&
                within_confidence_interval(decompression_times, iteration);
        }
    }
        
    free(data);
    cudaFree(d_data);
    cudaFree(d_compressed);
    cudaFree(d_decompressed);
    cudaStreamDestroy(stream);
    nvmlShutdown();
    PAPI_cleanup_eventset(EventSet);
    PAPI_destroy_eventset(&EventSet);
    PAPI_shutdown();

}

