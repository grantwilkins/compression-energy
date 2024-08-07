import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os

# Base path for the energy calculation data
base_path = os.path.expanduser(
    "~/compression-energy/analysis/data/abs-omp/omp-compress"
)


def load_timechart_data(base_path):
    data = {}

    # Walk through the directory structure
    for dataset in os.listdir(base_path):
        dataset_path = os.path.join(base_path, dataset)
        if os.path.isdir(dataset_path):
            data[dataset] = {}

            for item in os.listdir(dataset_path):
                parts = item.split("-")
                if len(parts) >= 4:
                    field = parts[0]
                    compressor = parts[1]
                    error_bound = parts[2]
                    threads = parts[3]

                    item_path = os.path.join(dataset_path, item)
                    if os.path.isdir(item_path):
                        # Find the AMDuProf directory
                        amd_prof_dir = next(
                            (
                                d
                                for d in os.listdir(item_path)
                                if d.startswith("AMDuProf")
                            ),
                            None,
                        )

                        if amd_prof_dir:
                            timechart_path = os.path.join(
                                item_path, amd_prof_dir, "timechart.csv"
                            )

                            if os.path.exists(timechart_path):
                                # Load the CSV file, skipping the first 150 rows
                                df = pd.read_csv(timechart_path, skiprows=150)

                                # Store the dataframe in the nested dictionary
                                if field not in data[dataset]:
                                    data[dataset][field] = {}
                                if compressor not in data[dataset][field]:
                                    data[dataset][field][compressor] = {}
                                if error_bound not in data[dataset][field][compressor]:
                                    data[dataset][field][compressor][error_bound] = {}
                                data[dataset][field][compressor][error_bound][
                                    threads
                                ] = df

    return data


# Load all timechart data
timechart_data = load_timechart_data(base_path)

# Constants
time_width = 0.5  # 500ms

# Load compression stats
df_comp_metrics = pd.read_csv(
    "/Users/grantwilkins/compression-energy/analysis/data/abs-omp/compression_metrics_omp.csv"
)


def calculate_cpu_energy(df, cpu_cores):
    return sum(df[f"core{core}-power"].sum() for core in cpu_cores) * time_width


# Calculate and add energy data to compression_metrics_omp.csv
for dataset, dataset_data in timechart_data.items():
    for field, field_data in dataset_data.items():
        for compressor, compressor_data in field_data.items():
            for error_bound, error_bound_data in compressor_data.items():
                for threads, df in error_bound_data.items():
                    energy = calculate_cpu_energy(df, range(0, 127))
                    total_time = df["RecordId"].iloc[-1] * time_width

                    matching_rows = df_comp_metrics[
                        (df_comp_metrics["Dataset"] == dataset)
                        & (df_comp_metrics["Field"] == field)
                        & (df_comp_metrics["Compressor"] == compressor)
                        & (df_comp_metrics["ABS Error Bound"] == float(error_bound))
                        & (df_comp_metrics["Threads"] == int(threads))
                    ]

                    if not matching_rows.empty:
                        for _, row in matching_rows.iterrows():
                            compression_energy = (
                                energy * row["Compression Runtime (s)"] / (total_time)
                            )
                            decompression_energy = (
                                energy * row["Decompression Runtime (s)"] / (total_time)
                            )

                            # Write energy results to the dataframe
                            df_comp_metrics.loc[row.name, "Compression Energy (J)"] = (
                                compression_energy
                            )
                            df_comp_metrics.loc[
                                row.name, "Decompression Energy (J)"
                            ] = decompression_energy
                    else:
                        print(
                            f"No matching data found for {dataset} - {field} - {compressor} - Error Bound: {error_bound} - Threads: {threads}"
                        )

# Save the updated dataframe
df_comp_metrics.to_csv("compression_metrics_omp_with_energy.csv", index=False)
