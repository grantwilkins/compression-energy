import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os

import os
import pandas as pd

# Base path for the energy calculation data
base_path = os.path.expanduser(
    "~/compression-energy/analysis/energy-calculation/serial-compress"
)


def load_timechart_data(base_path):
    data = {}

    # Walk through the directory structure
    for dataset in os.listdir(base_path):
        dataset_path = os.path.join(base_path, dataset)
        if os.path.isdir(dataset_path):
            data[dataset] = {}

            for field_compressor in os.listdir(dataset_path):
                field_compressor_path = os.path.join(dataset_path, field_compressor)
                if os.path.isdir(field_compressor_path):
                    # Find the AMDuProf directory
                    amd_prof_dir = next(
                        (
                            d
                            for d in os.listdir(field_compressor_path)
                            if d.startswith("AMDuProf")
                        ),
                        None,
                    )

                    if amd_prof_dir:
                        timechart_path = os.path.join(
                            field_compressor_path, amd_prof_dir, "timechart.csv"
                        )

                        if os.path.exists(timechart_path):
                            # Load the CSV file, skipping the first 150 rows
                            df = pd.read_csv(timechart_path, skiprows=150)

                            # Store the dataframe in the nested dictionary
                            field = field_compressor.split("-")[0]
                            compressor = field_compressor.split("-")[-1]
                            if field not in data[dataset]:
                                data[dataset][field] = {}
                            data[dataset][field][compressor] = df

    return data


# Load all timechart data
timechart_data = load_timechart_data(base_path)

# Constants
time_width = 0.1

total_times = {}
for dataset, dataset_data in timechart_data.items():
    for field, field_data in dataset_data.items():
        for compressor, df in field_data.items():
            total_times[(dataset, field, compressor)] = (
                df["RecordId"].iloc[-1] * 1e3 * time_width
            )


def calculate_cpu_energy(df, cpu_cores):
    return sum(df[f"core{core}-power"].sum() for core in cpu_cores) * time_width


def plot_power_over_time(df, dataset, field, compressor):
    # Calculate the sum of power for cores 0 to 16
    df["total_power"] = df[[f"core{i}-power" for i in range(17)]].sum(axis=1)

    # Create the plot
    plt.figure(figsize=(12, 6))
    plt.plot(df["RecordId"] * time_width, df["total_power"])

    plt.title(f"Total Power Over Time for {dataset} - {field} - {compressor}")
    plt.xlabel("Time (s)")
    plt.ylabel("Total Power (W)")
    plt.grid(True)

    # Save the plot
    plt.savefig(f"{dataset}_{field}_{compressor}_power_plot.png")
    plt.close()


# Plot for NYX temperature.f32 with SZ compressor
nyx_sz_df = timechart_data["nyx"]["temperature.f32"]["sz"]
plot_power_over_time(nyx_sz_df, "NYX", "temperature.f32", "SZ")


df_comp_metrics = pd.read_csv("compression_metrics.csv")

for dataset, dataset_data in timechart_data.items():
    for field, field_data in dataset_data.items():
        for compressor, df in field_data.items():
            energy = calculate_cpu_energy(df, range(0, 16))
            matching_rows = df_comp_metrics[
                (df_comp_metrics["Dataset"] == dataset)
                & (df_comp_metrics["Field"].str.contains(field, case=False))
                & (df_comp_metrics["Compressor"] == compressor)
            ]
            if not matching_rows.empty:
                for _, row in matching_rows.iterrows():
                    compression_energy = (
                        energy * row["Compression Time (ms)"]
                    ) / total_times[(dataset, field, compressor)]
                    decompression_energy = (
                        energy * row["Decompression Time (ms)"]
                    ) / total_times[(dataset, field, compressor)]

                    # Write energy results to the dataframe
                    df_comp_metrics.loc[row.name, "Compression Energy (J)"] = (
                        compression_energy
                    )
                    df_comp_metrics.loc[row.name, "Decompression Energy (J)"] = (
                        decompression_energy
                    )
            else:
                print(f"No matching data found for {dataset} - {field} - {compressor}")

# Save the updated dataframe
df_comp_metrics.to_csv("compression_metrics_with_energy.csv", index=False)
