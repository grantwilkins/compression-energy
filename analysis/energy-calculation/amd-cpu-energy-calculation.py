import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os

# Base path for the energy calculation data
base_path = os.path.expanduser(
    "~/compression-energy/analysis/data/abs-serial/serial-compress"
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
                if len(parts) >= 3:
                    field = "-".join(parts[:-2])
                    compressor = parts[-2]
                    error_bound = parts[-1]

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
                                data[dataset][field][compressor][error_bound] = df

    return data


# Load all timechart data
timechart_data = load_timechart_data(base_path)

# Constants
time_width = 0.5  # Updated to 500ms

# Load compression stats
df_comp_metrics = pd.read_csv(
    "/Users/grantwilkins/compression-energy/analysis/data/abs-serial/compression_metrics.csv"
)


def calculate_cpu_energy(df, cpu_cores):
    return sum(df[f"core{core}-power"].sum() for core in cpu_cores) * time_width


def plot_power_over_time(df, dataset, compressor, error_bound):
    # Calculate the sum of power for cores 0 to 16
    df["total_power"] = df[[f"core{i}-power" for i in range(17)]].sum(axis=1)

    # Create the plot
    plt.figure(figsize=(12, 6))
    plt.plot(df["RecordId"] * time_width, df["total_power"])

    plt.title(
        f"Total Power Over Time for {dataset} - {compressor} - Error Bound: {error_bound}"
    )
    plt.xlabel("Time (s)")
    plt.ylabel("Total Power (W)")
    plt.grid(True)

    # Save the plot
    plt.savefig(f"{dataset}_{compressor}_{error_bound}_power_plot.png")
    plt.close()


# # Example plot (you may need to adjust the exact keys based on your data structure)
# example_dataset = list(timechart_data.keys())[0]
# example_compressor = list(timechart_data[example_dataset].keys())[0]
# example_error_bound = list(timechart_data[example_dataset][example_compressor].keys())[
#     0
# ]
# example_df = timechart_data[example_dataset][example_compressor][example_error_bound]
# plot_power_over_time(
#     example_df, example_dataset, example_compressor, example_error_bound
# )

# Calculate and add energy data to compression_stats.csv
for dataset, dataset_data in timechart_data.items():
    for field, field_data in dataset_data.items():
        for compressor, compressor_data in field_data.items():
            for error_bound, df in compressor_data.items():
                energy = calculate_cpu_energy(df, range(0, 127))
                total_time = df["RecordId"].iloc[-1] * time_width

                matching_rows = df_comp_metrics[
                    (df_comp_metrics["Dataset"] == dataset)
                    & (df_comp_metrics["Field"] == field)
                    & (df_comp_metrics["Compressor"] == compressor)
                    & (df_comp_metrics["ABS Error Bound"] == float(error_bound))
                ]

                if not matching_rows.empty:
                    for _, row in matching_rows.iterrows():
                        compression_energy = (
                            energy
                            * row["Compression Runtime (ms)"]
                            / (total_time * 1000)
                        )
                        decompression_energy = (
                            energy
                            * row["Decompression Runtime (ms)"]
                            / (total_time * 1000)
                        )

                        # Write energy results to the dataframe
                        df_comp_metrics.loc[row.name, "Compression Energy (J)"] = (
                            compression_energy
                        )
                        df_comp_metrics.loc[row.name, "Decompression Energy (J)"] = (
                            decompression_energy
                        )
                else:
                    print(
                        f"No matching data found for {dataset} - {field} - {compressor} - Error Bound: {error_bound}"
                    )


# Save the updated dataframe
with open("compression_metrics_with_energy.csv", "w") as f:
    f.write(
        "Compressor,Dataset,Field,REL Error Bound,Iteration,Compression Rate,Decompression Rate,Avg Difference,Avg Error,Difference Range,Error Range,Max Error,Max PW-REL Error,Max REL Error,Min Error,Min PW-REL Error,Min REL Error,MSE,Number of Elements,PSNR,RMSE,NRMSE,Max Value,Mean Value,Min Value,Value Range,Std Dev,Bit Rate,Compressed Size,Compression Ratio,Decompressed Size,Uncompressed Size,Compression Runtime (ms),Decompression Runtime (ms),CPU Core,Compression Energy (J),Decompression Energy (J)\n"
    )
df_comp_metrics.to_csv("compression_metrics_with_energy.csv", index=False)
