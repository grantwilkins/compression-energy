import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os

# Base path for the energy calculation data
base_path = os.path.expanduser(
    "~/compression-energy/analysis/energy-calculation/serial-compress"
)


def load_timechart_data(base_path):
    data = {}
    for dataset in os.listdir(base_path):
        dataset_path = os.path.join(base_path, dataset)
        if os.path.isdir(dataset_path):
            data[dataset] = {}
            for field_compressor in os.listdir(dataset_path):
                field_compressor_path = os.path.join(dataset_path, field_compressor)
                if os.path.isdir(field_compressor_path):
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
                            df = pd.read_csv(timechart_path, skiprows=150)
                            df["Timestamp"] = pd.to_datetime(
                                df["Timestamp"], format="%H:%M:%S:%f"
                            )
                            df["TotalPower"] = df[
                                [f"core{i}-power" for i in range(128)]
                            ].sum(axis=1)
                            field = field_compressor.split("-")[0]
                            compressor = field_compressor.split("-")[-1]
                            if field not in data[dataset]:
                                data[dataset][field] = {}
                            data[dataset][field][compressor] = df
    return data


# Load all timechart data
timechart_data = load_timechart_data(base_path)


def calculate_energy(start_time, duration_ms, power_data):
    end_time = start_time + pd.Timedelta(milliseconds=duration_ms)
    relevant_power = power_data[
        (power_data["Timestamp"] >= start_time) & (power_data["Timestamp"] < end_time)
    ]

    if len(relevant_power) == 0:
        return 0

    # Calculate energy: Power * Time
    energy = np.trapz(
        relevant_power["TotalPower"], relevant_power["Timestamp"].astype(int) / 10**9
    )
    return energy


def plot_power_over_time(df, dataset, field, compressor):
    plt.figure(figsize=(12, 6))
    plt.plot(df["Timestamp"], df["TotalPower"])
    plt.title(f"Total Power Over Time for {dataset} - {field} - {compressor}")
    plt.xlabel("Time")
    plt.ylabel("Total Power (W)")
    plt.grid(True)
    plt.savefig(f"{dataset}_{field}_{compressor}_power_plot.pdf")
    plt.close()


df_comp_metrics = pd.read_csv("compression_metrics.csv")

for dataset, dataset_data in timechart_data.items():
    for field, field_data in dataset_data.items():
        for compressor, power_data in field_data.items():
            matching_rows = df_comp_metrics[
                (df_comp_metrics["Dataset"] == dataset)
                & (df_comp_metrics["Field"].str.contains(field, case=False))
                & (df_comp_metrics["Compressor"] == compressor)
            ]

            if not matching_rows.empty:
                start_time = power_data["Timestamp"].iloc[0]
                cumulative_time = 0
                for _, row in matching_rows.iterrows():
                    compression_start = start_time + pd.Timedelta(
                        milliseconds=cumulative_time
                    )
                    compression_end = compression_start + pd.Timedelta(
                        milliseconds=row["Compression Time (ms)"]
                    )
                    decompression_start = compression_end
                    decompression_end = decompression_start + pd.Timedelta(
                        milliseconds=row["Decompression Time (ms)"]
                    )

                    compression_energy = calculate_energy(
                        compression_start, row["Compression Time (ms)"], power_data
                    )
                    decompression_energy = calculate_energy(
                        decompression_start, row["Decompression Time (ms)"], power_data
                    )

                    df_comp_metrics.loc[row.name, "Compression Energy (J)"] = (
                        compression_energy
                    )
                    df_comp_metrics.loc[row.name, "Decompression Energy (J)"] = (
                        decompression_energy
                    )

                    cumulative_time += (
                        row["Compression Time (ms)"] + row["Decompression Time (ms)"]
                    )

                # Plot power over time
                plot_power_over_time(power_data, dataset, field, compressor)
            else:
                print(f"No matching data found for {dataset} - {field} - {compressor}")

# Save the updated dataframe
df_comp_metrics.to_csv("compression_metrics_with_energy.csv", index=False)

# Plot for NYX temperature.f32 with SZ compressor (as in your original script)
if (
    "nyx" in timechart_data
    and "temperature.f32" in timechart_data["nyx"]
    and "sz" in timechart_data["nyx"]["temperature.f32"]
):
    nyx_sz_df = timechart_data["nyx"]["temperature.f32"]["sz"]
    plot_power_over_time(nyx_sz_df, "NYX", "temperature.f32", "SZ")
