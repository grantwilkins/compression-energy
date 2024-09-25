import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["font.family"] = "Times New Roman"

# List of CSV files
csv_files = [
    "compression_metrics.csv",
    "compression_metrics_szx_serial.csv",
    "compression_metrics_mgard_x.csv",
]

# Read and combine all CSV files
dfs = []
for file in csv_files:
    try:
        df = pd.read_csv(file)
        dfs.append(df)
    except FileNotFoundError:
        print(f"Warning: File {file} not found. Skipping.")

# Combine all dataframes
df = pd.concat(dfs, ignore_index=True)

# Fill missing values with 0
df = df.fillna(0)

# If there are any columns that should be specifically handled, add logic here
if "Compressor" not in df.columns:
    df["Compressor"] = "Unknown"

# Reset the index of the combined dataframe
df = df.reset_index(drop=True)

df["Total Energy (J)"] = df["Compression Energy (J)"] + df["Decompression Energy (J)"]
df["REL Error Bound"] = df["REL Error Bound"].apply(lambda x: f"{x:.0E}")
df["Dataset"] = df["Dataset"].apply(lambda x: x.upper())
df["Energy per Bit (J/bit)"] = df["Total Energy (J)"] / (df["Number of Elements"] * 32)
df = df[df["REL Error Bound"] != "1E-06"]

sns.set(style="whitegrid", context="talk", font_scale=1.5, font="Times New Roman")
df = df[(df["Compressor"] != "MGARD") & (df["Compressor"] != "MGARD-X")]

# Define the desired order of compressors
compressor_order = ["SZ2", "SZ3", "ZFP", "QoZ", "SZx", "MGARD-X"]


# Function to create and save a stacked bar plot
def create_stacked_energy_plot(data, chip, dataset):
    # Group by REL Error Bound and Compressor, then calculate mean energies
    grouped = (
        data.groupby(["REL Error Bound", "Compressor"])
        .agg({"Compression Energy (J)": "mean", "Decompression Energy (J)": "mean"})
        .reset_index()
    )

    # Get unique error bounds and compressors
    error_bounds = sorted(grouped["REL Error Bound"].unique())
    compressors = [c for c in compressor_order if c in grouped["Compressor"].unique()]

    # Check if there are any compressors in the data
    if not compressors:
        print(f"No compressors found for {chip} - {dataset}. Skipping this plot.")
        return

    plt.figure(figsize=(6, 6))
    if dataset == "S3D":
        plt.figure(figsize=(8, 6))

    # Set up the plot
    x = np.arange(len(error_bounds))
    width = 0.8 / len(compressors)
    fig, ax = plt.subplots(figsize=(6, 6))
    if dataset == "S3D":
        fig, ax = plt.subplots(figsize=(8, 6))

    # Use a colorblind-friendly palette
    palette = sns.color_palette("colorblind", n_colors=len(compressor_order))
    color_dict = dict(zip(compressor_order, palette))

    # Plot bars for each compressor
    for i, compressor in enumerate(compressors):
        comp_data = grouped[grouped["Compressor"] == compressor]

        # Ensure data aligns with error bounds
        comp_energies = []
        decomp_energies = []
        for err in error_bounds:
            err_data = comp_data[comp_data["REL Error Bound"] == err]
            if not err_data.empty:
                comp_energies.append(err_data["Compression Energy (J)"].values[0])
                decomp_energies.append(err_data["Decompression Energy (J)"].values[0])
            else:
                comp_energies.append(0)
                decomp_energies.append(0)

        ax.bar(
            x + i * width,
            comp_energies,
            width,
            label=compressor,
            color=color_dict[compressor],
            alpha=0.7,
        )
        ax.bar(
            x + i * width,
            decomp_energies,
            width,
            bottom=comp_energies,
            color=color_dict[compressor],
            alpha=1,
        )

    # Customize the plot
    ax.set_ylabel("Energy (J)")
    ax.set_xlabel("REL Error Bound")
    ax.set_title(f"{chip}")
    ax.set_xticks(x + width * (len(compressors) - 1) / 2)
    ax.set_xticklabels(error_bounds, rotation=45, ha="right")

    # Create a more succinct legend
    if dataset == "S3D":
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(
            handles,
            labels,
            title="Compressor",
            bbox_to_anchor=(1.05, 1),
            loc="upper left",
            borderaxespad=0.0,
        )

    # ax.set_yscale("log")

    # Set y-axis limits
    data_max = data["Total Energy (J)"].max()
    data_min = data["Total Energy (J)"].min()
    # ax.set_ylim(data_min / 10, data_max * 10)
    ax.set_ylim(0, data_max * 1.1)

    plt.tight_layout()
    plt.savefig(f"stacked_energy_plot_{chip}_{dataset}.pdf", bbox_inches="tight")
    plt.close()


# Iterate over each unique combination of Chip and Dataset
for chip in df["Chip"].unique():
    for dataset in df["Dataset"].unique():
        # Filter data for current chip and dataset
        subset = df[(df["Chip"] == chip) & (df["Dataset"] == dataset)]

        if subset.empty:
            print(f"No data for {chip} - {dataset}. Skipping.")
            continue

        # Create stacked bar plot for Compression and Decompression Energy
        create_stacked_energy_plot(subset, chip, dataset)

print("Stacked bar plots have been created and saved.")
