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
df["Total Runtime (s)"] = (
    df["Compression Runtime (s)"] + df["Decompression Runtime (s)"]
)


sns.set(style="whitegrid", context="talk", font_scale=1.5, font="Times New Roman")
df = df[(df["Compressor"] != "MGARD") & (df["Compressor"] != "MGARD-X")]

# Define the desired order of compressors
compressor_order = [
    "SZ2",
    "SZ3",
    "ZFP",
    "QoZ",
    "SZx",
]


def create_common_legend():
    fig, ax = plt.subplots(figsize=(2, 2))
    palette = sns.color_palette("colorblind", n_colors=len(compressor_order))
    for compressor, color in zip(compressor_order, palette):
        ax.bar(0, 0, color=color, label=compressor)
    handles, labels = ax.get_legend_handles_labels()
    plt.close(fig)
    return handles, labels


# Create common legend
common_handles, common_labels = create_common_legend()


# Function to create and save a stacked bar plot
def create_stacked_energy_plot(data, chip, dataset, common_handles, common_labels):
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
        ax.legend(
            common_handles,
            common_labels,
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


def create_stacked_runtime_plot(data, chip, dataset, common_handles, common_labels):
    # Group by REL Error Bound and Compressor, then calculate mean energies
    grouped = (
        data.groupby(["REL Error Bound", "Compressor"])
        .agg({"Compression Runtime (s)": "mean", "Decompression Runtime (s)": "mean"})
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
                comp_energies.append(err_data["Compression Runtime (s)"].values[0])
                decomp_energies.append(err_data["Decompression Runtime (s)"].values[0])
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
    ax.set_ylabel("Runtime (s)")
    ax.set_xlabel("REL Error Bound")
    ax.set_title(f"{chip}")
    ax.set_xticks(x + width * (len(compressors) - 1) / 2)
    ax.set_xticklabels(error_bounds, rotation=45, ha="right")

    # Create a more succinct legend
    if dataset == "S3D":
        ax.legend(
            common_handles,
            common_labels,
            title="Compressor",
            bbox_to_anchor=(1.05, 1),
            loc="upper left",
            borderaxespad=0.0,
        )

    # ax.set_yscale("log")

    # Set y-axis limits
    data_max = data["Total Runtime (s)"].max()
    data_min = data["Total Runtime (s)"].min()
    # ax.set_ylim(data_min / 10, data_max * 10)
    ax.set_ylim(0, data_max * 1.1)

    plt.tight_layout()
    plt.savefig(f"stacked_runtime_plot_{chip}_{dataset}.pdf", bbox_inches="tight")
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
        create_stacked_energy_plot(subset, chip, dataset, common_handles, common_labels)
        create_stacked_runtime_plot(
            subset, chip, dataset, common_handles, common_labels
        )

print("Stacked bar plots have been created and saved.")

# Calculate and print mean metrics for each compressor, dataset, chip, and error bound
print("\nMean Metrics for each Compressor, Dataset, Chip, and Error Bound:")
print("=" * 100)

metrics = [
    "Compression Ratio",
    "PSNR",
    "Compression Runtime (s)",
    "Decompression Runtime (s)",
    "Bit Rate",
]

df = df[df["Chip"] == "Intel Xeon CPU Max 9480"]

df_s3d = df[df["Dataset"] == "NYX"]
# Create a scatter plot of Compression Ratio vs Total Energy
plt.figure(figsize=(10, 5))


sns.scatterplot(
    x="Compression Ratio",
    y="Total Energy (J)",
    style="Compressor",
    hue="Compressor",
    data=df_s3d,
    markers=["o", "s", "D", "^", "v"],
    s=200,
)
plt.xlabel("Compression Ratio")
plt.ylabel("Total Energy (J)")
plt.legend(
    bbox_to_anchor=(1.02, 0.5), loc="center left", title="Compressor", frameon=False
)
plt.yscale("log")
plt.xscale("log")
plt.grid(True, which="both", ls="--", alpha=0.5)
plt.ylim(1e2, 1e4)

plt.tight_layout()
plt.savefig("compression_ratio_vs_energy.pdf", bbox_inches="tight")
plt.close()

print("Compression Ratio vs Total Energy plot has been created and saved.")

# Calculate correlation between Compression Ratio and Total Energy
correlation = df["Compression Ratio"].corr(df["Total Energy (J)"])
print(f"\nCorrelation between Compression Ratio and Total Energy: {correlation:.4f}")

plt.clf()
plt.figure(figsize=(10, 5))


sns.lineplot(
    x="PSNR",
    y="Total Energy (J)",
    style="Compressor",
    hue="Compressor",
    data=df_s3d,
    markers=["o", "s", "D", "^", "v"],
    # s=200,
)
plt.xlabel("PSNR")
plt.ylabel("Total Energy (J)")
plt.legend(
    bbox_to_anchor=(1.02, 0.5), loc="center left", title="Compressor", frameon=False
)
# plt.yscale("log")
# plt.xscale("log")
# plt.xlim(10, 500)
plt.grid(True, which="both", ls="--", alpha=0.5)

plt.tight_layout()
plt.savefig("psnr_vs_energy.pdf", bbox_inches="tight")
plt.close()

print("Compression Ratio vs Total Energy plot has been created and saved.")

# Calculate correlation between Compression Ratio and Total Energy
correlation = df["Compression Ratio"].corr(df["Total Energy (J)"])
print(f"\nCorrelation between Compression Ratio and Total Energy: {correlation:.4f}")


for dataset in df["Dataset"].unique():
    for error_bound in df["REL Error Bound"].unique():
        print(
            f"\nChip: Intel Xeon CPU Max 9480, Dataset: {dataset}, Error Bound: {error_bound}"
        )
        print("-" * 80)
        for compressor in df["Compressor"].unique():
            subset = df[
                (df["Dataset"] == dataset)
                & (df["Compressor"] == compressor)
                & (df["REL Error Bound"] == error_bound)
            ]
            if not subset.empty:
                print(f"  {compressor}:")
                for metric in metrics:
                    mean_value = subset[metric].mean()
                    print(f"    {metric}: {mean_value:.4f}")
            else:
                print(f"  {compressor}: No data available")
        print()

print("Mean metrics calculation complete.")
