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
    "io_results.csv",
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

df["REL Error Bound"] = df["REL Error Bound"].apply(lambda x: f"{x:.0E}")
df["Dataset"] = df["Dataset"].apply(lambda x: x.upper())
df = df[df["REL Error Bound"] != "1E-06"]


# Set the style for all plots
sns.set(style="whitegrid", context="talk", font_scale=1.2, font="Times New Roman")

# Get unique combinations of I/O Method and Dataset
chips = df["Chip"].unique()
datasets = df["Dataset"].unique()
error_bounds = sorted(df["REL Error Bound"].unique())
io_methods = df["I/O Method"].unique()
print(df["Compressor"].unique())

# Create a separate plot for each combination
for io_method in io_methods:
    for chip in chips:
        for dataset in datasets:
            plt.figure(figsize=(5, 5))

            # Filter data for current I/O Method and Dataset
            data = df[
                (df["Chip"] == chip)
                & (df["Dataset"] == dataset)
                & (df["I/O Method"] == io_method)
            ]

            # Calculate mean of uncompressed energy
            uncompressed_mean = data[data["Compressor"] == "Uncompressed"][
                "I/O Energy (J)"
            ].mean()

            # Create the bar plot with reversed x-axis
            ax = sns.barplot(
                x="REL Error Bound",
                y="I/O Energy (J)",
                hue="Compressor",
                hue_order=["SZ2", "SZ3", "ZFP", "QoZ", "SZx"],
                data=data,
            )
            ax.set_xticklabels(error_bounds, rotation=30, ha="right")

            # Add dashed line for uncompressed mean
            ax.axhline(
                y=uncompressed_mean,
                color="r",
                linestyle="--",
                label="Uncompressed Mean",
            )

            # Customize the plot
            plt.title(f"{io_method} - {chip}")
            plt.xlabel("REL Error Bound")
            plt.ylabel("I/O Energy (J)")
            plt.ylim(0, uncompressed_mean * 1.1)

            # Create a legend outside the plot
            handles, labels = ax.get_legend_handles_labels()

            # Add the uncompressed line to the legend
            handles.append(plt.Line2D([0], [0], color="r", linestyle="--"))
            labels.append("Uncompressed Mean")

            # Place the legend outside the plot
            plt.legend(handles, labels, loc="center left", bbox_to_anchor=(1, 0.5))

            # Remove the legend
            ax.get_legend().remove()

            # Adjust layout and save the plot
            plt.tight_layout()
            plt.savefig(f"io_results_{io_method}_{chip}_{dataset}.pdf")
            plt.close()
# plt.legend(loc="center left", bbox_to_anchor=(1, 0.5))
# plt.savefig("io_results.pdf", bbox_inches="tight")

print("Stacked bar plots have been created and saved.")
