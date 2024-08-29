import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
# plt.rcParams.update(
#     {
#         "text.usetex": False,  # Use LaTeX to write all text
#         "font.family": "Times New Roman",
#         "font.serif": ["Times New Roman"],  # Specify Times New Roman
#         "axes.labelsize": 10,  # LaTeX default is 10pt font.
#         "font.size": 12,
#         "legend.fontsize": 12,
#         "xtick.labelsize": 10,
#         "ytick.labelsize": 10,
#         "text.latex.preamble": r"\usepackage{times}",  # Load Times font for LaTeX
#     }
# )
import pandas as pd

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
# For example, if 'Compressor' column is missing in some files:
if "Compressor" not in df.columns:
    df["Compressor"] = "Unknown"

# Reset the index of the combined dataframe
df = df.reset_index(drop=True)

# Print the shape of the resulting dataframe
print(f"Combined dataframe shape: {df.shape}")

# Optionally, print the column names to verify all expected columns are present
print("Columns in the combined dataframe:")
print(df.columns.tolist())

# df = pd.read_csv("compression_metrics_szx_serial.csv")
df["Total Energy (J)"] = df["Compression Energy (J)"] + df["Decompression Energy (J)"]
df["REL Error Bound"] = df["REL Error Bound"].apply(lambda x: f"{x:.0E}")
df["Dataset"] = df["Dataset"].apply(lambda x: x.upper())
df["Energy per Bit (J/bit)"] = df["Total Energy (J)"] / (df["Number of Elements"] * 32)
# print(df)

sns.set(style="whitegrid", context="talk", font_scale=1.5)
# sns.set_palette("colorblind")
df = df[df["Compressor"] != "MGARD"]

# sns.relplot(
#     kind="line",
#     x="Bit Rate",
#     y="Energy per Bit (J/bit)",
#     hue="Compressor",
#     col="Dataset",
#     marker="o",
#     markersize=12,
#     data=df,
# )
# plt.legend(loc="center left", bbox_to_anchor=(1.1, 0.5))
# # plt.clf()
# plt.show()

# plt.clf()
g = sns.catplot(
    x="REL Error Bound",
    y="Total Energy (J)",
    hue="Compressor",
    col="Dataset",
    data=df,
    kind="bar",
    col_wrap=2,
)

plt.legend(loc="center left", bbox_to_anchor=(1.1, 0.5))

# Rotate x-axis labels for each subplot
for ax in g.axes.flat:
    ax.set_ylim(1e2, 1e6)
    ax.set_yscale("log")
    ax.set_xticklabels(
        ["1e-1", "1e-2", "1e-3", "1e-4", "1e-5", "1e-6"],
        rotation=45,
    )

plt.tight_layout()
plt.show()
plt.clf()

g = sns.catplot(
    x="REL Error Bound",
    y="Compression Runtime (s)",
    hue="Compressor",
    col="Dataset",
    data=df,
    kind="point",
    sharey=False,
    col_wrap=2,
)

plt.legend(loc="center left", bbox_to_anchor=(1.1, 0.5))

# Rotate x-axis labels for each subplot
# s

plt.tight_layout()
plt.show()

plt.clf()
g = sns.catplot(
    x="REL Error Bound",
    y="Compression Ratio",
    hue="Compressor",
    col="Dataset",
    data=df,
    kind="bar",
    sharey=False,
    col_wrap=2,
)

plt.legend(loc="center left", bbox_to_anchor=(1.1, 0.5))

# Rotate x-axis labels for each subplot
# s
for ax in g.axes.flat:
    # ax.set_ylim(1e2, 1e6)
    ax.set_yscale("log")
    ax.set_xticklabels(
        ax.get_xticklabels(),
        rotation=45,
    )

plt.tight_layout()
plt.show()

# sns.relplot(
#     kind="line",
#     x="Bit Rate",
#     y="PSNR",
#     hue="Compressor",
#     col="Dataset",
#     # marker="o",
#     # markersize=12,
#     data=df,
#     # legend=False,
# )
# plt.show()
