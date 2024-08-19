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

df = pd.read_csv("compression_metrics_omp.csv")
df["Total Energy (J)"] = df["Compression Energy (J)"] + df["Decompression Energy (J)"]
df["Compressor"] = df["Compressor"].apply(lambda x: x.upper())
df["Dataset"] = df["Dataset"].apply(lambda x: x.upper())
df["Energy per Bit (J/bit)"] = df["Total Energy (J)"] / (df["Number of Elements"] * 32)
# print(df)

sns.set(style="whitegrid", context="talk", font_scale=1.5)
sns.set_palette("colorblind")
df = df[df["Compressor"] != "SZ_OMP"]

# plt.figure(figsize=(10, 6))
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
# plt.legend(loc="center left", bbox_to_anchor=(1, 0.5))
# # plt.clf()
# plt.show()

plt.figure(figsize=(12, 6))  # Increased figure width to accommodate legend
g = sns.catplot(
    x="Threads",
    y="Total Energy (J)",
    hue="Compressor",
    col="Dataset",
    data=df,
    kind="bar",
    sharey=False,
)

# Rotate x-axis labels for each subplot
for ax in g.axes.flat:
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
