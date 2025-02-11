import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

data = {
    "Dataset": ["QMCPack", "ISABEL", "CESM-ATM", "EXAFEL"],
    "SZ2": [14.35, 58.3, 30.75, 23.4],
    "ZFP": [14.60, 44.7, 31.50, 4.73],
    "zstd": [1.20, 1.15, 1.66, 1.96],
    "C-Blosc2": [1.01, 1.00, 1.40, 1.11],
    "fpzip": [1.75, 2.11, 2.32, 1.11],
    "FPC": [1.09, 1.10, 1.54, 1.00],
}

df = pd.DataFrame(data)

# Melt to long form for Seaborn
df_melted = df.melt(
    id_vars=["Dataset"],
    var_name="Compressor",
    value_name="Compression Ratio",
)

plt.figure(figsize=(8, 4))
sns.set(style="whitegrid", context="talk", font_scale=1.2, font="Times New Roman")

barplot = sns.barplot(
    data=df_melted,
    x="Dataset",
    y="Compression Ratio",
    hue="Compressor",
    palette="colorblind",
)

# Post-process each bar.  If it’s SZ2 or ZFP, apply a hatch pattern:
for patch, compressor in zip(barplot.patches, df_melted["Compressor"]):
    if compressor in ["SZ2", "ZFP"]:
        patch.set_hatch("//")  # Diagonal lines

plt.ylim(0, 60)
plt.tight_layout()

# Add a second legend entry explaining the hatch (optional)
handles, labels = barplot.get_legend_handles_labels()
hatch_patch = mpatches.Patch(
    facecolor="white", hatch="//", edgecolor="black", label="Lossy (SZ2/ZFP)"
)
non_hatch_patch = mpatches.Patch(facecolor="white", edgecolor="black", label="Lossless")
handles.append(hatch_patch)
labels.append("EBLC")


plt.legend(
    handles,
    labels,
    title="Compressor",
    loc="upper left",
    bbox_to_anchor=(1, 1.1),
    frameon=False,
)
plt.savefig("lossless_compression.pdf", bbox_inches="tight")
