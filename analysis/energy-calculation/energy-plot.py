import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
plt.rcParams.update(
    {
        "text.usetex": True,  # Use LaTeX to write all text
        "font.family": "serif",
        "font.serif": ["Linux Libertine O"],  # Specify the Libertine font
        "axes.labelsize": 10,  # LaTeX default is 10pt font.
        "font.size": 12,
        "legend.fontsize": 12,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "text.latex.preamble": r"\usepackage{libertine} \usepackage[libertine]{newtxmath}",  # Load libertine font
    }
)

df = pd.read_csv("compression_metrics_with_energy.csv")
df["Total Energy (J)"] = df["Compression Energy (J)"] + df["Decompression Energy (J)"]

sns.set(style="whitegrid", context="talk", font_scale=1.5)
sns.set_palette("colorblind")

sns.lineplot(
    x="Bit Rate",
    y="Total Energy (J)",
    hue="Dataset",
    marker="o",
    markersize=12,
    data=df,
    # legend=False,
)

plt.show()

sns.barplot(x="Error Bound", y="Total Energy (J)", hue="Dataset", data=df)
plt.show()
