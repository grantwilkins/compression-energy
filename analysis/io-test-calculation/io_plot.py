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

plt.figure(figsize=(10, 6))
sns.set(style="whitegrid", context="talk", font_scale=1.5, font="Times New Roman")

sns.catplot(
    kind="bar",
    x="REL Error Bound",
    y="I/O Energy (J)",
    hue="Compressor",
    row="I/O Method",
    col="Dataset",
    hue_order=["SZ2", "SZ3", "ZFP", "QoZ"],
    sharey=False,
    data=df,
)
# plt.legend(loc="center left", bbox_to_anchor=(1, 0.5))
plt.savefig("io_results.pdf", bbox_inches="tight")

print("Stacked bar plots have been created and saved.")
