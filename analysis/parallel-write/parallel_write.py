import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import math

df = pd.read_csv("parallel_write_experiment.csv")

compressor_order = ["SZ2", "SZ3", "ZFP", "QoZ", "Original"]


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
sns.set(style="whitegrid", context="talk", font_scale=1.2, font="Times New Roman")

print(df["Compressor"].unique())


# Function to create and save a stacked bar plot
def create_stacked_energy_plot(data, dataset):
    # Group by REL Error Bound and Compressor, then calculate mean energies
    grouped = (
        data.groupby(["Cores", "Compressor"])
        .agg({"Compression Energy (J)": "mean", "I/O Energy (J)": "mean"})
        .reset_index()
    )
    print(grouped)

    # Get unique error bounds and compressors
    cores = sorted(grouped["Cores"].unique())
    compressors = [c for c in compressor_order if c in grouped["Compressor"].unique()]

    # Set up the plot
    x = np.arange(len(cores))
    width = 0.8 / len(compressors)
    fig, ax = plt.subplots(figsize=(10, 4.5))

    # Use a colorblind-friendly palette
    palette = sns.color_palette("colorblind", n_colors=len(compressor_order))
    color_dict = dict(zip(compressor_order, palette))

    # Plot bars for each compressor
    for i, compressor in enumerate(compressors):
        comp_data = grouped[grouped["Compressor"] == compressor]

        # Ensure data aligns with error bounds
        comp_energies = []
        io_energies = []
        for core in cores:
            core_data = comp_data[comp_data["Cores"] == core]
            if not core_data.empty:
                num_nodes = math.ceil(core_data["Cores"].values[0] / 48)
                if compressor == "Original":
                    comp_energies.append(0)
                    io_energies.append(
                        core_data["I/O Energy (J)"].values[0] * num_nodes
                    )
                else:
                    comp_energies.append(
                        core_data["Compression Energy (J)"].values[0] * num_nodes
                    )
                    io_energies.append(
                        core_data["I/O Energy (J)"].values[0] * num_nodes
                    )
            else:
                comp_energies.append(0)
                io_energies.append(0)

        ax.bar(
            x + i * width,
            comp_energies,
            width,
            label=compressor,
            color=color_dict[compressor],
            alpha=0.8,
        )
        ax.bar(
            x + i * width,
            io_energies,
            width,
            bottom=comp_energies,
            color=color_dict[compressor],
            alpha=1,
        )

    # Customize the plot
    ax.set_ylabel("Energy (J)")
    ax.set_xlabel("Cores")

    # ax.set_yscale("log")

    # # Set y-axis limits
    # data_max = data["Total Energy (J)"].max()
    # data_min = data["Total Energy (J)"].min()
    # # ax.set_ylim(data_min / 10, data_max * 10)
    # ax.set_ylim(0, data_max * 1.1)
    ax.legend(
        common_handles,
        common_labels,
        title="Compressor",
        bbox_to_anchor=(1.02, 0.5),
        loc="center left",
        borderaxespad=0.0,
        frameon=False,
    )
    ax.set_xticks(x + width * (len(compressors) - 1) / 2)
    ax.set_xticklabels(cores, rotation=45, ha="right")

    plt.tight_layout()
    plt.savefig(f"parallel_write.pdf", bbox_inches="tight")
    plt.close()


# Iterate over each unique combination of Chip and Dataset
create_stacked_energy_plot(df[df["I/O Method"] == "HDF5"], "nyx")
