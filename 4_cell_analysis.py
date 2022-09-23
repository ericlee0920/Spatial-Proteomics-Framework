import argparse
from pathlib import Path
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def normalize_expression(df, columns, method="arcsinh"):
    def log_transformation(arr):
        return np.log(arr + 1)

    def arcsinh_transformation(arr, cofactor=0.8):
        # unbounded transformation
        arr = np.arcsinh(arr.T / cofactor)
        return arr.T

    def clip_percentile(arr, percentile=99.5):
        # clip the highest value to be the value of the percentile
        cap = np.percentile(arr, percentile)
        return np.where(arr > cap, cap, arr)

    if method == "arcsinh":
        for col in columns:
            temp = arcsinh_transformation(df[col], cofactor=0.8)
            temp = clip_percentile(temp)
            if temp.sum() == 0.0:
                continue
            df[col] = temp
    else:   # log transform
        for col in columns:
            df[col] = log_transformation(df[col])
    return 0


def plot_cell_area(cell_seg, roi, out_dir):
    fig, ax = plt.subplots(1, 1, figsize=(3.8, 2.5), dpi=500)
    # sns.violinplot(x="seg_method", y="cell_area", data=cell_seg, ax=ax, linewidth=0.5, color="white")
    # sns.stripplot(cell_seg["seg_method"], cell_seg["cell_area"], jitter=True, size=2, ax=ax)
    # plt.xlabel("Segmentation Method", fontsize=6)
    # plt.ylabel("Cell Area", fontsize=6)
    # plt.grid(alpha=0.1)
    sns.histplot(data=cell_seg, x="cell_area", binwidth=3)
    plt.xlabel("Cell Area", fontsize=6)
    plt.ylabel("Count", fontsize=6)
    plt.xticks(fontsize=6, rotation=35, ha="right")
    plt.yticks(fontsize=6)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.2, box.width, box.height * 0.8])
    plt.title(f"Cell areas ROI{roi}", weight='bold', fontsize=8)
    plt.savefig(out_dir / f"{roi}_cell_area_hist.png", dpi=500)
    plt.clf()


def plot_cell_stats(cell_df, roi, out_dir):
    fig, ax = plt.subplots(1, 1, figsize=(4.5, 2.5), dpi=500)
    sns.stripplot(cell_df["Marker"], cell_df["MeanIntensity"], jitter=True, size=1, ax=ax)
    # sns.violinplot(x="Marker", y="MeanIntensity", data=cell_df, ax=ax, linewidth=1, color="white")
    plt.xlabel("Marker", fontsize=6)
    plt.ylabel("Intensity", fontsize=6)
    plt.xticks(fontsize=5, rotation=35, ha="right")
    plt.yticks(fontsize=6)
    plt.grid(alpha=0.1)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.2, box.width, box.height * 0.8])
    plt.title(f"Marker intensity ROI{roi}", weight='bold', fontsize=8)
    plt.savefig(out_dir / f"{roi}_cell_marker_intensity.png", dpi=500)
    plt.clf()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Visualization of results from analysis.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-d', '--dirs', nargs='+', default=[], help="Directories to process", required=True)
    args = parser.parse_args()
    data_dirs = args.dirs

    # INPUTS example
    # data_dirs = ["test2_analysis"]

    """
    Visualization
    """
    for str_dir in data_dirs:
        # get path
        path_dir = Path(str_dir)

        # output dirs
        out_dir = path_dir / "viz"
        os.makedirs(out_dir, exist_ok=True)

        # get data
        data = pd.read_csv(path_dir / "cpout" / "cell.csv")
        channel_map = pd.read_csv(path_dir / "cpout" / "panel.csv")
        channel_map = channel_map[channel_map.full == 1].Target.values


        # plot
        plot = True
        range_of_runs = pd.unique(data.ImageNumber)
        roi = range_of_runs[0]

        for roi in range_of_runs:
            df = data[data["ImageNumber"] == roi]

            """
            Fig 1: Method : Cell Area
            """
            cell_seg = pd.DataFrame({"cell_area": df["AreaShape_Area"], "seg_method": ["Ilastik"] * len(df)})

            if plot:
                plot_cell_area(cell_seg, roi, out_dir)

            """
            Fig 2: Markers : Cell Mean Intensity
            """
            grep_col = [i for i in df.columns if i.startswith("Intensity_MeanIntensity_FullStack_")]

            normalize_expression(df, grep_col, method="arcsinh")

            cell_mean = pd.melt(df, id_vars=["ObjectNumber"], value_vars=grep_col, var_name="Marker",
                                value_name="MeanIntensity")
            channel_dict = dict(
                zip(["c" + str(i) for i in range(1, len(channel_map) + 1)], channel_map))
            cell_mean = cell_mean.assign(Marker=[channel_dict[i] for i in [i.split("_")[-1] for i in cell_mean.Marker]])

            if plot:
                plot_cell_stats(cell_mean, roi, out_dir)
