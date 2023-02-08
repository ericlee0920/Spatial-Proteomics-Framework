import argparse
import os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import phenograph
from itertools import combinations


def process_data(cell_df, rel_df, panel_arr):
    """
    Processing output from the modified ImcSegmentationPipeline
    """
    # Subset expression columns
    exp_cols = [i for i in cell_df.columns if i.startswith("Intensity_MeanIntensity_FullStack_")]
    exp_df = cell_df[exp_cols]

    # Renaming according to panel
    cols_dict = dict(zip(["c" + str(i) for i in range(1, len(panel_arr) + 1)], panel_arr))
    exp_df.columns = [cols_dict[i] for i in [i.split("_")[-1] for i in exp_df.columns]]

    # Scaling up and normalize
    exp_df = exp_df * 1E5
    exp_df = normalize_expression(exp_df)
    exp_df = pd.concat([cell_df.iloc[:, 0], exp_df], axis=1)

    # Subset location columns
    loc_cols = ['ImageNumber', 'AreaShape_Area', 'AreaShape_Center_X', 'AreaShape_Center_Y']
    loc_df = cell_df[loc_cols]
    loc_df.columns = ["ImageNumber", "Area", "X", "Y"]

    # Subset relation columns
    rel_cols = ['First Image Number', 'First Object Number', 'Second Object Number']
    rel_df = rel_df[rel_cols]
    rel_df.columns = ["ImageNumber", "Cell", "Neighbour"]

    return exp_df, loc_df, rel_df


def normalize_expression(df, method="arcsinh"):
    """
    Normalize expression data using arcsinh or log transformation
    """
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
        for col in df.columns:
            temp = arcsinh_transformation(df[col], cofactor=0.8)
            temp = clip_percentile(temp)
            if temp.sum() == 0.0:
                continue
            df[col] = temp
    elif method == "log":   # log transform
        for col in df.columns:
            df[col] = log_transformation(df[col])
    else:
        raise Exception("Input arcsinh or log")

    return df


def df_to_csv(exp_df, loc_df, rel_df, run_name, out_dir):
    """
    Output files to csv
    """
    exp_df.to_csv(out_dir / f"{run_name}_exp.csv", index=False)
    loc_df.to_csv(out_dir / f"{run_name}_loc.csv", index=False)
    rel_df.to_csv(out_dir / f"{run_name}_rel.csv", index=False)


def plot_cell_area(cell_area, run_name, out_dir, roi=None):
    fig, ax = plt.subplots(1, 1, figsize=(3.8, 2), dpi=500)
    # sns.violinplot(x="seg_method", y="cell_area", data=cell_area, ax=ax, linewidth=0.5, color="white")
    # sns.stripplot(cell_area["seg_method"], cell_area["cell_area"], jitter=True, size=2, ax=ax)
    sns.histplot(data=cell_area, x="cell_area", binwidth=3)

    plt.xlabel("Cell Area", fontsize=6)
    plt.ylabel("Cell Count", fontsize=6)
    plt.xticks(fontsize=6, rotation=35, ha="right")
    plt.yticks(fontsize=6)

    plt.grid(alpha=0.1)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.2, box.width, box.height * 0.8])

    if roi:
        plt.title(f"{run_name} ROI{roi} cell area histogram", weight='bold', fontsize=8)
        plt.savefig(out_dir / f"{roi}_cell_area_hist.png")
    else:
        plt.title(f"{run_name} cell area histogram", weight='bold', fontsize=8)
        plt.savefig(out_dir / f"{run_name}_cell_area_hist.png")
    plt.clf()


def plot_cell_mean_intensity(cell_mean, run_name, out_dir, roi=None):
    fig, ax = plt.subplots(1, 1, figsize=(3.8, 2), dpi=500)
    # sns.stripplot(cell_mean["Marker"], cell_mean["Mean_Intensity"], jitter=True, size=1, ax=ax)
    sns.violinplot(x="Marker", y="Mean_Intensity", data=cell_mean, linewidth=1, color="white", ax=ax)

    plt.xlabel("Marker", fontsize=6)
    plt.ylabel("Mean Intensity", fontsize=6)
    plt.xticks(fontsize=6, rotation=35, ha="right")
    plt.yticks(fontsize=6)

    plt.grid(alpha=0.1)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.2, box.width, box.height * 0.8])

    if roi:
        plt.title(f"{run_name} ROI{roi} mean marker intensity", weight='bold', fontsize=8)
        plt.savefig(out_dir / f"{roi}_cell_marker_intensity.png")
    else:
        plt.title(f"{run_name} mean marker intensity", weight='bold', fontsize=8)
        plt.savefig(out_dir / f"{run_name}_cell_marker_intensity.png")
    plt.clf()


def plot_cluster_mean_heatmap(cluster_mean, run_name, out_dir):
    fig, ax = plt.subplots(1, 1, figsize=(5, 5), dpi=500)

    sns.heatmap(cluster_mean, cmap="coolwarm", ax=ax)

    plt.xlabel("Marker", fontsize=6)
    plt.ylabel("Label", fontsize=6)
    plt.xticks(fontsize=6)
    plt.yticks(fontsize=6)

    cax = ax.figure.axes[-1]
    cax.tick_params(labelsize=6)

    box = ax.get_position()
    ax.set_position([box.x0 + box.height * 0.025, box.y0 + box.height * 0.1, box.width, box.height - box.height * 0.1])
    plt.title(f"{run_name} mean cluster heatmap", weight='bold', fontsize=8)

    plt.savefig(out_dir / f"{run_name}_mean_cluster_heatmap.png")
    plt.clf()


def plot_cell_count(cell_count, run_name, out_dir):
    fig, ax = plt.subplots(1, 1, figsize=(3.8, 2.5), dpi=500)

    sns.barplot(data=cell_count, x="ROI", y="Count")

    plt.xlabel("ROI", fontsize=6)
    plt.ylabel("Cell Count", fontsize=6)
    plt.xticks(fontsize=6)
    plt.yticks(fontsize=6)

    box = ax.get_position()
    ax.set_position(
        [box.x0 + box.height * 0.025, box.y0 + box.height * 0.1, box.width, box.height - box.height * 0.1])
    plt.title(f"{run_name} cell count", weight='bold', fontsize=8)

    plt.savefig(out_dir / f"{run_name}_mean_cell_count.png")
    plt.clf()


def plot_cell_mean_intensity_zoomin(cell_mean, run_name, out_dir):
    pairs = [(1, 9), (9, 17), (17, 23)]
    for a, b in pairs:
        cell_subset = cell_mean[cell_mean.Marker.isin(exp_df.columns[a:b])]
        fig, ax = plt.subplots(1, 1, figsize=(7, 4), dpi=500)
        sns.stripplot(cell_subset["Marker"], cell_subset["Mean_Intensity"], jitter=True, size=1, ax=ax,
                      alpha=0.08)
        sns.violinplot(x="Marker", y="Mean_Intensity", data=cell_subset, linewidth=1, color="black", ax=ax)

        plt.xlabel("Marker", fontsize=6)
        plt.ylabel("Mean Intensity", fontsize=6)
        plt.xticks(fontsize=6, rotation=35, ha="right")
        plt.yticks(fontsize=6)
        plt.ylim(0, cell_mean["Mean_Intensity"].max())

        plt.grid(alpha=0.1)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        box = ax.get_position()
        ax.set_position([box.x0, box.y0 + box.height * 0.2, box.width, box.height * 0.8])

        plt.title(f"{run_name} mean marker intensity", weight='bold', fontsize=8)
        plt.savefig(out_dir / f"{run_name}_cell_marker_{a}_{b}.png")
        plt.clf()


if __name__ == '__main__':
    # parser = argparse.ArgumentParser(description="Visualization of results from analysis.",
    #                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # parser.add_argument('-d', '--dirs', nargs='+', default=[], help="Directories to process", required=True)
    # args = parser.parse_args()
    # data_dirs = args.dirs

    # INPUTS example
    data_dirs = ["XP4570_analysis", "XP4571_analysis", "XP4572_analysis"]

    # # Step 1: Process the data
    # for str_dir in data_dirs:
    #     run_name = str_dir.split("_")[0]
    #     print(f"Processing {run_name}")
    #
    #     # Get path
    #     path_dir = Path(str_dir)
    #
    #     # Output dirs
    #     out_dir = path_dir / "viz"
    #     os.makedirs(out_dir, exist_ok=True)
    #
    #     # Get data
    #     cell_df = pd.read_csv(path_dir / "cpout" / "cell.csv")
    #     rel_df = pd.read_csv(path_dir / "cpout" / "Object relationships.csv")
    #     panel_df = pd.read_csv(path_dir / "cpout" / "panel.csv")
    #     panel_arr = panel_df[panel_df.full == 1].Target.values
    #
    #     # Normalize data and output to csv
    #     exp_df, loc_df, rel_df = process_data(cell_df, rel_df, panel_arr)
    #     df_to_csv(exp_df, loc_df, rel_df, run_name, out_dir)
    #
    # # Step 2: Cluster the data
    # for str_dir in data_dirs:
    #     run_name = str_dir.split("_")[0]
    #     print(f"Clustering {run_name}")
    #
    #     # Get path
    #     path_dir = Path(str_dir)
    #
    #     # Output dirs
    #     out_dir = path_dir / "viz"
    #
    #     # Get data
    #     exp_df = pd.read_csv(out_dir / f"{run_name}_exp.csv")
    #
    #     # Cluster
    #     communities, graph, Q = phenograph.cluster(exp_df.iloc[:, 1:], k=200)
    #     pd.concat([exp_df.iloc[:, 0],
    #                pd.Series(communities, name="label")], axis=1)\
    #         .to_csv(out_dir / f"{run_name}_labels.csv", index=False)

    # Step 3: Visualize the data
    for str_dir in data_dirs:
        run_name = str_dir.split("_")[0]
        print(f"Visualizing {run_name}")

        # Get path
        path_dir = Path(str_dir)

        # Output dirs
        out_dir = path_dir / "viz"

        # Get data
        exp_df = pd.read_csv(out_dir / f"{run_name}_exp.csv")
        exp_df = exp_df[[exp_df.columns[0]] + list(exp_df.iloc[:, 1:].mean().sort_values().index)]
        loc_df = pd.read_csv(out_dir / f"{run_name}_loc.csv")
        rel_df = pd.read_csv(out_dir / f"{run_name}_rel.csv")
        label_df = pd.read_csv(out_dir / f"{run_name}_labels.csv")

        # Fig 1: Plot cell area
        cell_area = pd.DataFrame({"cell_area": loc_df.iloc[:, 1:]["Area"], "seg_method": ["Ilastik"] * len(loc_df)})
        plot_cell_area(cell_area, run_name, out_dir)

        # Fig 2: Plot cell mean intensity
        cell_mean = pd.melt(exp_df.iloc[:, 1:], var_name="Marker", value_name="Mean_Intensity")
        plot_cell_mean_intensity(cell_mean, run_name, out_dir)

        # Fig 3: Plot summary cluster heat map
        cluster_mean = pd.concat(
            [exp_df.iloc[np.where(label_df.label == i)[0], 1:].mean() for i in range(label_df.label.max()+1)], axis=1)
        plot_cluster_mean_heatmap(cluster_mean, run_name, out_dir)

        # Fig 4: Plot cell count
        cell_count = exp_df.ImageNumber.value_counts(sort=False).reset_index()
        cell_count.columns = ["ROI", "Count"]
        plot_cell_count(cell_count, run_name, out_dir)

        # Fig 5: Plot subset violins of Fig 2
        plot_cell_mean_intensity_zoomin(cell_mean, run_name, out_dir)

        # Fig 6: Co-expression
        for i in [1, 2, 3]:
            superset_exp = exp_df[exp_df.ImageNumber == i].iloc[:, 1:]
            for exp in ["coexp", "nonexp"]:
                if exp == "coexp":
                    subset_exp = superset_exp.iloc[:, [3, 12, 14, 15, 16, 17, 18, 19, 20, 21]]
                else:
                    subset_exp = superset_exp.iloc[:, [0, 1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 13]]

                subset_col = subset_exp.columns
                subset_exp_mean = subset_exp.mean()
                subset_loc = loc_df[loc_df.ImageNumber == i].iloc[:, 2:]

                fig, ax = plt.subplots(len(subset_col), len(subset_col), figsize=(20, 20), dpi=500)

                for a in combinations(range(0, len(subset_col)), 2):
                    ax[a[0], a[1]].scatter(subset_exp.iloc[:, a[0]], subset_exp.iloc[:, a[1]], s=0.1)

                    idx0 = np.where(subset_exp.iloc[:, a[0]] > subset_exp_mean.iloc[a[0]])[0]
                    idx1 = np.where(subset_exp.iloc[:, a[1]] > subset_exp_mean.iloc[a[1]])[0]
                    idx_coexp = list(set(idx0) & set(idx1))
                    idx_nonexp = list(set(range(len(subset_exp))) - (set(idx0) & set(idx1)))
                    subset_loc0 = subset_loc.iloc[idx_coexp, :]
                    subset_loc1 = subset_loc.iloc[idx_nonexp, :]

                    ax[a[1], a[0]].scatter(subset_loc0.X, subset_loc0.Y, s=0.1, c="blue")
                    ax[a[1], a[0]].scatter(subset_loc1.X, subset_loc1.Y, s=0.1, c="white")
                    if a[0] == 0:
                        ax[a[1], a[0]].set_ylabel(f"{subset_col[a[1]]}", fontsize=20)
                    if a[1] == len(subset_col) - 1:
                        ax[a[1], a[0]].set_xlabel(f"{subset_col[a[0]]}", fontsize=20)
                ax[0, 0].set_ylabel(f"{subset_col[0]}", fontsize=20)
                ax[len(subset_col) - 1, len(subset_col) - 1].set_xlabel(f"{subset_col[len(subset_col) - 1]}",
                                                                        fontsize=20)
                plt.savefig(out_dir / f"square_{exp}{i}.png")
