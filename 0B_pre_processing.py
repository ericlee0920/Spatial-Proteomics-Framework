import argparse
from pathlib import Path
import pandas as pd
from os import PathLike
from typing import Optional, Sequence, Union, List
import numpy as np
import shutil
import tifffile
from scipy.ndimage import maximum_filter
import re


def create_analysis_stacks(
    acquisition_dir: Union[str, PathLike],
    analysis_dir: Union[str, PathLike],
    analysis_channels: Sequence[str],
    suffix: Optional[str] = None,
    hpf: Optional[float] = None,
) -> None:
    Path(analysis_dir).mkdir(exist_ok=True)
    for acquisition_img_file in Path(acquisition_dir).glob("*.tiff"):
        acquisition_channels_file = acquisition_img_file.with_name(
            acquisition_img_file.name[:-5] + ".csv"
        )
        acquisition_img = tifffile.imread(acquisition_img_file)
        assert acquisition_img.ndim == 3
        acquisition_channels: pd.DataFrame = pd.read_csv(acquisition_channels_file)
        assert len(acquisition_channels.index) == acquisition_img.shape[0]
        analysis_channel_indices = [
            acquisition_channels["channel_name"].tolist().index(channel_name)
            for channel_name in analysis_channels
        ]
        analysis_img = acquisition_img[analysis_channel_indices]
        analysis_img_file = Path(analysis_dir) / (
            acquisition_img_file.name[:-5] + ".tiff"
        )
        if suffix is not None:
            analysis_img_file = analysis_img_file.with_name(
                analysis_img_file.name[:-5] + suffix + ".tiff"
            )
        analysis_channels_file = analysis_img_file.with_suffix(".csv")
        if hpf is not None:
            analysis_img = filter_hot_pixels(analysis_img, hpf)
        tifffile.imwrite(
            analysis_img_file, data=analysis_img.astype(np.uint16), imagej=True
        )
        with analysis_channels_file.open("w") as f:
            f.write("\n".join(analysis_channels))


def filter_hot_pixels(img: np.ndarray, thres: float) -> np.ndarray:
    kernel = np.ones((1, 3, 3), dtype=bool)
    kernel[0, 1, 1] = False
    max_neighbor_img = maximum_filter(img, footprint=kernel, mode="mirror")
    return np.where(img - max_neighbor_img > thres, max_neighbor_img, img)


def sort_channels_by_mass(channels: Sequence[str]) -> List[str]:
    return sorted(channels, key=lambda channel: int(re.sub("[^0-9]", "", channel) or 0))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Preprocess MCD files from Fluidigm IMC.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-d', '--dirs', nargs='+', default=[], help="Directories to process", required=True)
    args = parser.parse_args()
    str_dir = args.dirs

    # INPUTS example
    # str_dir = "metabric1"

    """
    Pre-processing
    """
    # get path
    path_dir = Path(str_dir)

    # panel info
    panel_file = f"{str_dir}/panel.csv"
    panel_channel_col = "Metal Tag"
    panel_keep_col = "full"
    panel_ilastik_col = "ilastik"

    # working directory
    work_dir = Path(f"{str_dir}_analysis")
    work_dir.mkdir(exist_ok=True)

    # output dirs
    acquisitions_dir = work_dir / "ometiff"
    ilastik_dir = work_dir / "ilastik"
    crops_dir = work_dir / "crops"
    cellprofiler_input_dir = work_dir / "cpinp"
    cellprofiler_output_dir = work_dir / "cpout"
    final_images_dir = cellprofiler_output_dir / "images"
    final_masks_dir = cellprofiler_output_dir / "masks"
    final_probabilities_dir = cellprofiler_output_dir / "probabilities"

    # create dirs
    acquisitions_dir.mkdir(exist_ok=True)
    crops_dir.mkdir(exist_ok=True)
    ilastik_dir.mkdir(exist_ok=True)
    cellprofiler_input_dir.mkdir(exist_ok=True)
    cellprofiler_output_dir.mkdir(exist_ok=True)
    final_images_dir.mkdir(exist_ok=True)
    final_masks_dir.mkdir(exist_ok=True)
    final_probabilities_dir.mkdir(exist_ok=True)

    # processing tiffs
    tiff_files = list(path_dir.rglob("*.tiff"))
    file_name = str(tiff_files[0]).split('\\')[1]
    file_dir = acquisitions_dir / file_name[:-5]
    file_dir.mkdir(exist_ok=True)
    shutil.copy2(tiff_files[0], file_dir / file_name)
    panel_read = pd.read_csv(panel_file).iloc[:, :2]
    panel_read.columns = ["channel_name", "channel_label"]
    panel_read.to_csv(file_dir / f"{file_name[:-5]}.csv", index=False)

    shutil.copy2(panel_file, cellprofiler_output_dir / "panel.csv")

    panel: pd.DataFrame = pd.read_csv(panel_file)

    """
    Fix
    """

    for acquisition_dir in acquisitions_dir.glob("*"):
        if acquisition_dir.is_dir():
            # Write full stack
            create_analysis_stacks(acquisition_dir, final_images_dir,
                                   sort_channels_by_mass(
                                       panel.loc[panel[panel_keep_col] == 1, panel_channel_col].tolist()),
                                   suffix="_full", hpf=50.0)
            # # Write ilastik stack
            create_analysis_stacks(acquisition_dir, ilastik_dir,
                                   sort_channels_by_mass(
                                       panel.loc[panel[panel_ilastik_col] == 1, panel_channel_col].tolist()),
                                   suffix="_ilastik", hpf=50.0)

    # copy a file that contains the correct order of channels for the full stacks to the `cellprofiler_input_dir`.
    first_channel_order_file = next(final_images_dir.glob("*_full.csv"))
    shutil.copy2(first_channel_order_file, cellprofiler_input_dir / "full_channelmeta.csv")

    # generate channel metadata for the probability stack.
    probab_meta = ["CellCenter", "CellBorder", "Background"]
    with open(cellprofiler_input_dir / "probab_channelmeta_manual.csv", "w") as f:
        f.write("\n".join(probab_meta))
