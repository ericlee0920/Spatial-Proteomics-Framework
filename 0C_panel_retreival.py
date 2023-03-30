import argparse
from pathlib import Path
import pandas as pd
from readimc import MCDFile
from os import PathLike
from typing import Union


def extract_panel_from_mcd(mcd_file: Union[str, PathLike], path_dir):
    with MCDFile(mcd_file) as f_mcd:
        for slide in f_mcd.slides:
            panel_df = pd.DataFrame([slide.acquisitions[0].channel_names, slide.acquisitions[0].channel_labels]).T
            panel_df.columns = ["Metal Tag", "Target"]
            panel_df = panel_df.assign(full=[0] * panel_df.shape[0])
            panel_df = panel_df.assign(ilastik=[0] * panel_df.shape[0])
            panel_df.to_csv(path_dir / "panel.csv", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Preprocess MCD files from Fluidigm IMC.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-d', '--dirs', nargs='+', default=[], help="Directories to process", required=True)
    args = parser.parse_args()
    data_dirs = args.dirs

    # INPUTS example
    # data_dirs = ["230213"]

    for str_dir in data_dirs:
        path_dir = Path(str_dir)
        mcd_files = list(path_dir.rglob("*.mcd"))
        for mcd in mcd_files:
            print(mcd)
            extract_panel_from_mcd(mcd, path_dir)
