import argparse
from pathlib import Path
import pandas as pd
from readimc import MCDFile
from os import PathLike
import logging
from typing import Dict, Optional, Sequence, Union, List
from readimc.data import Acquisition, Panorama, Slide
import xtiff
import imageio
from dataclasses import dataclass
from xml.etree import ElementTree as ET
import numpy as np
import shutil
import tifffile
from scipy.ndimage import maximum_filter
import re


def extract_mcd_file(
    mcd_file: Union[str, PathLike],
    acquisition_dir: Union[str, PathLike],
    txt_files: Optional[Sequence[Union[str, PathLike]]] = None,
) -> pd.DataFrame:
    acquisition_origins = {}
    acquisition_is_valids = {}
    Path(acquisition_dir).mkdir(exist_ok=True)
    with MCDFile(mcd_file) as f_mcd:
        schema_xml_file = Path(acquisition_dir) / f"{Path(mcd_file).stem}_schema.xml"
        _extract_schema(f_mcd, schema_xml_file)
        for slide in f_mcd.slides:
            slide_stem = f"{Path(mcd_file).stem}_s{slide.id}"
            slide_img_file = Path(acquisition_dir) / f"{slide_stem}_slide.png"
            _extract_slide(f_mcd, slide, slide_img_file)
            for panorama in slide.panoramas:
                panorama_img_file = (
                    Path(acquisition_dir) / f"{slide_stem}_p{panorama.id}_pano.png"
                )
                _extract_panorama(f_mcd, panorama, panorama_img_file)
            for acquisition in slide.acquisitions:
                acquisition_img_file = (
                    Path(acquisition_dir)
                    / f"{slide_stem}_a{acquisition.id}_ac.ome.tiff"
                )
                acquisition_channels_file = acquisition_img_file.with_name(
                    acquisition_img_file.name[:-9] + ".csv"
                )
                acquisition_origin = "mcd"
                acquisition_is_valid = _extract_acquisition(
                    f_mcd, acquisition, acquisition_img_file, acquisition_channels_file
                )
                acquisition_origins[acquisition.id] = acquisition_origin
                acquisition_is_valids[acquisition.id] = acquisition_is_valid
        return _create_acquisition_metadata(
            f_mcd, acquisition_origins, acquisition_is_valids
        )

def _extract_schema(mcd_file_handle: MCDFile, schema_xml_file: Path) -> bool:
    try:
        with schema_xml_file.open("w") as f:
            f.write(mcd_file_handle.metadata)
        return True
    except Exception as e:
        logging.error(f"Error reading schema XML from file {mcd_file_handle.path.name}: {e}")
        return False


def _extract_slide(mcd_file_handle: MCDFile, slide: Slide, slide_img_file: Path) -> bool:
    try:
        slide_img = mcd_file_handle.read_slide(slide)
        if slide_img is not None:
            imageio.imwrite(slide_img_file, slide_img, compress_level=1)
        return True
    except Exception as e:
        logging.error(f"Error reading slide {slide.id} from file {mcd_file_handle.path.name}: {e}")


def _extract_panorama(mcd_file_handle: MCDFile, panorama: Panorama, panorama_img_file: Path) -> bool:
    try:
        panorama_img = mcd_file_handle.read_panorama(panorama)
        imageio.imwrite(panorama_img_file, panorama_img, compress_level=1)
        return True
    except Exception as e:
        logging.error(
            f"Error reading panorama {panorama.id} from file {mcd_file_handle.path.name}: {e}")
        return False


def _extract_acquisition(mcd_file_handle: MCDFile, acquisition: Acquisition,
                         acquisition_img_file: Path, acquisition_channels_file: Path) -> bool:
    try:
        acquisition_img = mcd_file_handle.read_acquisition(acquisition)
        _write_acquisition_image(
            mcd_file_handle, acquisition, acquisition_img, acquisition_img_file, acquisition_channels_file)
        return True
    except Exception as e:
        logging.error(f"Error reading acquisition {acquisition.id} from file {mcd_file_handle.path.name}: {e}")
        return False


def _write_acquisition_image(mcd_file_handle: MCDFile, acquisition: Acquisition,
                             acquisition_img: np.ndarray, acquisition_img_file: Path,
                             acquisition_channels_file: Path) -> None:
    xtiff.to_tiff(
        acquisition_img,
        acquisition_img_file,
        ome_xml_fun=get_acquisition_ome_xml,
        channel_names=acquisition.channel_labels,
        channel_fluors=acquisition.channel_names,
        xml_metadata=mcd_file_handle.metadata.replace("\r\n", ""),
    )
    pd.DataFrame(
        data={
            "channel_name": acquisition.channel_names,
            "channel_label": acquisition.channel_labels,
        }
    ).to_csv(acquisition_channels_file, index=False)


def _create_acquisition_metadata(mcd_file_handle: MCDFile, acquisition_origins: Dict[int, str],
                                 acquisition_is_valids: Dict[int, bool]) -> pd.DataFrame:
    return pd.DataFrame(
        data=[
            AcquisitionMetadata.from_mcd_file_acquisition(
                mcd_file_handle,
                acquisition,
                origin=acquisition_origins[acquisition.id],
                is_valid=acquisition_is_valids[acquisition.id],
            )
            for slide in mcd_file_handle.slides
            for acquisition in slide.acquisitions
        ]
    )


def create_analysis_stacks(
    acquisition_dir: Union[str, PathLike],
    analysis_dir: Union[str, PathLike],
    analysis_channels: Sequence[str],
    suffix: Optional[str] = None,
    hpf: Optional[float] = None,
) -> None:
    Path(analysis_dir).mkdir(exist_ok=True)
    for acquisition_img_file in Path(acquisition_dir).glob("*.ome.tiff"):
        acquisition_channels_file = acquisition_img_file.with_name(
            acquisition_img_file.name[:-9] + ".csv"
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
            acquisition_img_file.name[:-9] + ".tiff"
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


@dataclass
class AcquisitionMetadata:
    AcSession: str  # imctools
    slide_id: int  # imctools
    origin: str  # imctools
    source_path: str  # imctools
    id: int
    description: Optional[str]
    ablation_power: Optional[float]
    ablation_distance_between_shots_x: Optional[float]
    ablation_distance_between_shots_y: Optional[float]
    ablation_frequency: Optional[float]
    # AcquisitionROIID
    # OrderNumber
    signal_type: Optional[str]
    # DualCountStart
    # DataStartOffset
    # DataEndOffset
    start_timestamp: Optional[str]
    end_timestamp: Optional[str]
    # AfterAblationImageEndOffset
    # AfterAblationImageStartOffset
    # BeforeAblationImageEndOffset
    # BeforeAblationImageStartOffset
    roi_start_x_pos_um: Optional[float]
    roi_start_y_pos_um: Optional[float]
    roi_end_x_pos_um: Optional[float]
    roi_end_y_pos_um: Optional[float]
    movement_type: Optional[str]
    segment_data_format: Optional[str]
    # ValueBytes
    max_y: int
    max_x: int
    # PlumeStart
    # PlumeEnd
    template: Optional[str]
    profiling_type: Optional[str]
    has_before_ablation_image: bool  # imctools
    has_after_ablation_image: bool  # imctools
    is_valid: bool  # imctools

    def from_mcd_file_acquisition(
        f_mcd: MCDFile,
        acquisition: Acquisition,
        origin: str = "mcd",
        is_valid: bool = True,
    ) -> "AcquisitionMetadata":
        before_ablation_image_start_offset = acquisition.metadata.get(
            "BeforeAblationImageStartOffset"
        )
        before_ablation_image_end_offset = acquisition.metadata.get(
            "BeforeAblationImageEndOffset"
        )
        has_before_ablation_image = (
            before_ablation_image_start_offset is not None
            and before_ablation_image_end_offset is not None
            and int(before_ablation_image_start_offset)
            < int(before_ablation_image_end_offset)
        )
        after_ablation_image_start_offset = acquisition.metadata.get(
            "AfterAblationImageStartOffset"
        )
        after_ablation_image_end_offset = acquisition.metadata.get(
            "AfterAblationImageEndOffset"
        )
        has_after_ablation_image = (
            after_ablation_image_start_offset is not None
            and after_ablation_image_end_offset is not None
            and int(after_ablation_image_start_offset)
            < int(after_ablation_image_end_offset)
        )
        return AcquisitionMetadata(
            f_mcd.path.stem,
            acquisition.slide.id,
            origin,
            str(f_mcd.path.absolute()),
            acquisition.id,
            acquisition.description,
            float(acquisition.metadata.get("AblationPower") or "nan"),
            float(acquisition.metadata.get("AblationDistanceBetweenShotsX") or "nan"),
            float(acquisition.metadata.get("AblationDistanceBetweenShotsY") or "nan"),
            float(acquisition.metadata.get("AblationFrequency") or "nan"),
            acquisition.metadata.get("SignalType"),
            acquisition.metadata.get("StartTimeStamp"),
            acquisition.metadata.get("EndTimeStamp"),
            float(acquisition.metadata.get("ROIStartXPosUm") or "nan"),
            float(acquisition.metadata.get("ROIStartYPosUm") or "nan"),
            float(acquisition.metadata.get("ROIEndXPosUm") or "nan"),
            float(acquisition.metadata.get("ROIEndYPosUm") or "nan"),
            acquisition.metadata.get("MovementType"),
            acquisition.metadata.get("SegmentDataFormat"),
            int(acquisition.metadata["MaxY"]),
            int(acquisition.metadata["MaxX"]),
            acquisition.metadata.get("Template"),
            acquisition.metadata.get("ProfilingType"),
            has_before_ablation_image,
            has_after_ablation_image,
            is_valid,
        )


def get_acquisition_ome_xml(
    img: np.ndarray,
    image_name: Optional[str],
    channel_names: Optional[Sequence[str]],
    big_endian: bool,
    pixel_size: Optional[float],
    pixel_depth: Optional[float],
    channel_fluors: Optional[Sequence[str]] = None,
    xml_metadata: Optional[str] = None,
    **ome_xml_kwargs,
) -> ET.ElementTree:
    element_tree = xtiff.get_ome_xml(
        img,
        image_name,
        channel_names,
        big_endian,
        pixel_size,
        pixel_depth,
        **ome_xml_kwargs,
    )
    root_elem = element_tree.getroot()
    root_elem.set("Creator", "IMC Segmentation Pipeline")
    if channel_fluors is not None:
        assert len(channel_fluors) == img.shape[2]
        channel_elems = element_tree.findall("./Image/Pixels/Channel")
        assert channel_elems is not None and len(channel_elems) == img.shape[2]
        for channel_elem, channel_fluor in zip(channel_elems, channel_fluors):
            channel_elem.set("Fluor", channel_fluor)
    if xml_metadata is not None:
        structured_annot_elem = ET.SubElement(root_elem, "StructuredAnnotations")
        xml_annot_elem = ET.SubElement(structured_annot_elem, "XMLAnnotation")
        xml_annot_elem.set("ID", "Annotation:0")
        xml_annot_value_elem = ET.SubElement(xml_annot_elem, "Value")
        orig_metadata_elem = ET.SubElement(xml_annot_value_elem, "OriginalMetadata")
        orig_metadata_key_elem = ET.SubElement(orig_metadata_elem, "Key")
        orig_metadata_key_elem.text = "MCD-XML"
        orig_metadata_value_elem = ET.SubElement(orig_metadata_elem, "Value")
        orig_metadata_value_elem.text = xml_metadata
    return element_tree


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Preprocess MCD files from Fluidigm IMC.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-d', '--dirs', nargs='+', default=[], help="Directories to process", required=True)
    args = parser.parse_args()
    data_dirs = args.dirs

    # INPUTS example
    # data_dirs = ["data/20220330_chl_titration_slide_2", "data/20220728_chl_titration_slide_3"]

    """
    Pre-processing
    """
    for str_dir in data_dirs:
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

        # convert mcds to tiffs
        acquisition_metadatas = []
        mcd_files = list(path_dir.rglob("*.mcd"))
        for mcd in mcd_files:
            print(mcd)
            acquisition_metadata = extract_mcd_file(mcd, acquisitions_dir / mcd.stem)
            acquisition_metadatas.append(acquisition_metadata)
        acquisition_metadata = pd.concat(acquisition_metadatas, copy=False)
        acquisition_metadata.to_csv(cellprofiler_input_dir / "acquisition_metadata.csv")

        shutil.copy2(panel_file, cellprofiler_output_dir / "panel.csv")

        panel: pd.DataFrame = pd.read_csv(panel_file)

        for acquisition_dir in acquisitions_dir.glob("*"):
            if acquisition_dir.is_dir():
                # Write full stack
                create_analysis_stacks(acquisition_dir, final_images_dir,
                                       sort_channels_by_mass(
                                           panel.loc[panel[panel_keep_col] == 1, panel_channel_col].tolist()),
                                       suffix="_full", hpf=50.0)
                # Write ilastik stack
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
