# Spatial-Proteomics-Framework
A framework for processing spatial proteomics data from images to analysis. Modified from the ImcSegmentationPipeline from Bodenmiller Group.

## Installation
1. Download [CellProfiler](https://cellprofiler.org/releases)
    - configure by selecting *Files -> Preferences*, set *CellProfiler plugins directory* to `ImcPluginsCP/plugins` and restart CellProfiler.

2. Download [Ilastik](https://www.ilastik.org/download.html)

3. Run `conda install imageio numpy pandas readimc scipy tifffile xtiff`

## Process 0: Pre-process experiment directory
1. Input Directory
    - Create a directory with some `DIRECTORY_NAME`, e.g. MB0002_1_345
    - Inside the directory, we require two files:
      1. **Image** file of either:
        - `mcd` file: avoid spaces in naming. This is the output from the Hyperion Imaging System
        - `tiff` file: avoid spaces in naming. This is the an alternative output.
      2. **Panel** file: must be named `panel.csv`. This should have four columns like as the following matrix.
          - Metal tag is the identifier on the antibody. 
          - Target indicates a protein measured. 
          - In the full column, the markers of interest are denoted with 1, others are 0.
          - In the ilastik column, the markers used for segmentation are denoted with 1, others are 0.

            | Metal Tag | Target | full | ilastik|
            |-----------|--------|------|--------|
            | Dy162 | CD3 | 1 | 0 |
            | Ir193 | DNA | 0 | 1 |
            | Pt196 | Membrane | 0 | 1 |
            | ... | ... | ... | ... |

2. Analysis Directory
   - In lieu of setting up an analysis directory manually, we automate the creation here.
   Depending on your image file type, use the following commands in your terminal:
   - `mcd` file: `python 0_pre_processing.py -d DIRECTORY_NAME`
   - `tiff` file: `python 0B_pre_processing.py -d DIRECTORY_NAME`

3. Output
    - All files will be in a directory named `DIRECTORY_NAME_analysis`
    - Inside the analysis directory, there should be five folders: cpinp, cpout, crops, ilastik, and ometiff.

## Process 1: Train pixel classifier
1. CellProfiler
    - Open `1_prepare_ilastik.cppipe`
    - Drop the directory `DIRECTORY_NAME_analysis/ilastik` into the window under *Images*
    - In *Output Settings*: Set input to `DIRECTORY_NAME_analysis/cpinp`
    - In *Output Settings*: Set output to `DIRECTORY_NAME_analysis/crops`
    - Press *Analyze Images*
2. Ilastik
    - Open the `Ilastik` software.
    - Create a *Pixel Classification* project. Put it in your analysis directory and name it as you wish.
    - Press *Input Data -> Add New... -> Add separate Image(s)* and select all .h5 files in `analysis/crops`.
    - Press *Feature Selection* and select all features with σ >= 1
    - Press *Training*:
        - Rename labels: "Label 1 - Nucleus",  "Label 2 - Cytoplasm", "Label 3 - Background".
        - Scroll down to *Raw Input* in the *Group Visibility* box, press up and down to check the different Ilastik channels indicated in *panel.csv* 
        - Label the pixels by pressing the label and then the brush icon, switch back and forth between channels under *Raw Input*, you can also erase labels using the eraser icon
        - Press *Live Update* to update segmentation
        - Turn on/off segmentation and prediction masks displayed on the image to check your work by pressing the eye icon in each category in the *Group Visibility* box
    - Press *Prediction Export* and configure as the following:
        - Source: Probabilities
        - Choose Export Image Settings: Convert to Data Type: unsigned 16-bit, check Renormalize, Format: tiff
        - Press export all.
    - Press *Batch Processing* and press *Raw Data Files...* and select all _s2.h5 files in `analysis/ilastik` then click on *Process all files*

## Process 2 Cell segmentation
1. CellProfiler
    - Open `2_segment_ilastik.cppipe`
    - Drop the directory `DIRECTORY_NAME_analysis/ilastik` into the window under Images
    - In *Output Settings*: Set input to `DIRECTORY_NAME_analysis/cpinp`
    - In *Output Settings*: Set output to `DIRECTORY_NAME_analysis/cpout`
    - Press *Analyze Images*

2. Output
    - All files will be in a directory named `DIRECTORY_NAME_analysis/cpout`
    - Masks and probabilities are in `cpout/masks/` and `cpout/probabilities/`

## Process 3: Cell quantification
1. CellProfiler
    - Open `3_measure_mask.cppipe`
    - IMPORTANT: If you are using a `tiff` file OR you are receiving errors in image matching rules:
         - In *Metadata*: press *No* on the question *Extract metadata?*
         - In *NamesAndTypes*: set *Image set matching method* to *Order* instead of Metadata
    - Drop the directory `DIRECTORY_NAME_analysis/cpout` into the window under Images
    - In *Output Settings*: Set input to `DIRECTORY_NAME_analysis/cpinp`
    - In *Output Settings*: Set output to `DIRECTORY_NAME_analysis/cpout`
    - Set channels in the **first** *MeasureObjectIntensityMultichannel* to the number of full channels in *panel.csv* (default is 40)
    - Set channels in the **first** *MeasureImageIntensityMultichannel* to the number of full channels in *panel.csv* (default is 40)
    - Press *Analyze Images*

2. Output
    - All files will be in a directory named `DIRECTORY_NAME_analysis/cpout`
    - Expression features are stored in `cpout/cell.csv`
    - Spatial features are stored in `cpout/Object relationships.csv`

## Process 4: Visualize cell area and intensity level for each marker
1. Run `python 4_cell_analysis.py -d DIRECTORY_NAME_analysis`

2. Output
    - all files will be in a directory named *DIRECTORY_NAME_analysis/viz*
