# Spatial-Proteomics-Framework
A framework for processing spatial proteomics data from images to analysis. Modified from the ImcSegmentationPipeline from Bodenmiller Group.

## Installation
1. Download [CellProfiler](https://cellprofiler.org/releases)
    - configure by selecting *Files -> Preferences*, set *CellProfiler plugins directory* to `ImcPluginsCP/plugins` and restart CellProfiler.

2. Download [Ilastik](https://www.ilastik.org/download.html)

3. Run `conda install imageio numpy pandas readimc scipy tifffile xtiff`

## Process 0: Pre-process experiment directory
1. Input
    - create a directory with some `DIRECTORY_NAME`
    - mcd file: avoid spaces in naming. This is the output from the Hyperion Imaging System
    - panel file: must be named `panel.csv`. This should have four columns like:
      
      | Metal Tag | Target | full | ilastik|
      |-----------|--------|------|--------|
      | Dy162 | CD3 | 1 | 0 |
      | Ir193 | DNA | 0 | 1 |
      | Pt196 | Membrane | 0 | 1 |
      | ... | ... | ... | ... |

2. Run `python 0_pre_processing.py -d DIRECTORY_NAME`

3. Output
    - all files will be in a directory named `DIRECTORY_NAME_analysis` 

## Process 1: Train pixel classifier
1. CellProfiler
    - drop the directory `DIRECTORY_NAME_analysis/ilastik` into the window under *Images*
    - set *Output Settings* to `DIRECTORY_NAME_analysis/crops`
    - press *Analyze Images*
2. Ilastik
    - create a *Pixel Classification* project.
    - press *Input Data -> Add New... -> Add separate Image(s)* and select all .h5 files in `analysis/crops`.
    - press *Feature Selection* and select all features with σ >= 1
    - press *Training*:
        - rename labels: "Label 1 - Nucleus",  "Label 2 - Cytoplasm", add extra labels if appropriate e.g. background.
        - scroll down to *Raw Input* in the *Group Visibility* box, press up and down to check the different Ilastik channels indicated in *panel.csv* 
        - label the pixels by pressing the label and then the brush icon, switch back and forth between channels under *Raw Input*, you can also erase labels using the eraser icon
        - press *Live Update* to update segmentation
        - turn on/off segmentation and prediction masks displayed on the image to check your work by pressing the eye icon in each category in the *Group Visibility* box
    - press *Prediction Export* and configure as the following:
        - Source: Probabilities
        - Choose Export Image Settings: Convert to Data Type: unsigned 16-bit, check Renormalize, Format: tiff
    - press *Batch Processing* and press *Raw Data Files...* and select all _s2.h5 files in `analysis/ilastik` then click on *Process all files*

## Process 2 Cell segmentation
1. CellProfiler
    - drop the directory `DIRECTORY_NAME_analysis/ilastik` into the window under Images
    - set *Output Settings* to `DIRECTORY_NAME_analysis/cpout`
    - press *Analyze Images*

2. Output
    - all files will be in a directory named *DIRECTORY_NAME_analysis/cpout*
    - masks and probabilities are in `masks/` and `probabilities/`

## Process 3: Cell quantification
1. CellProfiler
    - drop the directory `DIRECTORY_NAME_analysis/cpout` into the window under Images
    - set *Input Settings* to `DIRECTORY_NAME_analysis/cpinp`
    - set *Output Settings* to `DIRECTORY_NAME_analysis/cpout`
    - set channels in the **first** *MeasureObjectIntensityMultichannel* to the number of full channels in *panel.csv* (default is 40)
    - set channels in the **first** *MeasureImageIntensityMultichannel* to the number of full channels in *panel.csv* (default is 40)
    - optionally set channels in the **second** *MeasureObjectIntensityMultichannel* to the number of labels made in Ilastik (e.g. nucleus, cytoplasm, background default is 3)
    - optionally set channels in the **second** *MeasureImageIntensityMultichannel* to the number of labels made in Ilastik (e.g. nucleus, cytoplasm, background default is 3)
    - press *Analyze Images*

2. Output
    - all files will be in a directory named *DIRECTORY_NAME_analysis/cpout*
    - spatial features are in `cell.csv`

## Process 4: Visualize cell area and intensity level for each marker
1. Run `python 4_cell_analysis.py -d DIRECTORY_NAME_analysis`

2. Output
    - all files will be in a directory named *DIRECTORY_NAME_analysis/viz*
