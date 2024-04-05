# Visitor Gene Pseudo Time

MatLab scripts for the reconstruction of pseudo-time lapses of gene-cluster visits from fluorescence microscopy data obtained from fixed (non-living) samples

Currently, the repository contains the scripts for extracting and analyzing Nikon .nd2 files with multiple XY positions in them.

# Requirements

You will need MatLab installed on your computer. You will also need the **Image Processing toolbox** and the **Statistic and Machine Learning toolbox**. The **Parallel toolbox** is helpful, as it can speed up your image processing.

You also first have to download bfmatlab, and add it to your MatLab path.

bfmatlab:
https://www.openmicroscopy.org/bio-formats/downloads/

Adding a directory to the MatLab path permanently:
1. Click the option "Set Path" in the Tab "Home"
2. Navigate to the place where you save bfmatlab, open the directory and add to path
3. Click the button "Preferences" in the Tab "Home"
4. Select the bfmatlab directory, then click the button "Save"

# Use instructions

1. Place all the imaging data you would like to process into directories that have a name that represents the according experimental conditions. These directory names lateron will serve as identifier names for the different conditions, so please choose well.
2. In the script MultiPosition_extraction_nd2.m specifiy which experimental conditions should be extracted from which folders. This search is recursive, and finds all .nd2 files that are stoored within a given directorz as well as its subdirectories. In this file, you also should specify which color channel represents which label. The correct order is Pol II Ser2Phos, Pol II Ser5Phos, and oligopaint signal.
3. you can use the script ReviewExtractedStacks_Oligopaint.m to review the etracted data. This is optional. Dat that are of insufficient technical quality can be removed from your hard drive in a targeted fashion based on this script, it indicates the file names for everything it shows.
4. The actual image processing then is executed by the script Segmentation_FeatureExtraction.m. You can also adjust some of the analysis parameters to refine the analysis. You can also comment out parfor and uncomment the regular for expression, and uncomment waitforbuttonpress to review the analysis performance.
5. Actual analysis results and figures can then be prudced with the downstream analysis scripts listed below.

# Script sets

Data extraction:

MultiPosition_extraction_nd2.m

ReviewExtractedStacks_Oligopaint.m

Image analysis: 

Segmentation_FeatureExtraction.m

Downstream analysis and figures:

PseudoTime_Trajectory.m

GeneMap_Interpolation.m

VisitorGene_PseudoTimeSorting.m

# Example data set

Example image data for oligopaint labeling of the gene klf2b is provided by the Hilbert lab via Zenodo:
[https://zenodo.org/record/5268833](https://zenodo.org/records/5268833)

First, create a folder called ImageData, which you place in the directory where you have saved the MatLab scripts. It is the easiest way to keep the scripts and data together in one and the same place.

Second, inside the ImageData folder, create folders for several developmental stages, at least Oblong and Sphere, if you like also Dome, Epi_30, and Epi_50.

[To be extended...]

# Acknowledgments

The scripts make use of the rdir.m scripts for recursive directory search

https://www.mathworks.com/matlabcentral/fileexchange/47125-rdir-m
