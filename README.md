# Mapping-spatial-organization-of-in-vitro-neuronal-networks-using-high-content-imaging
This repository contains all the scripts used in the project “Mapping Spatial Organization of In Vitro Neuronal Networks Using High-Content Imaging.” These scripts enable the analysis of complete in vitro neuronal networks.

The repository is organized into four main folders:

1. Analysis
   
Contains MATLAB scripts for the analysis of stitched wells using the Modular Image Analysis (MIA) framework.
MIA is divided into modules and uses parallel programming to analyze multiple wells simultaneously.
Nuclei Analysis module: differentiates between neuronal and non-neuronal cells.
Distance Analysis module: performs spatial analysis of cell distribution in the wells.

2. Cellpose
   
Contains scripts and models for segmentation of stitched wells.
NetworkDetection: segmentation of stitched wells using a pre-trained model (included). The model allows segmentation of DAPI and AnkG, enabling identification of neurons with axonal initial segments.
AIS_analysis: MATLAB script for analyzing axonal initial segment properties across the complete well.

3. Diameter_analysis
   
Contains Python scripts for calculating the diameter of each cell (based on the outputs from the Analysis folder).

4. Stitching
   
Contains ImageJ macros used for stitching tiled images into complete wells.
