# Tasmanian Deep Seated Landslide Probabilistic Runout Model tool

The primary aim of this project was to produce a geospatial probabilistic model for deep-seated landslide runout for Tasmania by deriving parameter distributions like reach angle, runout distance, and regression distance from the existing MRT Landslide Database, relying on mapped source areas. This is to address the lack of current probabilistic runout, required to consider variables currently not assessed, such as weathering and lithology.
This repository houses a QGIS plugin to make this approach accessible to geoscientists and planners for landslide hazard management across different landscapes. 
# Getting Started 
Prerequsites
QGIS version 4.0 and above and Python environment compatiable with QGIS
# How to use this tool
Download the "ProbabilisticLandslideRunout" folder as zip and install as zip into QGIS through QGIS's plugin repository. 
# Important notes
All rasters being placed into the model must have the same extent, coordinate system and resolution. This can be arranged through QGIS's built in 'Arrange Raster tool'. A machine with large memory or lots of patience is required for probabilistic runout modelling.
# Usage
A folder with a 2m DEM, Landslide Source Raster and Alpha angle (geology) layer is provided for testing from a case study from the West Tamar Area in Northern Tasmania. These can be downloaded and used

