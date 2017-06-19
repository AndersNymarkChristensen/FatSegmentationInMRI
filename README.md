# FatSegmentationInMRI
Matlab code for automatic segmentation of fat in MRI scans

A function for automatic formatting of DICOM to the form used på the functions are placed in 'auxiliary_functions' and named 'readAllFilesAsDICOM.m'. See the documentation in the file.

The bias-correction must be cloned/downloaded seperately from: https://github.com/cthla/iic to the 'biasCorr' folder
There is a seperate license and article that should be cited for this part

When the data has been formatted appropriately using the 'readAllFilesAsDICOM.m' the only function that need to be run is the one in the root: MRI_Abdominal_Segmentation.m

Further details are given in the header of the file.

Questions can be send to "anym AT dtu DOT dk"
