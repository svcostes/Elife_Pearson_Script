THIS SCRIPT WAS USED FOR THE GENOME-WIDE RNAi SCREEN (Swenson, J. et al., 2016, In prep)
S. Costes <svcostes@lbl.gov>, January 2011, LBNL


This script looks into the target_dir and assumes all subdirectories have images and must be analyzed. It segments the nucleus using the “nuc_chan” and then measures various imaging features of the “nuc_chan” and two other channels (called “chan1” and “chan2”). Additionally, it computes the Pearson correlation coefficient in a pair-wise manner between all channels within each nucleus. Nuclear segmentation is automatic. nuc_rad is the expected radius of the nucleus, if set too low, it will results in fragmented nuclei in the segmentation.

To run first install Matlab (http://www.mathworks.com/products/matlab/) and the Dipimage MATLAB library (http://biocomp.lbl.gov/download.html for instructions).

The script expects a colorized (RGB) tif image. The default is for the channel you want to use for nuclear segmentation (“nuc_chan”) to be blue. Note that the channel used for nuclear segmentation can be modified on line 30 of “Pearson_screening_for_pub.m”. “chan1” and “chan2” are also modifiable (lines 31-32) and the default is that chan1 is the red channel and chan2 is the green channel.


