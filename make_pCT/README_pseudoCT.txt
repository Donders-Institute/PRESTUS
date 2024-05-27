INSTRUCTIONS

This code creates a pseudoCT starting from a T1w and a PETRA UTE image. The code is divided into three scripts that need to be run one after the other. Two scripts are in bash and one script is in Matlab.

To run the code you will need to install:
- SimNIBS
- FSL
- ANTs

How to make a pseudoCT:

1) run the segmentation through SimNIBS using as an input a T1 and a PETRA UTE scan (the UTE is used instead than the T2). In order to perform this step you can either use the SimNIBS GUI or you can use PRESTUS, because PRESTUS has the feature of performing the segmentation using SimNIBS.

2) run in bash make_pseudoCT_part1.sh

3) run in Matlab soft_tissue_value.m

4) run in bash make_pseudoCT_part2.sh

The directory of steps 2-4 is the output folder of SimNIBS, that is the m2m folder. The outputs of the SimNIBS segmentation will be stored in that folder, as well as the pseudoCT and the other outputs of this code.

# developed by Eleonora Carpino
# eleonora.carpino@icloud.com
# 14 June 2023
