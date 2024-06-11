# PseudoCT

PRESTUS supports the use of UTE-based images as a source of pseudo-Hounsfield units. This is currently only implemented for layered setups.

### INSTRUCTIONS for creating a pseudoCT from UTE scans

```create_pseudoCT.sh``` creates a pseudoCT and an associated mask file starting from a T1w and a PETRA UTE image. 

To run the code you will need to install:
- SimNIBS
- FSL
- ANTs

How to make a pseudoCT:

1) run the segmentation through SimNIBS using as an input a T1 and a PETRA UTE scan (the UTE is used instead than the T2). In order to perform this step you can either use the SimNIBS GUI or you can use PRESTUS, because PRESTUS has the feature of performing the segmentation using SimNIBS.
2) run in bash create_pseudoCT.sh

create_pseudoCT.sh calls the MATLAB function find_soft_tissue_peak

The pseudoCT will be deposited in the m2m folder alongside the SimNIBS segmentation. 

The following is an example script that you can use for a create_pseudoCT.sh call.

```
#!/bin/bash

rootpath="$(pwd)/.."
rootpath=$(builtin cd $rootpath; pwd)

# Load necessary modules
module load ants
module load matlab
module load fsl

# Change directory to the scripts directory
scriptpath=${rootpath}/tools/PRESTUS/functions
cd "${scriptpath}" || { echo "Directory not found"; exit 1; }

# Source the script containing the function
source ${scriptpath}/create_pseudoCT.sh

subject_id="001" # 'sub-' will automatically be added
m2m_path=${rootpath}/data/simnibs

# Call the create_pseudoCT function
create_pseudoCT "$subject_id" "$m2m_path" "${scriptpath}"
```