#!/bin/bash

rootpath="$(pwd)/.." # move to functions directory
rootpath=$(builtin cd $rootpath; pwd)

# Load necessary modules
module load ants
module load matlab
module load fsl

# Change directory to the scripts directory
scriptpath=${rootpath}
cd "${scriptpath}" || { echo "Directory not found"; exit 1; }

# Source the script containing the function
source ${scriptpath}/pct/pct_create_pseudoCT.sh

for i in {1..1}; do # run only for sub-001
    subject_id=$(printf "%03d" ${i})
    m2m_path=${rootpath}/data/simnibs
    echo "$subject_id"

    # Call the create_pseudoCT function
    pct_create_pseudoCT "$subject_id" "$m2m_path" "/opt/matlab/R2022b/bin/matlab" "miscouridou" "1" # use miscouridou mapping
done
