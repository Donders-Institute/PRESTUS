#!/bin/bash
# Set up CUDA compilation on Donders HPC

cd $HOME

# First create the folder kwave-1.3-Donders
mkdir -p kwave-1.3-Donders/{binaries,source}

# Download CUDA source files (id=207) - name the OG download kwave-1.3-source
wget "http://www.k-wave.org/getfile.php?id=207" -O kwave-1.3-Donders/source/kwave-1.3-source.tar.gz

# Unpack to source dir
unzip kwave-1.3-Donders/source/kwave-1.3-source.tar.gz -d ~/kwave-1.3-Donders/source/
rm kwave-1.3-Donders/source/kwave-1.3-source.tar.gz

####################################################
############### Make MANUAL adjustments to MAKE file
####################################################

# # Everything will be linked dynamically
# LINKING = DYNAMIC

# # Set up paths
# CUDA_DIR = $(CUDA_PATH)
# HDF5_DIR = /usr
# ZLIB_DIR = /usr  
# SZIP_DIR = /usr

# # Select CPU architecture
# CPU_ARCH = AVX2

# # What CUDA GPU architectures to include in the binary
# Note that the upper limit will be determined by the CUDA version
# (set with module load cuda/xxx)
# CUDA_ARCH = --generate-code arch=compute_60,code=sm_60 \
#             --generate-code arch=compute_70,code=sm_70 \
#             --generate-code arch=compute_72,code=sm_72 \
#             --generate-code arch=compute_75,code=sm_75 \
#             --generate-code arch=compute_80,code=sm_80 \
#             --generate-code arch=compute_86,code=sm_86 \
#             --generate-code arch=compute_89,code=sm_89 \
#             --generate-code arch=compute_90,code=sm_90

