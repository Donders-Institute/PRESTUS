#!/bin/bash
# Compile kspaceFirstOrder-CUDA on Donders HPC  

# Donders HPC: load cuda module; version determines upper architecture limit
module load cuda/11.8

# Compile CUDA from source
cd $HOME/kwave-1.3-Donders/source/kspaceFirstOrder-CUDA
make clean
make -j$(nproc)

cp kspaceFirstOrder-CUDA $HOME/kwave-1.3-Donders/binaries/

# Ensure executable
chmod +x $HOME/kwave-1.3-Donders/binaries/kspaceFirstOrder-CUDA