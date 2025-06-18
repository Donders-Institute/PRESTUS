### Usage on high-performance computing clusters

The Donders provides access to a high-performance computing (HPC) cluster.High-performance computing (HPC) clusters are systems composed of interconnected computers (nodes) working together to solve complex computational problems. They enable parallel processing, allowing large-scale simulations, data analysis, and scientific research tasks to be completed more efficiently. Each node typically consists of multiple processors (CPUs/GPUs), memory, and storage.

PRESTUS is designed for HPC deployment: the most efficient way to run simulations is to start multiple simulations (e..g, for different individuals, target locations, or parameter setups) via parallel jobs. **PBS** (Portable Batch System) and **SLURM** ((Simple Linux Utility for Resource Management)) are tools that manage job scheduling, resource allocation, and execution in HPC systems. Both tools ensure efficient usage of cluster resources by managing multiple users and workloads. The Donders manages both PBS (historically) and SLURM (more recent) schedulers. To see current ressource usage, see https://grafana.dccn.nl/ (*intranet required*).

The following sections describe workflows at the Donders to start (interactive) jobs.

For more extensive documentation, see the [HPC wiki](https://hpc.dccn.nl/) (*intranet required*).

### PBS

1. Select an access node:

- `mentat001`
- `mentat002`
- `mentat003`
- `mentat004`

2. Login on access node:
`ssh abcxyz@mentat004.dccn.nl`

3. Start the VNC manager to connect GUI via VNC:
`vncmanager - 2`

4. Connect via VNC (see [here](https://intranet.donders.ru.nl/index.php?id=vnc00&no_cache=1&sword_list%5B%5D=VNC))
    TigerVNC: enter ip [e.g., mentat004.dccn.nl:12]

5. Load QSUB:
`module load qsub`

6. Start interactive job:
`qsub -I -l 'nodes=1:gpus=1,feature=cuda,walltime=05:00:00,mem=24gb,reqattr=cudacap>=5.0'`

7. Start MATLAB: 
`module load matlab/R2022b`
`matlab`

    *Note*: PBS only supports CUDA 11.2, which is dropped starting in R2023, see [this issue](https://github.com/Donders-Institute/PRESTUS/issues/50)

8. PRESTUS scripts: use `**_qsub*`
9. in `terminal` check `qstat`

### SLURM

1. Select an access node:

- `mentat005`
- `mentat006`
- `mentat007` (previously `mentat001s`)

2. Login on access node:
`ssh abcxyz@mentat001.dccn.nl`

3. Start the VNC manager to connect GUI via VNC:
vncmanager - 2

4. Connect via VNC (see [here](https://intranet.donders.ru.nl/index.php?id=vnc00&no_cache=1&sword_list%5B%5D=VNC))
    TigerVNC: enter ip [e.g., mentat007.dccn.nl:12]

5. Load slurm:
`module load slurm`

6. Start interactive job: 
    - 6a. If you do *not* need a MATLAB GUI:
    `srun --mem=8gb --time=01:00:00 --x11 -p interactive --pty bash -i`
    - 6b. **If you need a MATLAB GUI:**
    `srun --partition=gpu --gres=gpu:1 --mem=8G --time=01:00:00 --x11 --pty /bin/bash -i`

7. Start MATLAB: 
`module load matlab/R2024a`
`matlab`

    *Note*: SLURM supports CUDA 12.2, which is why recent MATLAB versions (at least up to R2024) should be supported, see [this issue](https://github.com/Donders-Institute/PRESTUS/issues/50)

7. Start MATLAB: `matlab`

8. PRESTUS scripts: use `**_slurm*`

9. in `terminal` check `squeue`; for migrating other commands see [this documentation](https://hpc.dccn.nl/docs/cluster_howto/compute_slurm.html#migrating-from-torque-pbs-to-slurm)

### Note on GPUs

TUS simulations are accelerated by GPUs, but requesting GPUs can lead to longer wait times as the current concurrent GPU limit per user is 4. To reduce wait times, it is possible to run acoustic simulations first (to confirm targeting) because these require less RAM) and then run thermal simulations with the final protocol. Avoid blanket simulations (e.g., circling through all participants with all permutations) especially for thermal simulations.

Given that PRESTUS is a MATLAB toolbox, it currently only supports *Nvidia GPUs*. When Nvidia GPUs are digitally partitioned, there appears to be an issue with identifying the assigned GPU in MATLAB R2024+. For SLURM jobs, MATLAB R2023b is currently deployed by default.

The following settings can be used to specify the HPC GPU setup.

| Field                           | Default | Explanation                  |
|-------------------------------|-------------|-----------------------------|
| parameters.hpc_gpu                   | "gpu:1"       | Specific GPUs could be requested here (e.g.,```"nvidia_a100-sxm4-40gb:1"```, but this is not recommended. ```scontrol show nodes \| egrep -o gres/gpu:.*=[0-9] \| egrep -o 'nvidia_.*=' \| sort \| uniq \| sed 's/=//'``` lists available GPU types.|
| parameters.hpc_partition                   | "gpu"       | The Donders HPC offers a ```gpu40g``` partition that should be used for the majority of thermal simulations.  It consists of nodes with GPU with vRAM > 40 GB.|
| parameters.hpc_reservation                   | ""       | By default do not use a reserved cue.|

#### Benchmark data Nvidia GPUs

Below are (potentially unrepresentative) benchmark data for different Nvidia GPUs. These simulations were run on a 256 by 216 by 192mm grid, with minor varitions depending on transducer placement. 

**Acoustic Simulations**

| GPU                           | Memory Used | Duration                  | Notes                          |
|-------------------------------|-------------|-----------------------------|--------------------------------|
| A100 80 GB                   | 12 GB       | 14 mins                     | Slightly smaller grid size     |
| A100 80 GB (partitioned in 2x40GB) | 12 GB       | 25 mins                     |                                |
| A100 40 GB                   | 12 GB       | 19 mins                     |                                |
| A16 16 GB                    | 12 GB       | 145 mins                    |                                |
| P100 16 GB                   | 12 GB       | 43 mins                     | Compiled, ~ L40s              |
| L40S 47 GB                   | 14 GB       | 28 mins                     | Compiled                      |

---

**Heating Simulations**

| GPU                           | Memory Used    | Duration                 | Notes                          |
|-------------------------------|----------------|-----------------------------|--------------------------------|
| A100 80 GB                   | ?? GB          | ~12s/trial (400 trials: ~90 mins) | Slightly smaller grid size     |
| A100 40 GB                   | ?? GB          | ~17s/trial (400 trials: ~115 mins) |                                |
| A16 16 GB                    | Out of RAM     | ???                         |                                |
| P100 16 GB                   | Out of RAM     | ???                         |                                |
| L40S 47 GB                   | ?? GB          | ~4s/trial (400 trials: ~45 mins)   |                                |
