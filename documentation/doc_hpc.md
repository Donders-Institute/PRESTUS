### Usage on high-performance computing clusters

The Donders provides access to a high-performance computing (HPC) cluster.High-performance computing (HPC) clusters are systems composed of interconnected computers (nodes) working together to solve complex computational problems. They enable parallel processing, allowing large-scale simulations, data analysis, and scientific research tasks to be completed more efficiently. Each node typically consists of multiple processors (CPUs/GPUs), memory, and storage.

PRESTUS is designed for HPC deployment: the most efficient way to run simulations is to start multiple simulations (e..g, for different individuals, target locations, or parameter setups) via parallel jobs. **PBS** (Portable Batch System) and **SLURM** ((Simple Linux Utility for Resource Management)) are tools that manage job scheduling, resource allocation, and execution in HPC systems. Both tools ensure efficient usage of cluster resources by managing multiple users and workloads. The Donders manages both PBS (historically) and SLURM (more recent) schedulers. To see current ressource usage, see https://grafana.dccn.nl/ (*intranet required*).

The following sections describe workflows at the Donders to start (interactive) jobs.

For more extensive documentation, see the [HCP wiki](https://hpc.dccn.nl/) (*intranet required*).

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

7. start MATLAB: `matlab`

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

7. Start MATLAB: `matlab`

8. PRESTUS scripts: use `**_slurm*`

9. in `terminal` check `squeue`; for migrating other commands see [this documentation](https://hpc.dccn.nl/docs/cluster_howto/compute_slurm.html#migrating-from-torque-pbs-to-slurm)