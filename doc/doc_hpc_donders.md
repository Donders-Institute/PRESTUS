# Donders HPC

The Donders cluster runs SLURM, and PRESTUS submits its jobs there itself ([Running on a Cluster](doc_hpc.md)). This page covers what is specific to the Donders cluster: starting an interactive session to submit from, writing your own MATLAB job scripts, and the PBS (qsub) queue.

On this page:

- [SLURM Workflow](#slurm-workflow): an interactive session from which you submit PRESTUS jobs.
- [Your Own MATLAB Job Scripts](#your-own-matlab-job-scripts): job scripts you submit yourself, and the launcher's argument pitfall.
- [PBS Workflow](#pbs-workflow): the same on the PBS (qsub) queue.

!!! note "At a glance"
    - **When you need it:** When you run PRESTUS on the Donders cluster.
    - **Before you start:** A Donders HPC account and access to the intranet.

Current resource usage is shown at https://grafana.dccn.nl/ (*intranet required*); the [HPC wiki](https://hpc.dccn.nl/) (*intranet required*) has the full documentation. The workflows below start an interactive session from which you submit PRESTUS jobs.

## SLURM Workflow

1. Log in to an access node (`mentat005`, `mentat006`, `mentat007`): `ssh abcxyz@mentat005.dccn.nl`
2. Start the VNC manager (`vncmanager - 2`) and connect via [VNC](https://intranet.donders.ru.nl/index.php?id=vnc00&no_cache=1&sword_list%5B%5D=VNC). [TigerVNC](https://tigervnc.org/): enter ip [e.g. `mentat007.dccn.nl:12`]
3. Load SLURM: `module load slurm`
4. Start an interactive job:
    - Without MATLAB GUI: `srun --mem=8gb --time=01:00:00 --x11 -p interactive --pty bash -i`
    - With MATLAB GUI: `srun --partition=gpu --gres=gpu:1 --mem=8G --time=01:00:00 --x11 --pty /bin/bash -i`
5. Start MATLAB: `module load matlab/R2023b` then `matlab`. SLURM supports [CUDA](https://developer.nvidia.com/cuda/toolkit) 12.2, see [PRESTUS issue #50](https://github.com/Donders-Institute/PRESTUS/issues/50). The submitted jobs load `matlab/R2023b` whatever version this session runs, because R2024+ can fail to identify partitioned GPUs ([GPU support](doc_hpc.md#gpu-support)).
6. Run PRESTUS with `platform: slurm` (or `auto`). Check job status with `squeue`. For PBS→SLURM command migration see [the Donders HPC wiki](https://hpc.dccn.nl/docs/cluster_howto/compute_slurm.html#migrating-from-torque-pbs-to-slurm).

## Your Own MATLAB Job Scripts

PRESTUS writes its own job scripts. For a script you submit yourself (e.g. a library rebuild), put all MATLAB code in a `.m` file and pass only its name:

```bash
#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --output=/project/<project>/logs/my_job_%j.log
module load matlab/R2023b
cd /project/<project>/code   # folder that holds my_job.m
matlab -batch my_job
```

The Donders `matlab` launcher splits its arguments at spaces, so `matlab -batch "a = 1; my_job(a)"` fails before MATLAB starts (`syntax error near unexpected token`, job exit code 2). Pass parameters as environment variables (`getenv`) or write one `.m` file per job. Use `/project/...` paths, and load `matlab/R2023b`: PRESTUS needs MATLAB's Java (YAML loader), and some `-batch` starts under other versions ran without it.

## PBS Workflow

1. Log in to an access node (`mentat001`, `mentat002`, `mentat003`, `mentat004`): `ssh abcxyz@mentat004.dccn.nl`
2. Start the VNC manager (`vncmanager - 2`) and connect via [VNC](https://intranet.donders.ru.nl/index.php?id=vnc00&no_cache=1&sword_list%5B%5D=VNC). TigerVNC: enter ip [e.g. `mentat004.dccn.nl:12`]
3. Load QSUB: `module load qsub`
4. Start an interactive job:

    ```
    qsub -I -l 'nodes=1:gpus=1,feature=cuda,walltime=05:00:00,mem=24gb,reqattr=cudacap>=5.0'
    ```

5. Start MATLAB: `module load matlab/R2022b` then `matlab`. PBS only supports CUDA 11.2, which is dropped from R2023 on; see [PRESTUS issue #50](https://github.com/Donders-Institute/PRESTUS/issues/50).
6. Run PRESTUS with `platform: qsub` (or `auto`). Check job status with `qstat`.
