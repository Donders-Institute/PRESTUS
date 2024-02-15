This is a set of functions for ultrasonic simulations including brain preprocessing, extracting data from Localite neuronavigation files, and more. 

The examples folder contains some examples of how to run the pipeline, including a tutorial.

# Changelog

19/01/23 - a big upgrade to support the latest versions of SimNIBS, MATLAB, and k-Wave. 

# Requirements

- MATLAB R2022b
- SimNIBS 4
- k-Wave 1.4 (should be installed in the 'toolbox' folder of PRESTUS or be added automatically on MATLAB startup if you use HPC)
- small toolboxes included in the repository

Tested on MATLAB 2022b, set up to work on Donders HPC (can work on a local PC as well). 

Before using the package, you need to have some libraries on your path. Three small ones are included in the toolbox folder, but you also would need to have [k-Wave](http://www.k-wave.org/download.php) and SimNIBS (https://simnibs.github.io/simnibs/build/html/index.html#simnibs-4) installed. **Please put k-Wave in the toolbox folder of PRESTUS or ensure that the paths to it are added on MATLAB startup if you are running PRESTUS on HPC**.

# How to run

The main function is single_subject_pipeline, which takes the subject ID and parameter structure as an input. The parameters are taken from the config files in the configs directory, see the configs in this folder for the information on configurable fields. Parameters include `results_filename_affix` field that can define a condition (e.g., stimulation site, protocol, or something else) in case several simulations are to be started for a given subject. 

For each subject, two files are needed to start, T1 and T2 scans. These files should allow identifying a subject based on the T1 & T2 filename templates set in the config files. These can include wildcards (\* - stands for any symbols repeated 0 or more times, like in bash or MATLAB dir command to be precise) and string substitution patterns (only `subject_id` is used for them for now). Put these files in the folder defined in data_path field in the configuration. If you want to use Localite data, you can either use an instrument marker file (put it in the same folder, add the subject id to it so that the `localite_instr_file_template` from the config file could be parsed) or trigger markers file, which you would need to preprocess for each subject (see example in `example_pipeline_sjoerd.m`). 

On the first run, it is recommended to start the scripts on a local computing node (i.e., without qsub) and without a GPU. Then (if things are good), the segmentation will be started for each subject. They take ~4 hours. After they are completed, start the jobs again, now remotely and with the GPU. This should ideally result in multiple files for each subject/condition in your `sim_outputs` folder. Intermediate plots and files are also created there. You can then run extra post-processing steps if needed.

# Processing steps

The main pipeline depends a bit on whether you set `water` or `layered` as the simulations medium. 

The first big part is to get the transducer and the focus positions and a segmentated head image if needed.

For `layered`:
1) the brain is segmented with SimNIBS
2) the segmented image is rotated so that the focal axis is aligned with the z-axis 
3) the segmented image is upsampled to match the simulations grid voxel size
4) the layers in the segmented image are smoothed, the gaps between skin and skull filled, the holes in the skull closed
5) the segmented image is cropped (based on the expanded cerebro-spinal fluid mask to get rid of part of the neck bones)

Steps 2&3 are done simultaneously to avoid the need to interpolate the image twice. 
The simulation grid is then created based on the segmented image size (with padding), and the transducer and the focus positions are computed based on the transformed coordinate system. 

For water-only simulations:
The simulation grid dimensions are taken from the parameter structure. The transducer and the focus positions are either taken from the same structure or, if they are missing there, computed so that transducer is at the z-axis and the focus is on the same axis with the distance based on the  expected_focal_distance_mm field in the parameters.

After that things are straightforward: k-Wave medium, source, grid, and sensor are set up and the simulations are started. When they are completed, the maximum ISPPA map is computed and the maximum pressure and intensity points are estimated for different masks. The heating simulations can be enabled by setting `run_heating_simulations` to 1. If `run_posthoc_water_simulations` is set in parameters, then after the layered simulations, the same simulation parameters are used to run the simulations for water only.

# Troubleshooting

In case of problems, a) look at the plots in the data folder, do they look fine? b) if the simulations use the cluster, look for error logs in `batch_job_logs` folder under your `data_path`; c) re-run the job locally to find errors. Sometimes there are errors related (I think) to differences in GPU versions on the cluster, so if the job fails during the simulations, you can either try to simply restart it or run it in an interactive GPU-enabled session (qsub -I with GPU, then load MATLAB/R2022b). If the grid is too big, the GPU might run into memory issues. In this case, adjust csf_mask_expansion_factor to reduce the grid size (but make sure that the whole skull fits in).

Finally, if the simulations on the GPU run more than two hours, it is likely that something is wrong. Most likely, you forgot to switch interactive flag to zero and they are stuck because a confirmation is required from the user to start. 
