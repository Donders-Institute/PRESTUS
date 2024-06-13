# How to run a simulation

The main function is single_subject_pipeline, which takes the subject ID and parameter structure as an input. The parameters are taken from the config files in the configs directory, see the configs in this folder for the information on configurable fields. Parameters include `results_filename_affix` field that can define a condition (e.g., stimulation site, protocol, or something else) in case several simulations are to be started for a given subject. 

For each subject, two files are needed to start, T1 and T2 scans. These files should allow identifying a subject based on the T1 & T2 filename templates set in the config files. These can include wildcards (\* - stands for any symbols repeated 0 or more times, like in bash or MATLAB dir command to be precise) and string substitution patterns (only `subject_id` is used for them for now). Put these files in the folder defined in data_path field in the configuration. If you want to use Localite data, you can either use an instrument marker file (put it in the same folder, add the subject id to it so that the `localite_instr_file_template` from the config file could be parsed) or trigger markers file, which you would need to preprocess for each subject (see example in `example_pipeline_sjoerd.m`). 

On the first run, it is recommended to start the scripts on a local computing node (i.e., without qsub) and without a GPU. Then (if things are good), the segmentation will be started for each subject. They take ~4 hours. After they are completed, start the jobs again, now remotely and with the GPU. This should ideally result in multiple files for each subject/condition in your `sim_outputs` folder. Intermediate plots and files are also created there. You can then run extra post-processing steps if needed.

# Processing steps

The main pipeline depends on whether you set `water` or `layered` as the simulations medium. 

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

## Supported simulation setups

- `layered`: the simulation medium varies in terms of water, brain, skin and skull, with potential subdivision into cortical bone (`skull_cortical`) and trabecular bone (`skull_trabecular`)
- `layered`+`pseudoCT`: acoustic properties of bone are based on pseudo-Hounsfield units
- `water`: homogeneous water medium
- `brain`: homogeneous brain medium
- `water_and_skull`: homogeneous water medium + skull
- `brain_and_skull`: homogeneous water medium + brain; note that brain medium is also assumed outside the skull