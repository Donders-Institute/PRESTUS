clear
cd /project/3015999.02/andche_sandbox/orca-lab/project/tuSIM/ % change path here

% add paths
addpath('functions')
addpath(genpath('toolboxes')) 
addpath('/home/common/matlab/fieldtrip/qsub') % uncomment if you are using Donders HPC


%% create images in MNI space
parameters = load_parameters('sjoerd_config_opt_CTX500-026_73.5mm.yaml');
parameters.simulation_medium = 'layered';

thresholds = [75, 85];

for watts = [40, 50]
    for hertz = [250, 500]
        outputs_path = sprintf('/project/3023001.06/Simulations/all_files_with_output/%iKHz/%iW_per_cm2/', hertz, watts);
        files = dir(outputs_path);

        subject_list = [];
        for i = 1:length(files)
            fname = files(i).name;
            if regexp(fname, '^sub-\d+$') 
                subject_list = [subject_list str2num(fname(5:7))];
            end
        end
        for threshold = thresholds
            for target = ["left_amygdala","right_amygdala"]
                fprintf('Target %s\n', target)
                if strcmp(target,"left_amygdala") % labels are swapped
                    target_lbl = 'L';
                else
                    target_lbl = 'R';
                end
                parameters.results_filename_affix = regexprep(sprintf('_target_%s', target),'[^a-zA-Z0-9_]','_');
                label = sprintf('juelich_prob_GM_Amygdala_laterobasal_group%s_thr%i_bin.nii.gz', target_lbl(1), threshold);
                target_mask = logical(niftiread(fullfile(outputs_path,'../../../', label)));
                s = regionprops3(target_mask,"Centroid", "BoundingBox");
                center_roi_slice = round(s.Centroid(1));

                create_group_MNI_plots(outputs_path, subject_list, parameters, 'plot_max_intensity', 1, 'outputs_suffix', sprintf('_max_focus_slice_roi_thresh_%i', threshold), 'ROI_MNI_mask', target_mask)
                create_group_MNI_plots(outputs_path, subject_list, parameters, 'plot_max_intensity', 0, 'slice_to_plot', center_roi_slice, 'outputs_suffix', sprintf('_roi_based_slice_roi_thresh_%i', threshold), 'ROI_MNI_mask', target_mask)

            end
        end

    end
end
            
