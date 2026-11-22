function nifti_transducer(parameters, planimg)
% NIFTI_TRANSDUCER  Export a model of the planned transducer in T1 space
%
% Writes the planned transducer bowl(s) as a NIfTI on the voxel grid of the
% planning T1 (planimg.t1_header), the header every T1w output carries, with
% value 1 on the emitting face and 0.5 on the display thickness behind it
% (WRITE_TRANSDUCER_NIFTI), so
% it overlays on the T1 and the acoustic maps for placement checks (e.g. in
% FSLeyes). The positions are the planned T1-voxel positions that
% PRESTUS_PIPELINE records before the grid stage
% (parameters.state.planned_trans_t1_vox / planned_focus_t1_vox), the same
% positions TARGETING_ANALYSIS uses. trans_pos is the rear centre of the bowl.
%
% Written for layered media with annular transducers, when
% io.save_transducer_nii is set and NIfTI creation is enabled. Matrix arrays
% are skipped: TRANSDUCER_SETUP aligns them with the grid z axis, which holds
% only in the simulation grid. No MNI copy is written, because the transducer
% lies outside the head that the SimNIBS warp covers.
%
% Output: <io.dir_nii_T1w>/sub-XXX_<medium>_T1w<affix>_transducer.nii.gz
%
% Use as:
%   nifti_transducer(parameters, planimg)
%
% Input:
%   parameters - PRESTUS config struct (after the grid stage); reads
%                io.save_transducer_nii and io.transducer_nii_back_mm
%   planimg    - struct with field t1_header (planning T1, RAS+)
%
% See also: WRITE_TRANSDUCER_NIFTI, TRANSDUCER_SETUP, NIFTI_MEDIUM, TARGETING_ANALYSIS

arguments
    parameters (1,1) struct
    planimg    (1,1) struct
end

    if isfield(parameters.modules, 'run_nifti_creation') && ...
            parameters.modules.run_nifti_creation == 0
        return
    end
    if ~strcmp(parameters.simulation.medium, 'layered')
        return
    end
    if isfield(parameters.grid, 'axisymmetric') && parameters.grid.axisymmetric == 1
        return
    end
    if ~should_save_output(parameters.io, 'save_transducer_nii')
        return
    end
    if ~isfield(planimg, 't1_header') || isempty(planimg.t1_header)
        warn('nifti_transducer: no planning T1 header; transducer NIfTI not written.');
        return
    end

    % Planned positions in T1 voxels (PRESTUS_PIPELINE, after placement;
    % re-recorded in the RAS+ T1 grid by GRID_TISSUE_SETUP)
    tpos = [];
    fpos = [];
    if isfield(parameters, 'state') && isfield(parameters.state, 'planned_trans_t1_vox') && ...
            isfield(parameters.state, 'planned_focus_t1_vox')
        tpos = parameters.state.planned_trans_t1_vox;
        fpos = parameters.state.planned_focus_t1_vox;
    end
    n = numel(parameters.transducer);
    if isempty(tpos) || isempty(fpos) || size(tpos, 1) ~= n || size(fpos, 1) ~= n
        warn('nifti_transducer: planned T1 positions unavailable; transducer NIfTI not written.');
        return
    end
    is_annular = arrayfun(@(t) strcmp(t.type, 'annular'), parameters.transducer);
    if ~all(is_annular)
        warn('nifti_transducer: only annular transducers are supported; transducer NIfTI not written.');
        return
    end

    out_file = fullfile(parameters.io.dir_nii_T1w, sprintf('sub-%03d_%s_T1w%s_transducer.nii.gz', ...
        parameters.subject_id, parameters.simulation.medium, parameters.io.output_affix));
    if ~confirm_overwriting(out_file, parameters)
        return
    end

    back_mm = 3;
    if isfield(parameters.io, 'transducer_nii_back_mm') && ~isempty(parameters.io.transducer_nii_back_mm)
        back_mm = double(parameters.io.transducer_nii_back_mm);
    end

    % WRITE_TRANSDUCER_NIFTI reads the positions from parameters.transducer,
    % which after the grid stage holds simulation-grid indices
    params_t1 = parameters;
    for ti = 1:n
        params_t1.transducer(ti).trans_pos = tpos(ti, :);
        params_t1.transducer(ti).focus_pos = fpos(ti, :);
    end
    write_transducer_nifti(params_t1, planimg.t1_header, out_file, back_mm);
end
