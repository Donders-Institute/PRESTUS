function nifti_medium(parameters, planimg, medium_masks, kwave_medium, pseudoCT)
% NIFTI_MEDIUM  Export medium-stage data to NIfTI
%
% Coordinates four distinct export steps:
%   1. Pipeline outputs  — medium_masks (see NIFTI_MEDIUM_MASKS)
%   2. Pipeline outputs  — pseudoCT, only when parameters.pct.enabled == 1 (see NIFTI_PCT)
%   3. Property maps     — per-tissue acoustic/thermal properties (see NIFTI_MEDIUM_PROPERTIES)
%                          only written when parameters.io.save_property_maps is true
%   4. T1-in-MNI         — registers the subject T1 to MNI if not yet done (see SIMNIBS_T1_TO_MNI)
%
% Use as:
%   nifti_medium(parameters, planimg, medium_masks, kwave_medium)
%   nifti_medium(parameters, planimg, medium_masks, kwave_medium, pseudoCT)
%
% Input:
%   parameters   - PRESTUS config struct
%   planimg      - struct with fields: t1_image_orig, inv_transf, t1_header, transf
%   medium_masks - (:,:,:) uint8, simulation-grid tissue labels
%   kwave_medium - struct with medium property fields
%   pseudoCT     - (:,:,:) numeric, simulation-grid Hounsfield values (optional)
%
% See also: NIFTI_MEDIUM_MASKS, NIFTI_PCT, NIFTI_MEDIUM_PROPERTIES, SIMNIBS_T1_TO_MNI, NIFTI_ACOUSTIC, NIFTI_THERMAL

arguments
    parameters   (1,1) struct
    planimg      (1,1) struct
    medium_masks
    kwave_medium (1,1) struct
    pseudoCT               = []
end

    if isfield(parameters.modules, 'run_nifti_creation') && ...
            parameters.modules.run_nifti_creation == 0
        disp('No nifti creation requested...')
        return
    end

    if ~contains(parameters.simulation.medium, {'layered'; 'phantom'})
        return
    end

    if isfield(parameters.grid, 'axisymmetric') && parameters.grid.axisymmetric == 1
        disp('Axisymmetric mode: skipping medium NIfTI export.')
        return
    end

    is_layered = strcmp(parameters.simulation.medium, 'layered');
    m2m_folder = fullfile(parameters.path.seg, sprintf('m2m_sub-%03d', parameters.subject_id));

    nifti_medium_masks(parameters, planimg, medium_masks, is_layered, m2m_folder);

    if is_layered && isfield(parameters, 'pct') && isfield(parameters.pct, 'enabled') && ...
            parameters.pct.enabled == 1 && ~isempty(pseudoCT)
        nifti_pct(parameters, planimg, pseudoCT, m2m_folder);
    end

    if is_layered && should_save_output(parameters.io, 'save_property_maps')
        nifti_medium_properties(parameters, planimg, kwave_medium, m2m_folder);
    end

    simnibs_t1_to_mni(parameters, is_layered, m2m_folder);
end
