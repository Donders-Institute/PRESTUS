function defaced_folder = preproc_deface_m2m(segmentation_folder, parameters)
% PREPROC_DEFACE_M2M  Create a defaced copy of a SimNIBS m2m folder
%
% Copies the m2m folder to <m2m_folder>_defaced/, then runs pydeface on
% every anatomical NIfTI in the copy. Non-anatomical files (tissue labels,
% surfaces, warps, meshes) are copied unchanged. The original folder is
% never modified.
%
% A sentinel file (.defaced) is written on success. If it already exists
% the function returns immediately so repeated pipeline runs are cheap.
%
% Use as:
%   defaced_folder = preproc_deface_m2m(segmentation_folder, parameters)
%
% Input:
%   segmentation_folder - full path to an existing m2m_sub-NNN/ folder
%   parameters          - (1,1) PRESTUS parameters struct; must contain
%                           parameters.startup.simnibs_bin_path
%
% Output:
%   defaced_folder - path to the created m2m_sub-NNN_defaced/ folder,
%                    or '' if defacing was skipped or failed
%
% See also: PREPROC_SEGMENTATION, SEGMENTATION_RUN

    arguments
        segmentation_folder (1,:) char
        parameters          (1,1) struct
    end

    defaced_folder = '';

    %% Resolve pydeface binary
    if ~isfield(parameters, 'startup') || ~isfield(parameters.startup, 'simnibs_bin_path') ...
            || isempty(parameters.startup.simnibs_bin_path)
        warning('preproc_deface_m2m: parameters.startup.simnibs_bin_path not set; skipping defacing.');
        return;
    end
    pydeface_bin = fullfile(parameters.startup.simnibs_bin_path, 'pydeface');
    if ~isfile(pydeface_bin)
        warning('preproc_deface_m2m: pydeface not found at %s; skipping defacing.', pydeface_bin);
        return;
    end

    %% Paths
    defaced_folder = [segmentation_folder, '_defaced'];
    sentinel       = fullfile(defaced_folder, '.defaced');

    if isfile(sentinel)
        fprintf('Defaced m2m folder already exists, skipping:\n  %s\n', defaced_folder);
        return;
    end

    %% Copy entire m2m folder
    fprintf('Copying m2m folder for defacing...\n  %s\n  -> %s\n', ...
        segmentation_folder, defaced_folder);
    if isfolder(defaced_folder)
        % Partial copy from a previous interrupted run - remove and restart
        rmdir(defaced_folder, 's');
    end
    [status, out] = system(sprintf('cp -r "%s" "%s"', segmentation_folder, defaced_folder));
    if status ~= 0
        warning('preproc_deface_m2m: cp failed:\n%s', out);
        defaced_folder = '';
        return;
    end

    %% Anatomical NIfTIs to deface (relative paths within the m2m folder)
    % Only T1/T2 derived images - label maps and warp fields are not identifiable.
    anat_files = { ...
        'T1.nii.gz'; ...
        'UTE_reg.nii.gz'; ...
        fullfile('label_prep',   'T1_upsampled.nii.gz'); ...
        fullfile('label_prep',   'T2_upsampled.nii.gz'); ...
        fullfile('segmentation', 'T1_bias_corrected.nii.gz'); ...
        fullfile('segmentation', 'T2_bias_corrected.nii.gz'); ...
        fullfile('toMNI',        'T1_to_MNI_post-hoc.nii.gz'); ...
    };

    failed = {};
    for i = 1:numel(anat_files)
        fpath = fullfile(defaced_folder, anat_files{i});
        if ~isfile(fpath)
            continue;
        end
        fprintf('  Defacing %s...\n', anat_files{i});
        cmd = sprintf('"%s" "%s" --outfile "%s" --force 2>&1', pydeface_bin, fpath, fpath);
        [status, out] = system(cmd);
        if status ~= 0
            warning('preproc_deface_m2m: pydeface failed for %s:\n%s', anat_files{i}, out);
            failed{end+1} = anat_files{i}; %#ok<AGROW>
        end
    end

    if ~isempty(failed)
        warning('preproc_deface_m2m: defacing incomplete (%d file(s) failed). Defaced folder may not be suitable for sharing.', numel(failed));
    end

    %% Write sentinel
    fclose(fopen(sentinel, 'w'));
    fprintf('Defaced m2m folder ready:\n  %s\n', defaced_folder);
end
