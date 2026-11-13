function convert_final_to_MNI_simnibs(path_to_input_img, m2m_folder, path_to_output_img, parameters, options)
% CONVERT_FINAL_TO_MNI_SIMNIBS  Transform one or more images to MNI space
%
% Warps NIfTI image(s) from subject-specific (T1 conform) space to MNI space
% using either the SimNIBS subject2mni CLI tool or a native MATLAB
% implementation selected via options.method (default: 'simnibs').
%
% Native method ('native'):
%   Reads MNI2Conform_nonl.nii.gz from the m2m toMNI folder directly.
%   The field stores absolute world coordinates (mm) in subject conform
%   space for every voxel of the warp grid. The warp is loaded once and
%   reused across all images in a batch. Out-of-FOV voxels (outside the
%   warp grid or input image) are set to 0, then FillValues are applied.
%   Uses griddedInterpolant (base MATLAB, no toolbox required).
%   Equivalent to SimNIBS nonl transformation at order 0 (nearest) or 1
%   (linear); higher orders are not supported and silently fall back to 1.
%
%   Contrast with CONVERT_FINAL_TO_MNI_MATLAB, which applies the 12-DOF
%   affine registration (MNI2conform_12DOF.txt) to an in-memory array and
%   returns the transformed volume, composite affine, and MNI header.
%   Use that function when image data is already in memory or the affine
%   matrix is needed for downstream calculations. Use this function (either
%   method) when writing final result NIfTIs to disk via nonlinear warp.
%
%   Set parameters.analysis.mni_warp_method = 'native' to activate the
%   native path globally; PRESTUS callers propagate this via options.method.
%
% SimNIBS method ('simnibs'):
%   Calls the subject2mni CLI. An optional LD_LIBRARY_PATH export is
%   prepended when parameters.hpc.ld_library_path is set. Because SimNIBS
%   appends '_MNI' to the output filename, each result is renamed to the
%   requested path. Inputs are converted one after another in the foreground.
%   Each exit status and output file is checked, and any failure raises
%   PRESTUS:MNIConversionFailed listing the volumes and subject2mni output.
%
% Use as:
%   convert_final_to_MNI_simnibs(path_to_input_img, m2m_folder, path_to_output_img, parameters)
%   convert_final_to_MNI_simnibs({in1,in2,...}, m2m_folder, {out1,out2,...}, parameters)
%   convert_final_to_MNI_simnibs(..., 'interpolation_order', 0)
%   convert_final_to_MNI_simnibs(..., 'FillValues', vals)
%   convert_final_to_MNI_simnibs(..., 'method', 'native')
%
% Input:
%   path_to_input_img  - path or cell array of paths to input NIfTI(s) in subject space
%   m2m_folder         - path to SimNIBS m2m folder (must contain toMNI subdirectory)
%   path_to_output_img - desired output path(s) in MNI space (.nii.gz)
%   parameters         - PRESTUS config; used by SimNIBS method for
%                        startup.simnibs_bin_path and hpc.ld_library_path
%   options.interpolation_order - 0=nearest, 1=linear (default: 1)
%   options.FillValues - scalar or vector (one per input) fill value for
%                        out-of-FOV voxels (default: [] = no fill)
%   options.method     - 'simnibs' (default) or 'native'
%
% See also: SIMULATION_NIFTI, CONVERT_FINAL_TO_MNI_MATLAB, SUBJECT2MNI_COORDS_LDFIX

    arguments
        path_to_input_img
        m2m_folder          string
        path_to_output_img
        parameters          struct  = struct()
        options.interpolation_order = 1
        options.FillValues          = []
        options.method      string  = 'simnibs'
    end

    % Normalise scalar inputs to cell arrays so the rest of the code is uniform
    if ischar(path_to_input_img) || isstring(path_to_input_img)
        input_list  = {char(path_to_input_img)};
        output_list = {char(path_to_output_img)};
    else
        input_list  = path_to_input_img;
        output_list = path_to_output_img;
    end
    n = numel(input_list);

    switch lower(options.method)
        case 'native'
            convert_native(input_list, output_list, n, m2m_folder, options);
        case 'simnibs'
            convert_simnibs(input_list, output_list, n, m2m_folder, parameters, options);
        otherwise
            error('PRESTUS:BadMethod', ...
                'convert_final_to_MNI_simnibs: unknown method ''%s''. Use ''simnibs'' or ''native''.', ...
                options.method);
    end
end

% -------------------------------------------------------------------------
function convert_simnibs(input_list, output_list, n, m2m_folder, parameters, options)
    % run_command's Env option sets LD_LIBRARY_PATH natively on the child
    % process -- no `export ...;` shell syntax needed at all.
    env_overrides = struct();
    if isfield(parameters, 'hpc') && isfield(parameters.hpc, 'ld_library_path') && ...
            ~isempty(parameters.hpc.ld_library_path) && ~ispc
        env_overrides.LD_LIBRARY_PATH = parameters.hpc.ld_library_path;
    end
    subject2mni_bin = fullfile(parameters.startup.simnibs_bin_path, 'subject2mni');

    % SimNIBS always appends '_MNI' to the output stem
    simnibs_names = cell(1, n);
    for k = 1:n
        simnibs_names{k} = strrep(output_list{k}, '.nii.gz', '_MNI.nii.gz');
    end

    % Each conversion runs in the foreground and is checked. Background runs
    % discard their output and cannot report a failure, so a missing volume used
    % to end the job as "completed". run_command invokes subject2mni directly via
    % argv, so input/output/m2m paths reach the process as literal arguments.
    failed = {};
    for k = 1:n
        [status, cmdout] = run_command({subject2mni_bin, '--in', input_list{k}, ...
            '--out', output_list{k}, '--m2mpath', m2m_folder, ...
            '--interpolation_order', num2str(options.interpolation_order)}, ...
            'Env', env_overrides);
        if status ~= 0 || ~isfile(simnibs_names{k})
            failed{end+1} = sprintf('%s (exit %d): %s', input_list{k}, status, ...
                strtrim(char(cmdout))); %#ok<AGROW>
        end
    end
    if ~isempty(failed)
        error('PRESTUS:MNIConversionFailed', ...
            'subject2mni failed for %d of %d volumes:\n  %s', numel(failed), n, ...
            strjoin(failed, sprintf('\n  ')));
    end

    for k = 1:n
        if ~strcmp(simnibs_names{k}, output_list{k}) && isfile(simnibs_names{k})
            movefile(simnibs_names{k}, output_list{k});
        end
    end

    % Post-hoc fill: SimNIBS writes 0 for out-of-FOV voxels
    apply_fill_values(output_list, n, options.FillValues);
end

% -------------------------------------------------------------------------
function convert_native(input_list, output_list, n, m2m_folder, options)
% Native MATLAB nonlinear warp using MNI2Conform_nonl.nii.gz.
%
% MNI2Conform_nonl stores absolute world coordinates (mm) in subject
% conform space for each voxel of the warp grid (pull convention).
% Algorithm per output voxel:
%   MNI vox → MNI world (ref affine)
%   → warp-field vox (inv warp affine)
%   → interpolate warp → subject world coords
%   → subject vox (inv input affine)
%   → interpolate input image

    % --- Load warp field (shared across all images) ---
    warp_info = niftiinfo(fullfile(m2m_folder, 'toMNI', 'MNI2Conform_nonl.nii.gz'));
    warp_data = single(niftiread(warp_info));   % [Wx Wy Wz 3], subject world coords (mm)
    % niftiinfo's Transform.T maps 0-based voxel indices to world; the
    % 1-based MATLAB grids below are shifted by one on the way in and out.
    warp_A    = warp_info.Transform.T';         % 4×4: warp vox (0-based) → world

    % --- Reference output space: final_tissues_MNI defines the MNI grid ---
    ref_info = niftiinfo(fullfile(m2m_folder, 'toMNI', 'final_tissues_MNI.nii.gz'));
    ref_A    = ref_info.Transform.T';           % 4×4: MNI vox (0-based) → MNI world
    ref_dims = ref_info.ImageSize(1:3);

    % --- Build subject-world-coordinate lookup for every MNI output voxel ---
    [xi, yi, zi] = ndgrid(1:ref_dims(1), 1:ref_dims(2), 1:ref_dims(3));
    N = numel(xi);

    % MNI vox → MNI world
    mni_h = [xi(:)' - 1; yi(:)' - 1; zi(:)' - 1; ones(1, N, 'single')];
    clear xi yi zi
    mni_world = ref_A * double(mni_h);
    clear mni_h

    % MNI world → warp-field vox
    warp_vox = warp_A \ mni_world;             % 4×N, 0-based
    warp_vox(1:3, :) = warp_vox(1:3, :) + 1;
    clear mni_world

    % Interpolate warp field → subject world coords (NaN outside warp FOV)
    Wd  = size(warp_data, [1 2 3]);
    qx  = reshape(warp_vox(1,:), ref_dims);
    qy  = reshape(warp_vox(2,:), ref_dims);
    qz  = reshape(warp_vox(3,:), ref_dims);
    clear warp_vox

    subj_world = zeros([ref_dims, 3], 'single');
    for c = 1:3
        Fw = griddedInterpolant({1:Wd(1), 1:Wd(2), 1:Wd(3)}, ...
            warp_data(:,:,:,c), 'linear', 'none');
        subj_world(:,:,:,c) = Fw(qx, qy, qz);
    end
    clear warp_data qx qy qz Fw

    % Flatten for reuse across images; rows are [x;y;z;1] per MNI voxel
    sw_h = [reshape(subj_world, [], 3)'; ones(1, N)];   % 4×N (NaN where outside warp)
    clear subj_world

    % Clamp interpolation order: griddedInterpolant supports 0 (nearest) or 1 (linear)
    if options.interpolation_order == 0
        interp_method = 'nearest';
    else
        interp_method = 'linear';
    end

    % --- Process each image ---
    for k = 1:n
        restore = set_writable_tmpdir(input_list{k});
        in_info = niftiinfo(input_list{k});
        in_data = single(niftiread(in_info));
        in_A    = in_info.Transform.T';
        in_dims = in_info.ImageSize(1:3);

        % Subject world → input voxel coords
        in_vox = in_A \ sw_h;                  % 4×N, 0-based
        in_vox(1:3, :) = in_vox(1:3, :) + 1;

        px = reshape(in_vox(1,:), ref_dims);
        py = reshape(in_vox(2,:), ref_dims);
        pz = reshape(in_vox(3,:), ref_dims);

        % Interpolate; griddedInterpolant returns NaN for out-of-bounds queries
        Fi      = griddedInterpolant({1:in_dims(1), 1:in_dims(2), 1:in_dims(3)}, ...
                    in_data, interp_method, 'none');
        out_vol = Fi(px, py, pz);
        out_vol(isnan(out_vol)) = 0;            % match SimNIBS zero-fill convention

        % Round before cast for integer output types to avoid truncation artefacts
        src_dtype = in_info.Datatype;
        if ~isfloat(cast(0, src_dtype))
            out_vol = round(out_vol);
        end
        out_vol = cast(out_vol, src_dtype);

        % Write using ref-space header with input datatype
        out_info             = ref_info;
        out_info.Datatype    = src_dtype;
        out_info.BitsPerPixel = in_info.BitsPerPixel;
        niftiwrite(out_vol, strrep(output_list{k}, '.nii.gz', ''), out_info, 'Compressed', true);
        restore();
    end

    apply_fill_values(output_list, n, options.FillValues);
end

% -------------------------------------------------------------------------
function apply_fill_values(output_list, n, fill_values)
    if isempty(fill_values)
        return
    end
    fills = fill_values;
    if isscalar(fills)
        fills = repmat(fills, 1, n);
    end
    for k = 1:n
        fv = fills(k);
        if fv == 0 || ~isfile(output_list{k})
            continue
        end
        restore = set_writable_tmpdir(output_list{k});
        info = niftiinfo(output_list{k});
        vol  = niftiread(info);
        vol(vol == 0) = cast(fv, class(vol));
        niftiwrite(vol, strrep(output_list{k}, '.nii.gz', ''), info, 'Compressed', true);
        restore();
    end
end

% -------------------------------------------------------------------------
function restore_fn = set_writable_tmpdir(ref_file)
% Point TMPDIR at the directory containing ref_file so that niftiinfo can
% unzip .nii.gz files even when the system tempdir is not writable (e.g. HPC).
    prev = getenv('TMPDIR');
    candidate = fileparts(ref_file);
    if ~isempty(candidate) && isfolder(candidate)
        setenv('TMPDIR', candidate);
    end
    restore_fn = @() setenv('TMPDIR', prev);
end
