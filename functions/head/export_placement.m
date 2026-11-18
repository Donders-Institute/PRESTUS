function files = export_placement(parameters)
% EXPORT_PLACEMENT  Write the placed transducer poses for other tools (placement.export_formats)
%
% Called at the end of the placement stage for every placement.mode. For each
% transducer and each format in placement.export_formats it writes, to
% io.dir_output:
%   'brainsight'  sub-XXX_placement<affix>[_txN]_brainsight.txt   (EXPORT_BRAINSIGHT_TRAJECTORY)
%   'slicer'      sub-XXX_placement<affix>[_txN]_slicer.tfm        (EXPORT_SLICER_TRANSFORM)
%   'kplan'       sub-XXX_placement<affix>[_txN]_kplan/            (EXPORT_KPLAN_PLAN): T1 and
%                 CT or pseudo-CT with the rotation removed, .kps, entry sheet
% The suffix _txN is added with several transducers. Never errors: a format
% that cannot be written is reported as a warning.
%
% Use as:
%   files = export_placement(parameters)
%
% Input:
%   parameters - PRESTUS parameters after placement (transducer(i).trans_pos /
%                focus_pos as 1-based T1 voxels; io.dir_output set)
%
% Output:
%   files - cell array of written paths (folders for 'kplan')
%
% See also: PREPROC_TRANSDUCER_PLACEMENT, WRITE_PLACEMENT_POSITIONS

    files = {};
    formats = get_field(parameters, {'placement', 'export_formats'}, {});
    if isempty(formats), return, end
    formats = lower(cellstr(string(formats)));
    known = {'brainsight', 'slicer', 'kplan'};
    for f = setdiff(formats, known)
        warning('prestus:placement:exportFormat', ...
            'placement.export_formats: unknown format ''%s'' (use %s).', f{1}, strjoin(known, ', '));
    end
    formats = intersect(known, formats, 'stable');
    if isempty(formats) || isempty(get_field(parameters, {'io', 'dir_output'}, ''))
        return
    end

    try
        t1_file = fullfile(parameters.path.anat, sprintf(parameters.path.t1_pattern, parameters.subject_id));
        T = niftiinfo(t1_file).Transform.T;            % 0-based voxel (row) -> RAS
        to_ras = @(v) [double(v(:)') - 1, 1] * T * [eye(3); 0 0 0];
    catch ME
        warning('prestus:placement:exportT1', 'Placement export skipped: %s', ME.message);
        return
    end
    affix = char(string(get_field(parameters, {'io', 'output_affix'}, '')));
    n_tr = numel(parameters.transducer);

    for ti = 1:n_tr
        tr = parameters.transducer(ti);
        if ~isfield(tr, 'trans_pos') || ~isfield(tr, 'focus_pos') || ...
                numel(tr.trans_pos) ~= 3 || numel(tr.focus_pos) ~= 3
            continue
        end
        apex = to_ras(tr.trans_pos);  focus = to_ras(tr.focus_pos);
        base = sprintf('sub-%03d_placement%s', parameters.subject_id, affix);
        if n_tr > 1, base = sprintf('%s_tx%d', base, ti); end
        stem = fullfile(parameters.io.dir_output, base);
        for f = formats
            try
                switch f{1}
                    case 'brainsight'
                        files{end+1} = export_brainsight_trajectory(apex, focus, base, [stem '_brainsight.txt'], ...
                            t1_file=t1_file); %#ok<AGROW>
                    case 'slicer'
                        files{end+1} = export_slicer_transform(apex, focus, [stem '_slicer.tfm']); %#ok<AGROW>
                    case 'kplan'
                        [offset, kname] = kplan_device(tr);
                        export_kplan_plan(t1_file, apex, focus, [stem '_kplan'], ...
                            sprintf('sub-%03d', parameters.subject_id), ct_file=planning_ct(parameters), ...
                            kplan_offset_mm=offset, transducer_name=kname);
                        files{end+1} = [stem '_kplan']; %#ok<AGROW>
                end
            catch ME
                warning('prestus:placement:exportFailed', 'Placement export ''%s'' (transducer %d) failed: %s', ...
                    f{1}, ti, ME.message);
            end
        end
    end
    if ~isempty(files)
        fprintf('[placement] exported: %s\n', strjoin(files, ', '));
    end
end

function ct = planning_ct(parameters)
% The measured CT (path.ct_pattern), else the pseudo-CT when it exists; '' otherwise.
    ct = '';
    try
        ct = pct_ct_file(parameters);
    catch
    end
    d = char(string(get_field(parameters, {'io', 'dir_pct'}, '')));
    if isempty(ct) && ~isempty(d) && isfile(fullfile(d, 'pseudoCT.nii.gz'))
        ct = fullfile(d, 'pseudoCT.nii.gz');
    end
end

function [offset, name] = kplan_device(tr)
% k-Plan library devices with a known apex -> exit plane distance (from k-Plan's
% device-parameter export of plan 69-928). Others: the entry sheet explains it.
    offset = [];  name = '';
    model = '';
    if isfield(tr, 'model'), model = char(string(tr.model)); end
    if isempty(model) && isfield(tr, 'serial'), model = char(string(tr.serial)); end
    if startsWith(model, 'CTX_500')
        offset = 10.82;  name = 'NEUROFUS-CTX-500-4-026 (k-Plan library)';
    end
end
