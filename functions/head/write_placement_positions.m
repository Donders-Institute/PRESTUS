function json_path = write_placement_positions(parameters, mode)
% WRITE_PLACEMENT_POSITIONS  Write the resolved transducer positions as a JSON sidecar
%
% Called at the end of the placement stage for every placement.mode, so that
% tools (e.g. the PRESTUS Nexus placement viewer) can show where a heuristic,
% Localite, MNI, PlanTUS or TUS_entry placement put the transducers without
% re-running it. Positions are 1-based T1 voxel indices, as in
% transducer(i).trans_pos / focus_pos after placement.
%
% File: <io.dir_output>/sub-XXX_placement<affix>.json
%   { "subject_id": 1, "mode": "heuristic", "space": "t1_voxel_1based",
%     "created": "2026-10-02T12:00:00",
%     "transducers": [ {"index": 1, "trans_pos": [x,y,z], "focus_pos": [x,y,z]} ] }
%
% Never errors: a failure to write is reported as a warning.
%
% Use as:
%   json_path = write_placement_positions(parameters, mode)
%
% Input:
%   parameters - PRESTUS parameters after placement (io.dir_output set)
%   mode       - placement mode used (char)
%
% Output:
%   json_path  - path written, or '' when nothing was written
%
% See also: PREPROC_TRANSDUCER_PLACEMENT, GENERATE_PLACEMENT_REPORT

    json_path = '';
    try
        if ~isfield(parameters, 'io') || ~isfield(parameters.io, 'dir_output') || ...
                isempty(parameters.io.dir_output)
            return
        end
        affix = '';
        if isfield(parameters.io, 'output_affix') && ischar(parameters.io.output_affix)
            affix = parameters.io.output_affix;
        end
        out = struct('subject_id', parameters.subject_id, 'mode', char(mode), ...
            'space', 't1_voxel_1based', ...
            'created', char(datetime('now', 'Format', 'yyyy-MM-dd''T''HH:mm:ss')));
        tr = struct('index', {}, 'trans_pos', {}, 'focus_pos', {});
        for ti = 1:numel(parameters.transducer)
            t = parameters.transducer(ti);
            tr(end+1) = struct('index', ti, ...
                'trans_pos', position_or_empty(t, 'trans_pos'), ...
                'focus_pos', position_or_empty(t, 'focus_pos')); %#ok<AGROW>
        end
        out.transducers = num2cell(tr);   % cell: a JSON list even for one transducer
        if ~isfolder(parameters.io.dir_output)
            mkdir(parameters.io.dir_output);
        end
        json_path = fullfile(parameters.io.dir_output, ...
            sprintf('sub-%03d_placement%s.json', parameters.subject_id, affix));
        fid = fopen(json_path, 'w');
        if fid < 0
            error('cannot open %s', json_path);
        end
        cleaner = onCleanup(@() fclose(fid));
        fprintf(fid, '%s\n', jsonencode(out, 'PrettyPrint', true));
    catch ME
        warning('prestus:placement:positionsNotWritten', ...
            'Could not write the placement positions JSON: %s', ME.message);
        json_path = '';
    end
end

function v = position_or_empty(t, field)
    v = [];
    if isfield(t, field) && isnumeric(t.(field)) && numel(t.(field)) == 3
        v = double(t.(field)(:)');
    end
end
