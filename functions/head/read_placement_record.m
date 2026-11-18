function [rec, record_file] = read_placement_record(run_dir)
% READ_PLACEMENT_RECORD  Read the run record of a placement run folder
%
% Reads placement_run.json (WRITE_PLACEMENT_RECORD), or else a
% placement_run.mat holding the record fields as variables (loaded as is),
% from a placement run folder (PLACEMENT_RUN_DIR) into one struct with N×3
% position matrices.
%
% Use as:
%   [rec, record_file] = read_placement_record(run_dir)
%
% Input:
%   run_dir - placement run folder
%
% Output:
%   rec         - struct: signature, trans_pos, focus_pos (T1 voxels, 1-based),
%                 trans_pos_ras, focus_pos_ras (RAS mm, [] when unknown),
%                 target_ras; [] when the folder holds no record
%   record_file - the record read (for messages)
%
% See also: WRITE_PLACEMENT_RECORD, PREPROC_TRANSDUCER_PLACEMENT, PLACEMENT_RUN_DIR

arguments
    run_dir {mustBeTextScalar}
end

    rec = [];
    record_file = fullfile(char(run_dir), 'placement_run.json');
    if isfile(record_file)
        sc = read_sidecar(record_file);
        rec = struct('signature', [], 'trans_pos', [], 'focus_pos', [], ...
            'trans_pos_ras', [], 'focus_pos_ras', [], 'target_ras', []);
        if isfield(sc, 'signature'); rec.signature = sc.signature; end
        if isfield(sc, 'target_ras'); rec.target_ras = as_row(sc.target_ras); end
        tr = sc.transducers;              % struct array, or cell when fields differ
        if ~iscell(tr); tr = num2cell(tr); end
        for ti = 1:numel(tr)
            for f = {'trans_pos', 'focus_pos', 'trans_pos_ras', 'focus_pos_ras'}
                if isfield(tr{ti}, f{1}) && ~isempty(as_row(tr{ti}.(f{1})))
                    rec.(f{1})(ti, :) = as_row(tr{ti}.(f{1}));
                end
            end
        end
        return
    end
    legacy = fullfile(char(run_dir), 'placement_run.mat');
    if isfile(legacy)
        record_file = legacy;
        rec = load(legacy);
    end
end

function v = as_row(v)
% Numeric row from a JSON array (decoded as column or cell); [] for empty/null.
    if iscell(v)
        if isempty(v); v = []; return; end
        v = cellfun(@double, v);
    end
    if isempty(v)
        v = [];
        return
    end
    v = double(v(:)');
end
