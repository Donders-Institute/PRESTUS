function record_file = write_placement_record(parameters, mode, run_dir, signature)
% WRITE_PLACEMENT_RECORD  Write the run record (placement_run.json) of a placement run
%
% A sidecar (WRITE_SIDECAR, kind 'placement') in the placement run folder
% holding the resolved positions (T1 voxels, 1-based; RAS mm when known), the
% target, and the settings they depend on (PLACEMENT_SIGNATURE). A later run
% reuses the placement from it according to placement.reuse
% (PREPROC_TRANSDUCER_PLACEMENT); other tools can read it with
% READ_PLACEMENT_RECORD.
%
% Use as:
%   record_file = write_placement_record(parameters, mode, run_dir, signature)
%
% Input:
%   parameters - PRESTUS parameters after placement (transducer(i).trans_pos ...)
%   mode       - placement mode ('heuristic' | 'plantus' | 'tus_entry')
%   run_dir    - placement run folder (PLACEMENT_RUN_DIR)
%   signature  - PLACEMENT_SIGNATURE, taken before the placement ran
%
% Output:
%   record_file - path written, or '' when writing failed (warning)
%
% See also: READ_PLACEMENT_RECORD, PLACEMENT_SIGNATURE, WRITE_SIDECAR

arguments
    parameters (1,1) struct
    mode       {mustBeTextScalar}
    run_dir    {mustBeTextScalar}
    signature  (1,1) struct
end

    mode = char(mode);
    record_file = '';
    if ~all(arrayfun(@(t) numel(t.trans_pos) == 3 && numel(t.focus_pos) == 3, parameters.transducer))
        warning('prestus:placement:noRecord', ...
            'Placement %s resolved no complete positions; no run record written.', mode);
        return
    end
    tr = cell(1, numel(parameters.transducer));
    for ti = 1:numel(parameters.transducer)
        t = parameters.transducer(ti);
        e = struct('index', ti, ...
            'trans_pos', double(t.trans_pos(:)'), 'focus_pos', double(t.focus_pos(:)'));
        if isfield(t, 'trans_pos_ras') && numel(t.trans_pos_ras) == 3 && ...
                isfield(t, 'focus_pos_ras') && numel(t.focus_pos_ras) == 3
            e.trans_pos_ras = double(t.trans_pos_ras(:)');
            e.focus_pos_ras = double(t.focus_pos_ras(:)');
        end
        tr{ti} = e;
    end
    content = struct('mode', mode, ...
        'space', struct('trans_pos', 't1_voxel_1based', 'trans_pos_ras', 'ras_mm'), ...
        'transducers', {tr}, 'target_ras', []);
    if isfield(parameters.placement, mode) && isfield(parameters.placement.(mode), 'target_ras')
        content.target_ras = double(parameters.placement.(mode).target_ras(:)');
    end
    content.signature = signature;
    record_file = write_sidecar(fullfile(char(run_dir), 'placement_run.json'), ...
        'placement', content, parameters);
end
