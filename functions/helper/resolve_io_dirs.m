function p = resolve_io_dirs(p)
% RESOLVE_IO_DIRS  Fill missing io path fields from path.sim + subject_id
%
% Mirrors the path-derivation logic in path_log_setup without any side
% effects: no directories are created, no logging is started.
%
% Call this at the top of any report generator that receives parameter
% structs built before the pipeline has run (e.g. params stored for
% sequential or uncertainty report generation).  Fields that are already
% set on the struct are left unchanged.
%
% Use as:
%   p = resolve_io_dirs(p)
%
% See also: PATH_LOG_SETUP, GET_OUTPUT_DIR

arguments
    p (1,1) struct
end

    % Derive dir_output from path.sim when absent
    if ~isfield(p, 'io') || ~isfield(p.io, 'dir_output') || isempty(p.io.dir_output)
        try
            p.io.dir_output = get_output_dir(p);
        catch
            return   % path.sim missing — nothing more we can do
        end
    end

    out = p.io.dir_output;

    % Subdirectory defaults — only fill gaps, never overwrite
    dirs = { ...
        'dir_nii',            fullfile(out, 'nii');    ...
        'dir_nii_T1w',        fullfile(out, 'nii');    ...
        'dir_nii_MNI',        fullfile(out, 'nii');    ...
        'dir_img',            fullfile(out, 'img');    ...
        'dir_tabular',        out;                      ...
        'dir_reports',        out;                      ...
        'dir_logs',           fullfile(out, 'log');    ...
        'dir_cache',          fullfile(out, 'cache');  ...
        'dir_debug',          fullfile(out, 'debug');  ...
    };
    for k = 1:size(dirs, 1)
        fname = dirs{k, 1};
        if ~isfield(p.io, fname) || isempty(p.io.(fname))
            p.io.(fname) = dirs{k, 2};
        end
    end

    % filename_table (CSV output path)
    if ~isfield(p.io, 'filename_table') || isempty(p.io.filename_table)
        if isfield(p.io, 'output_affix') && isfield(p, 'subject_id') && ...
                isfield(p, 'simulation') && isfield(p.simulation, 'medium')
            p.io.filename_table = fullfile(out, ...
                sprintf('sub-%03d_%s%s.csv', ...
                    p.subject_id, p.simulation.medium, p.io.output_affix));
        end
    end

end
