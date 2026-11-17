function sidecar_file = write_sidecar(sidecar_file, kind, content, parameters)
% WRITE_SIDECAR  Write a JSON sidecar describing a stage output later runs can start from
%
% Sidecars describe complete stage outputs that a later, separate run can
% start from (a placement run, a pseudo-CT, an acoustic pressure map). They
% are distinct from the within-run caches in cache/, which only the run that
% wrote them (or its re-runs) reloads. As in BIDS, a sidecar is a JSON file
% next to the data. Every sidecar starts with a common header, followed by
% the fields of CONTENT:
%
%   "sidecar_version": 1
%   "kind":            "placement" | "pseudoCT" | "acoustic_pressure"
%   "created":         ISO 8601 time
%   "generated_by":    {"name": "PRESTUS", "version": <git describe or VERSION>,
%                       "commit": <short hash>, "commit_date": <yyyy-mm-dd>,
%                       "dirty": "yes" | "no" | "unknown"}
%   "subject_id":      <parameters.subject_id>
%   ...                kind-specific fields from CONTENT
%
% A CONTENT field 'signature' holds the settings the output depends on; a
% reusing run compares it with its own (CHECK_SIDECAR_SIGNATURE). Never
% errors: a failure to write is reported as a warning and returns ''.
%
% Use as:
%   sidecar_file = write_sidecar(sidecar_file, kind, content, parameters)
%
% Input:
%   sidecar_file - path of the .json file to write
%   kind         - sidecar kind (char)
%   content      - (1,1) struct of kind-specific fields
%   parameters   - PRESTUS parameters (subject_id)
%
% Output:
%   sidecar_file - path written, or '' when nothing was written
%
% See also: READ_SIDECAR, FIND_SIDECAR, CHECK_SIDECAR_SIGNATURE

arguments
    sidecar_file {mustBeTextScalar}
    kind         {mustBeTextScalar}
    content      (1,1) struct
    parameters   (1,1) struct
end

    sidecar_file = char(sidecar_file);
    try
        % The run's own record (path_log_setup) when available, so the
        % sidecar names the code that produced the output.
        if isfield(parameters, 'io') && isfield(parameters.io, 'prestus_version') && ...
                isstruct(parameters.io.prestus_version)
            info = parameters.io.prestus_version;
        else
            info = struct('version', 'unknown', 'commit', 'unknown', 'commit_date', '', ...
                'branch', '', 'dirty', 'unknown');
            try
                [~, ~, info] = prestus_version();
            catch
            end
        end
        out = struct( ...
            'sidecar_version', 1, ...
            'kind',            char(kind), ...
            'created',         char(datetime('now', 'Format', 'yyyy-MM-dd''T''HH:mm:ss')), ...
            'generated_by',    struct('name', 'PRESTUS', 'version', info.version, ...
                                      'commit', info.commit, 'commit_date', info.commit_date, ...
                                      'dirty', info.dirty), ...
            'subject_id',      parameters.subject_id);
        fields = fieldnames(content);
        for i = 1:numel(fields)
            out.(fields{i}) = content.(fields{i});
        end
        folder = fileparts(sidecar_file);
        if ~isempty(folder) && ~isfolder(folder)
            mkdir(folder);
        end
        text = jsonencode(out, 'PrettyPrint', true);
        fid = fopen(sidecar_file, 'w', 'n', 'UTF-8');
        if fid < 0
            error('cannot open %s for writing', sidecar_file);
        end
        closer = onCleanup(@() fclose(fid));
        fprintf(fid, '%s\n', text);
    catch ME
        warning('prestus:sidecar:notWritten', 'Could not write the sidecar %s: %s', ...
            sidecar_file, ME.message);
        sidecar_file = '';
    end
end
