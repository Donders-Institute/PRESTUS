function write_axial_profiles(prof, csv_path, meta)
% WRITE_AXIAL_PROFILES  Write axial profiles as a profile CSV plus YAML sidecar
%
% Writes the layout READ_AXIAL_PROFILES reads (doc_transducer_library.md,
% Appendix): a title row, a row with the focal depths [mm from EP], then one
% row per axial position (distance from the exit plane [mm], intensity
% [W/cm²] per depth). Columns that are entirely NaN ('unavailable') are
% skipped. Next to it, {stem}.yaml records what the columns are, including
% for interpolated/extrapolated columns which measured depths they were
% built from (INTERPOLATE_AXIAL_PROFILES) so the profile a calibration was
% fitted to can be traced back to the measurement.
%
% Use as:
%   write_axial_profiles(prof, csv_path)
%   write_axial_profiles(prof, csv_path, meta)
%
% Input:
%   prof     - profile struct (READ_AXIAL_PROFILES or INTERPOLATE_AXIAL_PROFILES)
%   csv_path - output .csv; the sidecar is written next to it as .yaml
%   meta     - (optional) struct merged into the sidecar, e.g.
%              struct('tran_serial', ..., 'ds_serial', ...)
%
% See also: READ_AXIAL_PROFILES, INTERPOLATE_AXIAL_PROFILES

    if nargin < 3; meta = struct(); end
    keep = any(isfinite(prof.I), 1);
    d = prof.depths_ep_mm(keep);
    I = prof.I(:, keep);
    z = prof.z_ep_mm(:);
    if isempty(d); error('write_axial_profiles: no profile to write (all columns unavailable).'); end

    folder = fileparts(csv_path);
    if ~isempty(folder) && ~isfolder(folder); mkdir(folder); end
    fid = fopen(csv_path, 'w');
    if fid < 0; error('write_axial_profiles: cannot write %s.', csv_path); end
    closer = onCleanup(@() fclose(fid));
    fprintf(fid, ',Focal depths (mm)%s\n', repmat(',', 1, numel(d) - 1));
    fprintf(fid, 'Distance WRT Exit plane (mm)%s\n', sprintf(',%.10g', d));
    fprintf(fid, ['%.10g' repmat(',%.10g', 1, numel(d)) '\n'], [z, I]');
    clear closer

    % Sidecar
    [~, stem] = fileparts(csv_path);
    types = repmat({'measured'}, 1, numel(prof.depths_ep_mm));
    if isfield(prof, 'profile_type'); types = prof.profile_type; end
    types = types(keep);
    sc = struct();
    fn = fieldnames(meta);
    for k = 1:numel(fn); sc.(fn{k}) = meta.(fn{k}); end
    sc.quantity           = 'intensity';
    sc.units              = 'W/cm2';
    sc.distance_reference = 'exit_plane';
    sc.focal_label_is     = 'nominal';
    sc.depths_ep_mm       = d;
    sc.profile_type       = types;
    if any(~strcmp(types, 'measured'))
        cols = find(keep);
        sc.derived = struct('depth_ep_mm', {}, 'profile_type', {}, 'source_depths_ep_mm', {}, ...
                            'weights', {}, 'shift_mm', {});
        for k = find(~strcmp(types, 'measured'))
            c = cols(k); ok = isfinite(prof.source_depths(c, :));
            sc.derived(end+1) = struct('depth_ep_mm', d(k), 'profile_type', types{k}, ...
                'source_depths_ep_mm', prof.source_depths(c, ok), 'weights', prof.weights(c, ok), ...
                'shift_mm', prof.shift_mm(c, ok));
        end
        sc.derived_method = sprintf(['focus-aligned linear blend of the neighbouring measured profiles ' ...
            '(align: %s), interpolate_axial_profiles'], prof.align);
    end
    src = '';
    if isfield(prof, 'source'); src = char(prof.source); end
    sc.provenance = struct('source_csv', src, 'written_by', 'write_axial_profiles', ...
        'written_at', datestr(now, 'yyyy-mm-dd HH:MM')); %#ok<TNOW1,DATST>
    yaml.dumpFile(fullfile(folder, [stem '.yaml']), sc);
end
