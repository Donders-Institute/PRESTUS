function [file, prof_csv] = calibration_save_entry(lib, fit, out_dir, options)
% CALIBRATION_SAVE_ENTRY  Write a fitted library entry to a user library folder
%
% Writes <out_dir>/<combo>.yaml, where combo is the transducer serial plus
% '_<driving system>' if the entry has one; transducer.library_path pointed
% at out_dir makes studies use it. The entry is marked user_empirical and
% gets meta.fit: how it was fitted (CALIBRATION_FIT record: mode, depths,
% profile CSV, reference ISPPA, fit and profile options, k-Wave drive check),
% so the fit can be reproduced and CALIBRATION_GUI / the PRESTUS Nexus can
% start the next fit from the same settings. When the fit used interpolated
% or extrapolated profiles, those profiles are written next to the entry
% (<combo>_calibration_profiles.csv, as REFINE_TRANSDUCER_LIBRARY does).
%
% Use as:
%   file = calibration_save_entry(lib, fit, out_dir)
%   [file, prof_csv] = calibration_save_entry(lib, fit, out_dir, options)
%
% Input:
%   lib     - entry to write (fit.entry, or the k-Wave-corrected one)
%   fit     - CALIBRATION_FIT result (.mode, .prof, .record)
%   out_dir - folder to write to (created if missing)
%   options - struct, all optional:
%     .ShippedLibraryPath - shipped library folder (default config/calibration);
%                           writing there errors unless .AllowShipped
%     .AllowShipped       - default false
%     .Kwave              - k-Wave drive check settings to record in meta.fit.kwave
%
% Output:
%   file     - path of the written YAML
%   prof_csv - path of the written calibration profiles, or ''
%
% See also: CALIBRATION_FIT, CALIBRATION_FIT_ENTRY, WRITE_AXIAL_PROFILES,
%           RESOLVE_TRANSDUCER_FROM_SERIAL

    if nargin < 4; options = struct(); end
    assert(~isempty(out_dir), 'calibration_save_entry:noFolder', 'No folder to save the entry to.');
    shipped = opt(options, 'ShippedLibraryPath', fullfile(get_prestus_path(), 'config', 'calibration'));
    if isfolder(out_dir) && isfolder(shipped) && strcmp(canonical(out_dir), canonical(shipped)) ...
            && ~opt(options, 'AllowShipped', false)
        error('calibration_save_entry:shipped', ['%s is the shipped library. Save to a folder of ' ...
            'your own and set transducer.library_path to it.'], out_dir);
    end

    combo = char(lib.meta.tran_serial);
    if isfield(lib.meta, 'driving_system') && ~isempty(lib.meta.driving_system)
        combo = [combo '_' char(lib.meta.driving_system)];
    end
    lib.meta.calibration_type = 'user_empirical';
    lib.meta.refinement_mode  = fit.mode;
    lib.meta.last_updated     = char(datetime('today', 'Format', 'yyyy-MM-dd'));
    try; lib.meta.prestus_version = prestus_version(); catch; lib.meta.prestus_version = 'unknown'; end
    rec = fit.record;
    kw  = opt(options, 'Kwave', []);
    if ~isempty(kw); rec.kwave = kw; end
    lib.meta.fit = rec;

    if ~isfolder(out_dir); mkdir(out_dir); end
    prof_csv = '';
    if any(~strcmp(fit.prof.profile_type, 'measured'))
        % The profiles the entry was fitted to, next to it (as REFINE_TRANSDUCER_LIBRARY)
        prof_csv = fullfile(out_dir, [combo '_calibration_profiles.csv']);
        write_axial_profiles(fit.prof, prof_csv, struct('combo', combo));
        if strcmp(fit.mode, 'per_depth'); lib.refinement.calibration_profiles_csv = prof_csv; end
    end
    file = fullfile(out_dir, [combo '.yaml']);
    yaml.dumpFile(file, lib);
end

function p = canonical(p)
    p = char(java.io.File(char(p)).getCanonicalPath());
end

function v = opt(s, f, d)
    if isstruct(s) && isfield(s, f) && ~isempty(s.(f)); v = s.(f); else; v = d; end
end
