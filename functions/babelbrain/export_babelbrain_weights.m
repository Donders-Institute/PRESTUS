function export_babelbrain_weights(phases_rad, velocity, parameters, out_path, geo_phases_rad, elem_amp)
% EXPORT_BABELBRAIN_WEIGHTS  Write calibration result to BabelBrain HDF5 format
%
% Packs PRESTUS per-element phases (and, when non-uniform, per-ring relative
% amplitude shading) into the OptimizedWeightsFile format expected by
% BabelBrain (key 'CALIBRATION', complex64 phasors, shape [N_elem x 1]).
%
% BabelBrain applies these weights as:
%   u0 *= OptimizedWeights          (in AdjustWeightAmplitudes)
%
% where u0 already carries geometric (Rayleigh-steered) phases AND a uniform
% SourceAmpPa amplitude per element (babel_integration_focused_annular_array.py:
% u0[ring] = SourceAmpPa*exp(1i*AllPhi[ring])). The weights must
% therefore be:
%   - a delta PHASE correction relative to those geometric phases:
%       angle(w[n]) = phases_rad[n] - geo_phases_rad[n]
%   - a relative AMPLITUDE ratio (elem_amp[n] / elem_amp[reference]) when
%     elem_amp specifies non-uniform per-ring shading (apodization) —
%     otherwise unit magnitude, since SourceAmpPa (from the reference ring)
%     already supplies the absolute amplitude.
%
% Provide geo_phases_rad (geometric phases at the calibration depth) so the
% correct delta is written. Without it the absolute phases are written, which
% doubles the geometric steering inside BabelBrain and produces wrong results.
%
% Use as:
%   export_babelbrain_weights(phases_rad, velocity, parameters, out_path, geo_phases_rad, elem_amp)
%
% Input:
%   phases_rad         - [1 x N_elem] optimised element phases [rad]
%   velocity           - calibrated scalar particle velocity [m/s]
%   parameters         - PRESTUS config (used for n_elem, freq, metadata)
%   out_path           - output .h5 file path (created or overwritten)
%   geo_phases_rad     - [1 x N_elem] geometric phases at the calibration
%                        depth [rad], from set_real_phases * pi/180.
%                        Required for correct BabelBrain interoperability.
%   elem_amp           - (optional) [1 x N_elem] per-ring drive amplitude
%                        (any consistent unit — only ratios matter), e.g.
%                        tr.annular.elem_p_amp. Non-uniform values are encoded
%                        as relative magnitude on the weights (ratio to
%                        elem_amp(1), matching the reference ring used
%                        elsewhere for --source_amp_pa/elem_p_amp(1)).
%                        Omit or leave uniform for unit-magnitude (phase-only)
%                        weights.
%
% Output:
%   HDF5 file at out_path in the layout BabelViscoFDTD's SaveToH5py writes
%   and ReadFromH5py (babel_integration_base.py) reads: every root item carries
%   a string attribute 'type', and complex arrays are an h5py-style compound
%   {r, i}.
%     /CALIBRATION       - [N_elem] complex64 (compound r/i float32) phasors
%                          encoding hardware-correction delta from geometric
%                          phases (and relative amplitude shading, if elem_amp
%                          given); type='ndarray'
%     /DRIVE_phases_rad, /DRIVE_amp_ratio - [N_elem x 1] intended absolute
%                          drive, used by run_babelbrain_case.py to rebuild
%                          the weights against BabelBrain's own AllPhi
%     /metadata/...      - provenance fields (freq, velocity, n_elem, timestamp)
%
%   N_elem is one weight per RING. Unpatched BabelBrain expects one weight per
%   sub-element (AdjustWeightAmplitudes); run_babelbrain_case.py expands the
%   per-ring weights. Files written before this layout (/CALIBRATION/real +
%   /imag groups) are still read by IMPORT_BABELBRAIN_WEIGHTS and converted by
%   run_babelbrain_case.py.
%
% See also: IMPORT_BABELBRAIN_WEIGHTS, PERFORM_GLOBAL_SEARCH,
%           BabelBrain/TranscranialModeling/babel_integration/integration_templates/babel_integration_base.py

arguments
    phases_rad      (1,:) {mustBeNumeric}
    velocity        (1,1) {mustBeNumeric}
    parameters      (1,1) struct
    out_path        (1,:) char
    geo_phases_rad  (1,:) {mustBeNumeric} = []
    elem_amp        (1,:) {mustBeNumeric} = []
end

n_elem = parameters.transducer(1).annular.elem_n;
assert(numel(phases_rad) == n_elem, ...
    'phases_rad length (%d) must match transducer.annular.elem_n (%d)', ...
    numel(phases_rad), n_elem);

% Compute delta phases relative to geometric steering.
% BabelBrain applies weights on top of u0 that already has geometric phases,
% so the weights must encode only the hardware correction.
if isempty(geo_phases_rad)
    is_layered = isfield(parameters, 'simulation') ...
              && isfield(parameters.simulation, 'medium') ...
              && contains(parameters.simulation.medium, 'layered');
    if is_layered
        error(['export_babelbrain_weights: geo_phases_rad is required for layered (subject) mode.\n' ...
            'BabelBrain applies weights on top of geometric phases it computes internally.\n' ...
            'Writing absolute phases would double the geometric steering and produce wrong results.\n' ...
            'Supply geo_phases_rad from tr.annular.geo_phase_rad or compute_phases().']);
    else
        warning(['export_babelbrain_weights: geo_phases_rad not provided. ' ...
            'Writing absolute phases — BabelBrain will double the geometric steering. ' ...
            'Pass geo_phases_rad to get correct interoperability.']);
    end
    delta_rad = phases_rad(:);
else
    assert(numel(geo_phases_rad) == n_elem, ...
        'geo_phases_rad length (%d) must match n_elem (%d)', numel(geo_phases_rad), n_elem);
    % angle(exp(i*(opt-geo))) wraps correctly to (-pi, pi]
    delta_rad = angle(exp(1i * (phases_rad(:) - geo_phases_rad(:))));
end

% Relative amplitude shading: ratio of each ring's drive amplitude to the
% reference ring (elem_amp(1), matching the ring babelbrain_wrapper.m uses
% for --source_amp_pa/elem_p_amp(1) -- BabelBrain's SourceAmpPa supplies that
% ring's absolute amplitude, so only the RATIO needs to be carried here).
% Uniform (or absent) elem_amp collapses to all-ones, i.e. unit-magnitude
% phase-only weights.
if isempty(elem_amp)
    amp_ratio = ones(n_elem, 1);
else
    assert(numel(elem_amp) == n_elem, ...
        'elem_amp length (%d) must match n_elem (%d)', numel(elem_amp), n_elem);
    assert(elem_amp(1) ~= 0, 'export_babelbrain_weights: elem_amp(1) (the reference ring) must be nonzero.');
    amp_ratio = elem_amp(:) / elem_amp(1);
    if any(abs(amp_ratio - 1) > 1e-9)
        fprintf('[export_babelbrain_weights] Non-uniform per-ring amplitude shading encoded (ratios: %s).\n', ...
            mat2str(round(amp_ratio', 4)));
    end
end

% /CALIBRATION weight = the phase delta on top of BabelBrain's geometric steer.
% BabelBrain drives each ring with angle(u0) = AllPhi + angle(weight), and
% geo_phases_rad is computed at the same aim as AllPhi, so
%   weight = exp(+i*delta_rad)  =>  AllPhi + delta_rad = phases_rad.
% This is only the fallback: run_babelbrain_case.py normally rebuilds the weight
% from /DRIVE (below) against BabelBrain's own AllPhi.
weights = amp_ratio .* complex(cos(delta_rad), sin(delta_rad));   % [N_elem x 1]

% Write HDF5 — overwrite if exists.
% Write to a local temp file first to avoid HDF5 file-locking errors on
% network/Docker filesystems where repeated h5create open/close cycles fail.
if isfile(out_path)
    delete(out_path);
end

write_path = [tempname, '.h5'];
tmp_guard  = onCleanup(@() h5_tmp_cleanup(write_path));

% RAW DRIVE (for the AllPhi-referenced correction path in run_babelbrain_case.py):
% the intended per-ring absolute drive phases (phases_rad) and relative
% amplitude shading (amp_ratio). When present, run_babelbrain_case.py rebuilds
% the weight as amp_ratio .* exp(1i*(phases_rad - AllPhi)) using BabelBrain's
% OWN internally-computed AllPhi, so the realised per-ring phase equals
% phases_rad EXACTLY (effective = AllPhi + [phases_rad - AllPhi] = phases_rad),
% independent of any analytic reconstruction of BabelBrain's geometric steer.
% The /CALIBRATION dataset below stays as a fallback (analytic delta) for
% consumers that do not apply this correction. Flat root names (not a /DRIVE
% group): ReadFromH5py requires a 'type' attribute on every root item.
h5create(write_path, '/DRIVE_phases_rad', [n_elem 1], 'Datatype', 'double');
h5write(write_path,  '/DRIVE_phases_rad', double(phases_rad(:)));
h5create(write_path, '/DRIVE_amp_ratio',  [n_elem 1], 'Datatype', 'double');
h5write(write_path,  '/DRIVE_amp_ratio',  double(amp_ratio(:)));

% Metadata for provenance
h5create(write_path, '/metadata/velocity_m_s',  [1 1], 'Datatype', 'double');
h5write(write_path,  '/metadata/velocity_m_s',  velocity);

h5create(write_path, '/metadata/freq_hz',        [1 1], 'Datatype', 'double');
h5write(write_path,  '/metadata/freq_hz',        parameters.transducer(1).freq_hz);

h5create(write_path, '/metadata/n_elem',         [1 1], 'Datatype', 'int32');
h5write(write_path,  '/metadata/n_elem',         int32(n_elem));

timestamp = int32(posixtime(datetime('now', 'TimeZone', 'UTC')));
h5create(write_path, '/metadata/created_utc',    [1 1], 'Datatype', 'int32');
h5write(write_path,  '/metadata/created_utc',    timestamp);

% CALIBRATION dataset and 'type' attributes need the low-level API: h5create
% cannot write a compound type. h5py maps complex64 to a compound with float32
% fields 'r' and 'i' (1-D, as numpy writes res.x in TxCalibration.py), and
% ReadFromH5py returns it as a complex ndarray.
fid = H5F.open(write_path, 'H5F_ACC_RDWR', 'H5P_DEFAULT');
fid_guard = onCleanup(@() H5F.close(fid));

cplx_type = H5T.create('H5T_COMPOUND', 8);
H5T.insert(cplx_type, 'r', 0, 'H5T_IEEE_F32LE');
H5T.insert(cplx_type, 'i', 4, 'H5T_IEEE_F32LE');
space = H5S.create_simple(1, n_elem, []);
dset  = H5D.create(fid, 'CALIBRATION', cplx_type, space, 'H5P_DEFAULT');
H5D.write(dset, cplx_type, 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', ...
    struct('r', single(real(weights)), 'i', single(imag(weights))));
H5D.close(dset); H5S.close(space); H5T.close(cplx_type);

% 'type' attributes as SaveToH5py sets them (ProcType): arrays 'ndarray',
% Python scalars 'scalar', dict groups 'dict' with key type 'str'.
h5_write_type_attr(fid, '/CALIBRATION', 'ndarray');
h5_write_type_attr(fid, '/DRIVE_phases_rad', 'ndarray');
h5_write_type_attr(fid, '/DRIVE_amp_ratio', 'ndarray');
h5_write_type_attr(fid, '/metadata', 'dict', 'str');
for name = {'velocity_m_s', 'freq_hz', 'n_elem', 'created_utc'}
    h5_write_type_attr(fid, ['/metadata/' name{1}], 'scalar');
end
clear fid_guard   % close before moving the file

movefile(write_path, out_path);

fprintf('BabelBrain weights saved: %s\n', out_path);
fprintf('  N_elem=%d  velocity=%.4f m/s\n', n_elem, velocity);
fprintf('\n  Python import:\n');
fprintf('    from BabelViscoFDTD.H5pySimple import ReadFromH5py\n');
fprintf('    w = ReadFromH5py(''%s'')[''CALIBRATION'']\n', out_path);

end


function h5_write_type_attr(fid, obj_path, type_str, type_key)
% Write SaveToH5py's string attribute 'type' (and 'type_key' for groups) as a
% variable-length UTF-8 scalar, which h5py returns as str.
obj_id = H5O.open(fid, obj_path, 'H5P_DEFAULT');
obj_guard = onCleanup(@() H5O.close(obj_id));
str_type = H5T.copy('H5T_C_S1');
H5T.set_size(str_type, 'H5T_VARIABLE');
H5T.set_cset(str_type, H5ML.get_constant_value('H5T_CSET_UTF8'));
space = H5S.create('H5S_SCALAR');
attrs = {'type', type_str};
if nargin > 3
    attrs(end+1, :) = {'type_key', type_key};
end
for k = 1:size(attrs, 1)
    attr_id = H5A.create(obj_id, attrs{k, 1}, str_type, space, 'H5P_DEFAULT');
    H5A.write(attr_id, 'H5ML_DEFAULT', attrs(k, 2));
    H5A.close(attr_id);
end
H5S.close(space); H5T.close(str_type);
end


function h5_tmp_cleanup(tmp_path)
if isfile(tmp_path), delete(tmp_path); end
end
