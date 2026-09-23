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
% SourceAmpPa amplitude per element (BabelIntegrationANNULAR_ARRAY.py:
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
%   HDF5 file at out_path with:
%     /CALIBRATION/real, /CALIBRATION/imag - [N_elem x 1] float32 phasors
%                        encoding hardware-correction delta from geometric phases
%                        (and relative amplitude shading, if elem_amp given)
%     /DRIVE/phases_rad, /DRIVE/amp_ratio - [N_elem x 1] intended absolute drive,
%                        used by run_babelbrain_case.py to rebuild the weights
%                        against BabelBrain's own AllPhi
%     /metadata/...      - provenance fields (freq, velocity, n_elem, timestamp)
%
% See also: IMPORT_BABELBRAIN_WEIGHTS, PERFORM_GLOBAL_SEARCH,
%           BabelBrain/TranscranialModeling/BabelIntegrationBASE.py

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

% CALIBRATION dataset: real and imaginary parts stored as separate float32
% arrays to match NumPy complex64 on read.
% Reconstruct in Python: weights = data['real'] + 1j*data['imag']
h5create(write_path, '/CALIBRATION/real', [n_elem 1], 'Datatype', 'single');
h5create(write_path, '/CALIBRATION/imag', [n_elem 1], 'Datatype', 'single');
h5write(write_path,  '/CALIBRATION/real', single(real(weights)));
h5write(write_path,  '/CALIBRATION/imag', single(imag(weights)));

% RAW DRIVE (for the AllPhi-referenced correction path in run_babelbrain_case.py):
% the intended per-ring absolute drive phases (phases_rad) and relative
% amplitude shading (amp_ratio). When present, run_babelbrain_case.py rebuilds
% the weight as amp_ratio .* exp(1i*(phases_rad - AllPhi)) using BabelBrain's
% OWN internally-computed AllPhi, so the realised per-ring phase equals
% phases_rad EXACTLY (effective = AllPhi + [phases_rad - AllPhi] = phases_rad),
% independent of any analytic reconstruction of BabelBrain's geometric steer.
% The /CALIBRATION dataset above stays as a fallback (analytic delta) for
% consumers that do not apply this correction.
h5create(write_path, '/DRIVE/phases_rad', [n_elem 1], 'Datatype', 'double');
h5write(write_path,  '/DRIVE/phases_rad', double(phases_rad(:)));
h5create(write_path, '/DRIVE/amp_ratio',  [n_elem 1], 'Datatype', 'double');
h5write(write_path,  '/DRIVE/amp_ratio',  double(amp_ratio(:)));

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

movefile(write_path, out_path);

fprintf('BabelBrain weights saved: %s\n', out_path);
fprintf('  N_elem=%d  velocity=%.4f m/s\n', n_elem, velocity);
fprintf('\n  Python import:\n');
fprintf('    import h5py, numpy as np\n');
fprintf('    with h5py.File(''%s'', ''r'') as f:\n', out_path);
fprintf('        w = f[''/CALIBRATION/real''][:] + 1j*f[''/CALIBRATION/imag''][:]\n');

end


function h5_tmp_cleanup(tmp_path)
if isfile(tmp_path), delete(tmp_path); end
end
