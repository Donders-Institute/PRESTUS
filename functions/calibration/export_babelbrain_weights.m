function export_babelbrain_weights(phases_rad, velocity, parameters, out_path, geo_phases_rad)
% EXPORT_BABELBRAIN_WEIGHTS  Write calibration result to BabelBrain HDF5 format
%
% Packs PRESTUS per-element phases and particle velocity into the
% OptimizedWeightsFile format expected by BabelBrain (key 'CALIBRATION',
% complex64 phasors, shape [N_elem x 1]).
%
% BabelBrain applies these weights as:
%   u0 *= OptimizedWeights          (in AdjustWeightAmplitudes)
%
% where u0 already carries geometric (Rayleigh-steered) phases per element
% (BabelIntegrationANNULAR_ARRAY.py line 399). The weights must therefore
% be delta corrections relative to those geometric phases:
%   w[n] = exp(1i * (phases_rad[n] - geo_phases_rad[n]))
%
% Provide geo_phases_rad (geometric phases at the calibration depth) so the
% correct delta is written. Without it the absolute phases are written, which
% doubles the geometric steering inside BabelBrain and produces wrong results.
%
% Use as:
%   export_babelbrain_weights(phases_rad, velocity, parameters, out_path, geo_phases_rad)
%
% Input:
%   phases_rad         - [1 x N_elem] optimised element phases [rad]
%   velocity           - calibrated scalar particle velocity [m/s]
%   parameters         - PRESTUS config (used for n_elem, freq, metadata)
%   out_path           - output .h5 file path (created or overwritten)
%   geo_phases_rad     - [1 x N_elem] geometric phases at the calibration
%                        depth [rad], from set_real_phases * pi/180.
%                        Required for correct BabelBrain interoperability.
%
% Output:
%   HDF5 file at out_path with:
%     /CALIBRATION/real, /CALIBRATION/imag - [N_elem x 1] float32 phasors
%                        encoding hardware-correction delta from geometric phases
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
end

n_elem = parameters.transducer.annular.elem_n;
assert(numel(phases_rad) == n_elem, ...
    'phases_rad length (%d) must match transducer.annular.elem_n (%d)', ...
    numel(phases_rad), n_elem);

% Compute delta phases relative to geometric steering.
% BabelBrain applies weights on top of u0 that already has geometric phases,
% so the weights must encode only the hardware correction.
if isempty(geo_phases_rad)
    warning(['export_babelbrain_weights: geo_phases_rad not provided. ' ...
        'Writing absolute phases — BabelBrain will double the geometric steering. ' ...
        'Pass geo_phases_rad to get correct interoperability.']);
    delta_rad = phases_rad(:);
else
    assert(numel(geo_phases_rad) == n_elem, ...
        'geo_phases_rad length (%d) must match n_elem (%d)', numel(geo_phases_rad), n_elem);
    % angle(exp(i*(opt-geo))) wraps correctly to (-pi, pi]
    delta_rad = angle(exp(1i * (phases_rad(:) - geo_phases_rad(:))));
end

weights = complex(cos(delta_rad), sin(delta_rad));   % [N_elem x 1], unit amplitude

% Write HDF5 — overwrite if exists
if exist(out_path, 'file')
    delete(out_path);
end

% CALIBRATION dataset: real and imaginary parts stored as separate float32
% arrays to match NumPy complex64 on read.
% Reconstruct in Python: weights = data['real'] + 1j*data['imag']
h5create(out_path, '/CALIBRATION/real', [n_elem 1], 'Datatype', 'single');
h5create(out_path, '/CALIBRATION/imag', [n_elem 1], 'Datatype', 'single');
h5write(out_path,  '/CALIBRATION/real', single(real(weights)));
h5write(out_path,  '/CALIBRATION/imag', single(imag(weights)));

% Metadata for provenance
h5create(out_path, '/metadata/velocity_m_s',  [1 1], 'Datatype', 'double');
h5write(out_path,  '/metadata/velocity_m_s',  velocity);

h5create(out_path, '/metadata/freq_hz',        [1 1], 'Datatype', 'double');
h5write(out_path,  '/metadata/freq_hz',        parameters.transducer.freq_hz);

h5create(out_path, '/metadata/n_elem',         [1 1], 'Datatype', 'int32');
h5write(out_path,  '/metadata/n_elem',         int32(n_elem));

timestamp = int32(posixtime(datetime('now', 'TimeZone', 'UTC')));
h5create(out_path, '/metadata/created_utc',    [1 1], 'Datatype', 'int32');
h5write(out_path,  '/metadata/created_utc',    timestamp);

fprintf('BabelBrain weights saved: %s\n', out_path);
fprintf('  N_elem=%d  velocity=%.4f m/s\n', n_elem, velocity);
fprintf('\n  Python import:\n');
fprintf('    import h5py, numpy as np\n');
fprintf('    with h5py.File(''%s'', ''r'') as f:\n', out_path);
fprintf('        w = f[''/CALIBRATION/real''][:] + 1j*f[''/CALIBRATION/imag''][:]\n');

end
