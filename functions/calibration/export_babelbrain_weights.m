function export_babelbrain_weights(phases_rad, velocity, parameters, out_path)
% EXPORT_BABELBRAIN_WEIGHTS  Write calibration result to BabelBrain HDF5 format
%
% Packs PRESTUS per-element phases and particle velocity into the
% OptimizedWeightsFile format expected by BabelBrain (key 'CALIBRATION',
% complex64 phasors, shape [N_elem x 1]).
%
% BabelBrain applies these weights as:
%   u0 *= OptimizedWeights          (in AdjustWeightAmplitudes)
%
% The phasors are constructed as:
%   w = amplitude_scale * exp(1i * phases_rad)
%
% where amplitude_scale = velocity / baseline_velocity normalises to
% BabelBrain's SourceAmpPa convention. If baseline_velocity is omitted, the
% phasors have unit amplitude (|w| = 1), i.e. phase-only calibration.
%
% Use as:
%   export_babelbrain_weights(phases_rad, velocity, parameters, out_path)
%   export_babelbrain_weights(phases_rad, velocity, parameters, out_path, baseline_velocity)
%
% Input:
%   phases_rad         - [1 x N_elem] optimised element phases [rad]
%   velocity           - calibrated scalar particle velocity [m/s]
%   parameters         - PRESTUS config (used for n_elem, freq, metadata)
%   out_path           - output .h5 file path (created or overwritten)
%
% Output:
%   HDF5 file at out_path with:
%     /CALIBRATION   - [N_elem x 1] complex64 phasors
%     /metadata/...  - provenance fields (freq, velocity, n_elem, timestamp)
%
% See also: PERFORM_GLOBAL_SEARCH, PERFORM_JOINT_DEPTH_FIT,
%           BabelBrain/TranscranialModeling/BabelIntegrationBASE.py

arguments
    phases_rad  (1,:) {mustBeNumeric}
    velocity    (1,1) {mustBeNumeric}
    parameters  (1,1) struct
    out_path    (1,:) char
end

n_elem = parameters.transducer.annular.elem_n;
assert(numel(phases_rad) == n_elem, ...
    'phases_rad length (%d) must match transducer.annular.elem_n (%d)', ...
    numel(phases_rad), n_elem);

% Build complex phasors (unit amplitude — phase-only)
weights = complex(cos(phases_rad(:)), sin(phases_rad(:)));   % [N_elem x 1]

% Write HDF5 — overwrite if exists
if exist(out_path, 'file')
    delete(out_path);
end

% CALIBRATION dataset: real and imaginary parts stored as separate float32
% arrays interleaved in a struct to match NumPy complex64 on read.
% MATLAB's h5create/h5write stores complex natively as compound type,
% which Python's h5py reads as structured array. Write real/imag separately
% and reconstruct in Python: weights = data['real'] + 1j*data['imag']
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
