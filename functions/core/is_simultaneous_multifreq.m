function result = is_simultaneous_multifreq(parameters)
% IS_SIMULTANEOUS_MULTIFREQ  True for coherent transducers firing together at different frequencies
%
% Transducers that fire at the same time at more than one carrier frequency
% (simulation.transducer_coupling 'coherent') are simulated one at a time and
% combined by COMBINE_ASYNC_INTENSITY in 'simultaneous' mode: acoustics is
% linear, so the joint field is the sum of the complex fields. Transducers at
% one frequency are summed coherently (they interfere); the time-averaged
% intensity and heat source are sums over the frequencies, whose cross terms
% average out. ASYNC_TRANSDUCER_PIPELINE runs it (IS_ASYNC_MODE is true for it).
%
% Use as:
%   tf = is_simultaneous_multifreq(parameters)
%
% See also: IS_ASYNC_MODE, ASYNC_TRANSDUCER_PIPELINE, COMBINE_ASYNC_INTENSITY

result = false;
if ~isfield(parameters, 'transducer') || numel(parameters.transducer) < 2
    return
end
% stages of ASYNC_TRANSDUCER_PIPELINE run as ordinary jobs
if isfield(parameters, 'state') && isfield(parameters.state, 'async_subjob') && parameters.state.async_subjob
    return
end
% benchmark only: one joint k-Wave run with every tone, as a reference for the
% decomposition (benchmark script outside the toolbox). Its p_max_all
% is the peak over the recorded window, and its intensity, MI and heat source
% assume one frequency.
if isfield(parameters, 'state') && isfield(parameters.state, 'multifreq_joint_reference') && ...
        parameters.state.multifreq_joint_reference
    return
end
coupling = 'coherent';
if isfield(parameters, 'simulation') && isfield(parameters.simulation, 'transducer_coupling') && ...
        ~isempty(parameters.simulation.transducer_coupling)
    coupling = char(parameters.simulation.transducer_coupling);
end
freqs = [parameters.transducer.freq_hz];
result = strcmp(coupling, 'coherent') && numel(unique(freqs)) > 1;
end
