function s = library_static_model(lib)
% LIBRARY_STATIC_MODEL  Static model of a library entry, in FIT_TRANSDUCER_MODEL form
%
% Use as:
%   prior = library_static_model(lib)      % lib from yaml.loadFile
%
% Output:
%   s - struct with any of elem_phase_correction_deg, elem_amp_weight,
%       steering_sound_speed, z_shift_mm (fields absent in the entry are
%       omitted, i.e. geometric / unit gain / default speed / no shift)
%
% See also: FIT_TRANSDUCER_MODEL, REFINE_TRANSDUCER_LIBRARY

    s = struct();
    if ~isfield(lib, 'global_model') || ~isstruct(lib.global_model); return; end
    gm  = lib.global_model;
    map = {'elem_phase_correction_deg', 'elem_phase_correction_deg'; ...
           'elem_amp_weight',           'elem_amp_weight'; ...
           'steering_sound_speed',      'steering_sound_speed'; ...
           'fit_z_shift_mm',            'z_shift_mm'};
    for k = 1:size(map, 1)
        if isfield(gm, map{k, 1}) && ~isempty(gm.(map{k, 1})); s.(map{k, 2}) = double(gm.(map{k, 1})(:)'); end
    end
end
