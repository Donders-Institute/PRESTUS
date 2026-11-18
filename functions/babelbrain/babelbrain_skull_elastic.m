function skull = babelbrain_skull_elastic(parameters)
% BABELBRAIN_SKULL_ELASTIC  Skull elastic properties for BabelBrain, matched to kwave_elastic
%
% Evaluates the skull properties that acoustic_backend 'kwave_elastic' uses
% (MEDIUM_ELASTIC_PROPERTIES with the elastic: section of the config) for each
% skull layer in medium_properties (skull, skull_cortical, skull_trabecular),
% at the drive frequency, in the units of BabelBrain's material table
% (babelbrain.shear: 'prestus_elastic').
%
% BabelBrain's non-CT path has one material per label, so the properties are
% computed from the layer density and sound speed. This matches kwave_elastic
% on a label-based skull; with a pseudo-CT or smoothed property maps,
% kwave_elastic varies cS voxel-wise and BabelBrain does not.
%
% Units: medium_elastic_properties returns attenuation as a y = 2 coefficient
% a2 [dB/(cm MHz^2)] that is exact at the drive frequency f0
% (a2 * f0^2 = a * f0^y). BabelBrain's AddMaterial takes attenuation in Np/m at
% f0, so
%   att [Np/m] = a2 * (f0 / 1 MHz)^2 * 100 / (20 log10(e))
% (DB_CM_MHZ_TO_NP_M with power 2; 100 cm per m, 8.686 dB per neper).
%
% Use as:
%   skull = babelbrain_skull_elastic(parameters)
%
% Input:
%   parameters - PRESTUS config with medium_properties, transducer(1).freq_hz
%                and (optionally) elastic.*
%
% Output:
%   skull - struct with one field per skull layer present in medium_properties
%           (skull, skull_cortical, skull_trabecular), each a struct with
%           density [kg/m^3], sound_speed [m/s], shear_sound_speed [m/s],
%           long_att_np_m and shear_att_np_m [Np/m at f0]
%
% See also: MEDIUM_ELASTIC_PROPERTIES, BABELBRAIN_WRAPPER, DB_CM_MHZ_TO_NP_M

arguments
    parameters (1,1) struct
end

    mp     = parameters.medium_properties;
    labels = fieldnames(mp);
    n      = numel(labels);
    freq   = parameters.transducer(1).freq_hz;

    % One "voxel" per medium_properties label: label index i has the layer's
    % properties, so tissuemask_binary finds the skull layers by name.
    medium = struct('sound_speed', nan(1, n), 'density', nan(1, n), ...
                    'alpha_coeff', zeros(1, n), 'alpha_power', ones(1, n));
    for i = 1:n
        m = mp.(labels{i});
        if ~isstruct(m), continue; end
        if isfield(m, 'sound_speed') && ~isempty(m.sound_speed), medium.sound_speed(i) = m.sound_speed; end
        if isfield(m, 'density')     && ~isempty(m.density),     medium.density(i)     = m.density;     end
        if isfield(m, 'alpha_coeff') && ~isempty(m.alpha_coeff), medium.alpha_coeff(i) = m.alpha_coeff; end
        if isfield(m, 'alpha_power') && ~isempty(m.alpha_power), medium.alpha_power(i) = m.alpha_power; end
    end

    % No grid field: the shear points-per-wavelength check refers to the
    % k-Wave grid; BabelBrain sizes its own grid from the slowest wave.
    p = struct('medium_properties', mp, 'transducer', parameters.transducer(1));
    if isfield(parameters, 'elastic'), p.elastic = parameters.elastic; end
    em = medium_elastic_properties(p, medium, 1:n);
    if ~isfield(em, 'alpha_coeff_compression')
        em.alpha_coeff_compression = zeros(1, n);
        em.alpha_coeff_shear       = zeros(1, n);
    end

    skull = struct();
    for name = {'skull', 'skull_cortical', 'skull_trabecular'}
        i = find(strcmp(labels, name{1}), 1);
        if isempty(i), continue; end
        if ~isfinite(medium.density(i)) || ~isfinite(medium.sound_speed(i))
            error('babelbrain_skull_elastic:missingProperty', ...
                'medium_properties.%s needs density and sound_speed for babelbrain.shear ''prestus_elastic''.', name{1});
        end
        skull.(name{1}) = struct( ...
            'density',           medium.density(i), ...
            'sound_speed',       medium.sound_speed(i), ...
            'shear_sound_speed', em.sound_speed_shear(i), ...
            'long_att_np_m',     db_cm_mhz_to_np_m(em.alpha_coeff_compression(i), freq, 2), ...
            'shear_att_np_m',    db_cm_mhz_to_np_m(em.alpha_coeff_shear(i), freq, 2));
    end
    if isempty(fieldnames(skull))
        error('babelbrain_skull_elastic:noSkull', ...
            'babelbrain.shear ''prestus_elastic'' needs a skull, skull_cortical or skull_trabecular layer in medium_properties.');
    end
end
