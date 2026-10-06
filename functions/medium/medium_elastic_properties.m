function elastic_medium = medium_elastic_properties(parameters, kwave_medium, medium_masks)
% MEDIUM_ELASTIC_PROPERTIES  Build the k-Wave elastic medium (compressional + shear) from a fluid medium
%
% Converts the PRESTUS fluid medium (sound_speed, density, alpha_coeff,
% alpha_power) into the field set expected by pstdElastic2D/3D and adds a
% shear-wave sound speed and shear absorption inside the skull
% (skull, skull_cortical, skull_trabecular layers).
%
% Shear sound speed (elastic.shear_mapping), after Pichardo et al. (2017,
% Phys. Med. Biol. 62:6938), the viscoelastic model behind BabelBrain:
%   'pichardo'  cS = a(f0)*rho + b(f0) in skull voxels, with the published
%               fits a = 0.458, b = 545.52 at 270 kHz and a = 0.386,
%               b = 815.51 at 836 kHz, interpolated linearly in frequency
%               (extrapolated, with a warning, outside 270-836 kHz). Uses the
%               density map, so it applies to pseudo-CT and layer-based skulls
%               (layer densities) alike. (default)
%   'off'       cS = 0 everywhere (fluid skull through the elastic solver)
% Voxels with rho <= elastic.shear_density_min get no shear support. cS is
% capped so the Poisson ratio stays >= elastic.poisson_min:
%   cS <= cL * sqrt((1 - 2*nu_min) / (2*(1 - nu_min)))
%
% Absorption: pstdElastic uses a Kelvin-Voigt model with a fixed y = 2
% power law [dB/(MHz^2 cm)]. All coefficients are converted so that the
% attenuation at the source frequency f0 is exact: a*f0^y = a2*f0^2.
%   compressional - outside the skull (and inside, with
%                   elastic.long_attenuation 'medium'): the fluid
%                   alpha_coeff/alpha_power. In the skull with
%                   long_attenuation 'elastic' (default):
%                   elastic.alpha_coeff_long / alpha_power_long (default
%                   27.31 / 1.722: power law through the 33 / 231 Np/m at
%                   270 / 836 kHz that Pichardo et al. 2017 paired with their
%                   shear fit). The fluid skull values are effective values
%                   that already absorb mode-conversion losses.
%   shear         - elastic.alpha_coeff_shear / alpha_power_shear
%                   [dB/(MHz^y cm)] where cS > 0 (default 43.63 / 1.662:
%                   power law through Pichardo et al. 2017's 57 / 373 Np/m at
%                   270 / 836 kHz).
%
% Use as:
%   elastic_medium = medium_elastic_properties(parameters, kwave_medium, medium_masks)
%
% Input:
%   parameters   - PRESTUS config with elastic.*, transducer(1).freq_hz,
%                  medium_properties (layer labels), optionally
%                  grid.resolution_mm (shear points-per-wavelength check)
%   kwave_medium - fluid medium struct (sound_speed, density, alpha_coeff,
%                  alpha_power; other fields are ignored)
%   medium_masks - layer label map ([] or empty for a water-only medium)
%
% Output:
%   elastic_medium - struct with sound_speed_compression [m/s],
%                    sound_speed_shear [m/s], density [kg/m^3],
%                    alpha_coeff_compression and alpha_coeff_shear
%                    [dB/(MHz^2 cm)] (both omitted for a lossless medium)
%
% See also: ACOUSTIC_SIMULATION_ELASTIC, MEDIUM_SETUP, PSTDELASTIC3D

arguments
    parameters   (1,1) struct
    kwave_medium (1,1) struct
    medium_masks {mustBeNumericOrLogical} = []
end

    ec = get_elastic_config(parameters);
    f_MHz = parameters.transducer(1).freq_hz / 1e6;
    to_y2 = @(a, y) a .* f_MHz .^ (y - 2);   % exact at f0 under the y = 2 law

    density = kwave_medium.density;
    c_comp  = kwave_medium.sound_speed;

    %% Compressional absorption from the fluid medium
    if isfield(kwave_medium, 'alpha_coeff')
        alpha_power = 2;
        if isfield(kwave_medium, 'alpha_power'), alpha_power = kwave_medium.alpha_power; end
        alpha_comp = to_y2(kwave_medium.alpha_coeff, alpha_power);
    else
        alpha_comp = 0;   % lossless medium
    end

    %% Skull: shear speed, shear absorption, elastic longitudinal absorption
    c_shear     = zeros(size(density), 'like', density);
    alpha_shear = zeros(size(density), 'like', density);
    mapping     = 'off';
    if ~isempty(medium_masks)
        mask  = tissuemask_binary(parameters, medium_masks);
        skull = mask.skull;
        if ~any(skull(:))
            mapping = 'off';
        else
            if isscalar(density),   density = repmat(density, size(skull));   end
            if isscalar(c_comp),    c_comp  = repmat(c_comp,  size(skull));   end
            if isscalar(alpha_comp), alpha_comp = repmat(alpha_comp, size(skull)); end
            if ~isequal(size(skull), size(density))
                error('medium_elastic_properties:sizeMismatch', ...
                    'medium_masks size %s does not match the medium size %s.', ...
                    mat2str(size(skull)), mat2str(size(density)));
            end

            mapping = ec.shear_mapping;
            c_shear = zeros(size(skull), 'like', density);
            switch mapping
                case 'pichardo'
                    [a, b] = pichardo_shear_fit(parameters.transducer(1).freq_hz);
                    c_shear(skull) = a .* density(skull) + b;
                case 'off'
                otherwise
                    error('medium_elastic_properties:unknownMapping', ...
                        'Unknown elastic.shear_mapping "%s" (pichardo | off).', mapping);
            end
            c_shear(density <= ec.shear_density_min) = 0;

            % Poisson-ratio floor (bone: nu ~ 0.2-0.35)
            cs_max = c_comp .* sqrt((1 - 2*ec.poisson_min) / (2*(1 - ec.poisson_min)));
            capped = c_shear > cs_max;
            if any(capped(:))
                fprintf(['[medium_elastic] %d voxels (%.1f%% of shear voxels) capped to ' ...
                         'Poisson ratio >= %.2f\n'], nnz(capped), ...
                    100 * nnz(capped) / max(1, nnz(c_shear > 0)), ec.poisson_min);
                c_shear(capped) = cs_max(capped);
            end

            alpha_shear = zeros(size(skull), 'like', density);
            alpha_shear(c_shear > 0) = to_y2(ec.alpha_coeff_shear, ec.alpha_power_shear);

            switch ec.long_attenuation
                case 'elastic'
                    alpha_comp(skull) = to_y2(ec.alpha_coeff_long, ec.alpha_power_long);
                case 'medium'
                otherwise
                    error('medium_elastic_properties:unknownLongAttenuation', ...
                        'Unknown elastic.long_attenuation "%s" (elastic | medium).', ec.long_attenuation);
            end
        end
    end

    elastic_medium = struct( ...
        'sound_speed_compression', c_comp, ...
        'sound_speed_shear',       c_shear, ...
        'density',                 density, ...
        'alpha_coeff_compression', alpha_comp, ...
        'alpha_coeff_shear',       alpha_shear);
    % pstdElastic takes both absorption fields or neither; drop them for a
    % lossless medium to skip the Kelvin-Voigt terms.
    if ~any(alpha_comp(:)) && ~any(alpha_shear(:))
        elastic_medium = rmfield(elastic_medium, {'alpha_coeff_compression', 'alpha_coeff_shear'});
    end

    %% Report and resolution check
    if any(c_shear(:) > 0)
        cs_in = c_shear(c_shear > 0);
        db2np = 1 / (20 * log10(exp(1)));
        fprintf(['[medium_elastic] shear mapping %s: %d voxels, cS = %.0f–%.0f m/s (mean %.0f); ' ...
                 'at %.0f kHz: shear %.2f Np/cm, skull longitudinal %s\n'], ...
            mapping, numel(cs_in), min(cs_in), max(cs_in), mean(cs_in), f_MHz * 1e3, ...
            ec.alpha_coeff_shear * f_MHz^ec.alpha_power_shear * db2np, ...
            long_att_str(ec, f_MHz, db2np));
        if isfield(parameters, 'grid') && isfield(parameters.grid, 'resolution_mm')
            ppw = min(cs_in) / (parameters.transducer(1).freq_hz * parameters.grid.resolution_mm * 1e-3);
            if ppw < ec.min_ppw_shear
                warning('medium_elastic_properties:shearPPW', ...
                    ['Shear wave sampled at %.1f points per wavelength (cS = %.0f m/s, %.0f kHz, ' ...
                     '%.2f mm); elastic.min_ppw_shear is %g. Refine grid.resolution_mm to <= %.2f mm.'], ...
                    ppw, min(cs_in), f_MHz * 1e3, parameters.grid.resolution_mm, ec.min_ppw_shear, ...
                    min(cs_in) / (parameters.transducer(1).freq_hz * ec.min_ppw_shear) * 1e3);
            end
        end
    else
        fprintf('[medium_elastic] no shear support in medium (fluid run through elastic solver)\n');
    end
end

function [a, b] = pichardo_shear_fit(f0)
% Pichardo et al. 2017 (section 4.3): cS(rho) = a*rho + b, fitted at 270 and
% 836 kHz; linear in frequency in between.
    f_ref = [270e3 836e3];
    if f0 < f_ref(1) || f0 > f_ref(2)
        warning('medium_elastic_properties:shearSpeedExtrapolated', ...
            ['Source frequency %.0f kHz is outside the 270–836 kHz range of the ' ...
             'Pichardo et al. (2017) shear speed fits; extrapolating linearly.'], f0/1e3);
    end
    a = interp1(f_ref, [0.458 0.386],     f0, 'linear', 'extrap');
    b = interp1(f_ref, [545.52 815.51], f0, 'linear', 'extrap');
end

function s = long_att_str(ec, f_MHz, db2np)
    if strcmp(ec.long_attenuation, 'elastic')
        s = sprintf('%.2f Np/cm (elastic.alpha_coeff_long)', ec.alpha_coeff_long * f_MHz^ec.alpha_power_long * db2np);
    else
        s = 'from medium (fluid values)';
    end
end

function ec = get_elastic_config(parameters)
% Defaults mirror config_default.yaml (elastic section).
    ec = struct( ...
        'shear_mapping',     'pichardo', ...
        'shear_density_min', 1200, ...
        'poisson_min',       0.1, ...
        'alpha_coeff_shear', 43.63, ...
        'alpha_power_shear', 1.662, ...
        'alpha_coeff_long',  27.31, ...
        'alpha_power_long',  1.722, ...
        'long_attenuation',  'elastic', ...
        'min_ppw_shear',     6);
    if isfield(parameters, 'elastic') && isstruct(parameters.elastic)
        for f = fieldnames(ec)'
            if ~isfield(parameters.elastic, f{1}) || isempty(parameters.elastic.(f{1})), continue; end
            v = parameters.elastic.(f{1});
            if isstruct(ec.(f{1})) && isstruct(v)
                for g = fieldnames(v)'
                    ec.(f{1}).(g{1}) = v.(g{1});
                end
            else
                ec.(f{1}) = v;
            end
        end
    end
    ec.shear_mapping    = lower(char(ec.shear_mapping));
    ec.long_attenuation = lower(char(ec.long_attenuation));
end
