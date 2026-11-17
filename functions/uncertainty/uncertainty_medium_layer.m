function [layer, edits] = uncertainty_medium_layer(parameters, variant)
% UNCERTAINTY_MEDIUM_LAYER  Medium-property layer of one uncertainty variant
%
% The ITRUSST bounds of the variant (config/uncertainty/config_medium_<variant>.yaml)
% with the study's edits in simulation.uncertainty_overrides.<variant> merged
% on top. An edit that sets one absorption field of a tissue drops the other
% one from the bounds (as MEDIUM_ABSORPTION_LAYER does for config layers).
% Trans variants run acoustic only, so editing one of their thermal
% properties errors.
%
% Use as:
%   [layer, edits] = uncertainty_medium_layer(parameters, 'temp_high')
%
% Input:
%   parameters - PRESTUS config (reads simulation.uncertainty_overrides and,
%                to check the edited keys, medium_properties)
%   variant    - 'temp_high' | 'temp_low' | 'trans_high' | 'trans_low'
%
% Output:
%   layer - struct with a medium_properties field, merged onto the run's
%           config by UNCERTAINTY_PIPELINE
%   edits - the variant's edits (struct; no fields when the bounds are unchanged)
%
% See also: UNCERTAINTY_PIPELINE, UNCERTAINTY_VARIANTS, MEDIUM_ABSORPTION_LAYER

    variants = uncertainty_variants();
    variants = variants(~strcmp({variants.name}, 'default'));
    names = {variants.name};
    if ~any(strcmp(variant, names))
        error('uncertainty_medium_layer:badVariant', ...
            'Variant must be one of %s, not ''%s''.', strjoin(names, ', '), variant);
    end
    file = fullfile(get_prestus_path(), 'config', 'uncertainty', ['config_medium_' variant '.yaml']);
    if ~isfile(file)
        error('uncertainty_medium_layer:fileNotFound', 'Shipped bounds file not found: %s', file);
    end
    layer = yaml.loadFile(file, 'ConvertToArray', true);

    all_edits = struct();
    if isfield(parameters, 'simulation') && isfield(parameters.simulation, 'uncertainty_overrides') && ...
            isstruct(parameters.simulation.uncertainty_overrides)
        all_edits = parameters.simulation.uncertainty_overrides;
    end
    bad = setdiff(fieldnames(all_edits), names);
    if ~isempty(bad)
        error('uncertainty_medium_layer:badOverride', ...
            'simulation.uncertainty_overrides.%s: not a variant (use %s).', bad{1}, strjoin(names, ', '));
    end
    edits = struct();
    if isfield(all_edits, variant) && isstruct(all_edits.(variant))
        edits = all_edits.(variant);
    end
    if isempty(fieldnames(edits))
        return
    end

    heating = variants(strcmp(names, variant)).heating;
    thermal = {'thermal_conductivity', 'specific_heat_capacity', 'perfusion', ...
               'absorption_fraction', 'absorption_coeff'};
    absorption = {'absorption_fraction', 'absorption_coeff'};
    for t = fieldnames(edits)'
        tissue = t{1};
        key = sprintf('simulation.uncertainty_overrides.%s.%s', variant, tissue);
        if ~isfield(parameters, 'medium_properties') || ~isfield(parameters.medium_properties, tissue)
            error('uncertainty_medium_layer:badOverride', '%s: not a tissue of medium_properties.', key);
        end
        props = edits.(tissue);
        if ~isstruct(props)
            error('uncertainty_medium_layer:badOverride', '%s must list properties.', key);
        end
        for q = fieldnames(props)'
            prop = q{1};
            value = props.(prop);
            if ~isfield(parameters.medium_properties.(tissue), prop)
                error('uncertainty_medium_layer:badOverride', '%s.%s: not a medium property.', key, prop);
            end
            if ~(isnumeric(value) && isscalar(value) && isfinite(value))
                error('uncertainty_medium_layer:badOverride', '%s.%s must be a number.', key, prop);
            end
            if ~heating && any(strcmp(prop, thermal))
                error('uncertainty_medium_layer:badOverride', ...
                    '%s.%s: this variant runs acoustic only, so thermal properties are not used.', key, prop);
            end
        end
        set_abs = isfield(props, absorption);
        if all(set_abs)
            error('uncertainty_medium_layer:badOverride', ...
                '%s sets both absorption_fraction and absorption_coeff; set one of them.', key);
        end
        if any(set_abs) && isfield(layer.medium_properties, tissue)
            other = absorption{~set_abs};
            if isfield(layer.medium_properties.(tissue), other)
                layer.medium_properties.(tissue).(other) = [];
            end
        end
    end
    layer.medium_properties = mergestruct(layer.medium_properties, edits);
end
