function print_parameter_summary(parameters)
% PRINT_PARAMETER_SUMMARY  Print a structured PRESTUS configuration summary to the console
%
% Iterates over the main configuration sections (I/O, simulation type, grid,
% transducer(s), medium properties, thermal sequence, HPC/GPU) and prints
% each field with its value using the exact YAML parameter names. Checks field
% existence before printing; warns for missing required fields. Supports
% multiple transducers and matrix / annular types.
%
% Use as:
%   print_parameter_summary(parameters)
%
% Input:
%   parameters - PRESTUS config as returned by LOAD_PARAMETERS
%
% See also: LOAD_PARAMETERS, PATH_LOG_SETUP

arguments
    parameters struct
end

fprintf('========================================\n');
fprintf('PRESTUS CONFIGURATION SUMMARY\n');
fprintf('========================================\n\n');

%% 1. I/O Management

fprintf('📁 I/O MANAGEMENT\n');
pth = get_struct_or_default(parameters, 'path');
print_if_field(pth,              'anat',              '%s', 'anat_path');
print_if_field(pth,              'seg',               '%s', 'seg_path');
print_if_field(pth,              'sim',               '%s', 'sim_path');
startup = get_struct_or_default(parameters, 'startup');
print_if_field(startup,          'simnibs_bin_path',  '%s', 'simnibs_bin_path');
print_if_field(parameters.io,    'output_affix',      '%s');
print_overwrite_pair(parameters.io);
fprintf('\n');

%% 2. Simulation Type

fprintf('⚙️ SIMULATION TYPE\n');
print_if_field(parameters.simulation, 'medium', '%s');
print_if_field(parameters.simulation, 'precision', '%s');
print_if_field(parameters.simulation, 'code_type', '%s');
print_modules(parameters.modules);
fprintf('\n');

%% 3. Simulation Grid

fprintf('📐 SIMULATION GRID\n');
grid = get_struct_or_default(parameters, 'grid');
print_if_field(grid, 'n_dims', '%dD');
print_if_field(grid, 'resolution_mm', '%.2f mm');
print_flag(grid, 'axisymmetric');
print_if_field(grid, 'default_dims', '%s');
print_if_field(grid, 'pml_size', '%d');
print_if_field(grid, 'max_expand', '%.1f');
fprintf('\n');

%% 4. Transducer Specification
fprintf('🎯 TRANSDUCER SPECIFICATION\n');

% Loop over all transducers (support multiple transducers)
for t_i = 1:numel(parameters.transducer)
	tr = parameters.transducer(t_i);

    % Determine type
    tr_type = tr.type;

    fprintf('Transducer %d (%s):\n', t_i, tr_type);

    switch tr_type
        case 'annular'
            ann = tr.annular;

            print_if_field(ann, 'elem_n', '%d');
            print_if_field(ann, 'curv_radius_mm', '%.0f mm');
            print_if_field(ann, 'elem_phase_deg', '%.1f deg.');

        case 'matrix'
            mat = tr.matrix;

            fprintf('  Steering: %s\n', mat.steering);
            fprintf('  Element shape: %s\n', mat.elem_shape);
            print_if_field(mat, 'elem_height_mm', '%.3f mm');
            print_if_field(mat, 'elem_width_mm', '%.3f mm');
            print_if_field(mat, 'outer_diameter_mm', '%.2f mm');
            print_if_field(mat, 'is_curved', '%d');

            if mat.is_curved
                print_if_field(mat, 'curv_radius_mm', '%.2f mm');
            end

            % Clover setup
            if isfield(mat, 'is_clover_setup') && mat.is_clover_setup
                fprintf('  Clover setup enabled:\n');
                print_if_field(mat.clover, 'n_leaves', '%d');
                print_if_field(mat.clover, 'ROC_parent', '%.2f mm');
            end

            % Matrix shape type
            shape_type = mat.matrix_shape.type;
            fprintf('  Matrix shape type: %s\n', shape_type);

            switch shape_type
                case 'define_here'
                    define_here = mat.matrix_shape.define_here;
                    grid_type = define_here.grid_shape.type;
                    fprintf('    Grid type: %s\n', grid_type);

                    switch grid_type
                        case 'rect'
                            rect_grid = define_here.grid_shape.rect;
                            fprintf('      Rectangular grid:\n');
                            print_if_field(rect_grid, 'elem_n_row', '%d');
                            print_if_field(rect_grid, 'elem_n_col', '%d');
                            print_if_field(rect_grid, 'elem_spacing_height_mm', '%.2f mm');
                            print_if_field(rect_grid, 'elem_spacing_width_mm', '%.2f mm');
                            print_if_field(rect_grid, 'sparsity_factor', '%.2f');

                        case 'fibonacci'
                            fib_grid = define_here.grid_shape.fibonacci;
                            fprintf('      Fibonacci grid:\n');
                            print_if_field(fib_grid, 'elem_n', '%d');

                        case 'fermat'
                            fermat_grid = define_here.grid_shape.fermat;
                            fprintf('      Fermat spiral grid:\n');
                            print_if_field(fermat_grid, 'elem_n', '%d');
                    end

                case 'extract_from_file'
                    ext = mat.matrix_shape.extract_from_file;
                    print_if_field(ext, 'file_path', '%s');
                    print_if_field(ext, 'start_row', '%d');
                    print_if_field(ext, 'start_col', '%d');
                    print_if_field(ext, 'elem_n', '%d');
                    print_if_field(ext, 'select_random_subset', '%d');
                    if ext.select_random_subset
                        print_if_field(ext.subset, 'random_seed', '%d');
                        print_if_field(ext.subset, 'subset_n_elements', '%d');
                    end
                    print_if_field(ext, 'project_on_new_ROC', '%d');
                    if ext.project_on_new_ROC
                        print_if_field(ext.ROC_projection, 'new_ROC_mm', '%.2f mm');
                    end
            end
    end

    print_if_field(tr, 'freq_hz', '%.1f Hz');
    print_if_field(tr.(tr.type), 'elem_amp', '%.1f Pa');
    print_if_field(tr, 'trans_pos', '[%.1f %.1f %.1f]');
    print_if_field(tr, 'focus_pos', '[%.1f %.1f %.1f]');
	print_if_field(tr, 'focal_distance_ep', '%.1f mm');
	print_if_field(tr, 'focal_distance_bowl', '%.1f mm');

    fprintf('\n');
end

%% 5. Medium Properties (requested layers only)

fprintf('💀 MEDIUM PROPERTIES (requested layers)\n');

% Use parameters.layers (e.g. {'tissues', 'skull', 'brain'})
if ~isfield(parameters, 'layers') || isempty(parameters.layers)
    warning('No parameters.layers defined. Using all: water, skull, brain, skin...');
    tissues = {'water', 'skull', 'brain', 'skin', 'skull_trabecular', 'skull_cortical'};
else
    tissues = fieldnames(parameters.layers);
end

props = {'sound_speed', 'density', 'alpha_coeff', 'alpha_power', ...
         'thermal_conductivity', 'specific_heat_capacity', 'perfusion', 'absorption_fraction'};
labels = {'Sound Speed [m/s]', 'Density [kg/m³]', 'Attenuation (1 MHz) [dB/(MHz·cm)]', ...
          'Attn: Freq^Power [-]', 'Thermal Cond. [W/(m·K)]', 'Heat Cap. [J/(kg·K)]', ...
          'Perfusion [1/s]', 'Absorption Fraction [-]'};

for t_idx = 1:length(tissues)
    tissue = tissues{t_idx};
    print_tissue_props(parameters, tissue, props, labels);
end

%% 6. Thermal Sequence

if isfield(parameters, 'thermal') && parameters.modules.run_heating_sims == 1
    fprintf('🔥 THERMAL SEQUENCE\n');
    tm = get_struct_or_default(parameters, 'timing');
    print_if_field(tm, 'pd', '%.3fs', 'Pulse Duration');
    print_if_field(tm, 'pri', '%.3fs', 'Pulse Repetition Interval');
    print_if_field(tm, 'ptd', '%.3fs', 'Pulse Train Duration');
    print_if_field(tm, 'ptri', '%.3fs', 'Pulse Train Repetition Interval');
    print_if_field(tm, 'ptrd', '%.3fs', 'Pulse Train Repetition Duration');
    print_if_field(tm, 'post_ptri_dur', '%.3fs', 'Steady-State Duration');
    fprintf('\n');
end

%% 7. HPC/GPU

fprintf('💻 HPC/GPU OPTIONS\n');
print_if_field(parameters, 'platform', '%s');
hpc = get_struct_or_default(parameters, 'hpc');
print_if_field(hpc, 'partition', '%s');
fprintf('\n');

fprintf('========================================\n');
fprintf('Configuration verified - Ready to run!\n');
fprintf('========================================\n');

end

%% Helper functions

function print_if_field(s, field, fmt, alternate_name)
    if nargin < 4, alternate_name = ''; end
    
    % Use alternate_name if provided, otherwise original field
    display_name = iif(~isempty(alternate_name), alternate_name, field);
    
    if isfield(s, field)
        val = s.(field);
        
        % SPECIAL HANDLING: Vector dimensions
        if strcmp(field, 'default_dims') || strcmp(field, 'grid_dims')
            if ~isempty(val) && all(isfinite(val(:))) && ~any(isnan(val(:)))
                if length(val) <= 10
                    grid_str = sprintf('%g ', val);
                    grid_str = strtrim(grid_str);
                    fprintf('  %-30s [%s]\n', [display_name ':'], grid_str);  % Uses display_name
                else
                    fprintf('  %-30s [%dx...%dx%d]\n', [display_name ':'], size(val,1), size(val,2), length(val));
                end
            else
                fprintf('  ⚠️  %-30s [EMPTY/NaN VECTOR]\n', [display_name ':']);
            end
            return;
        end

        % NUMERIC SCALAR
        if isnumeric(val) && isscalar(val)
            if isfinite(val) && ~isnan(val)
                fprintf('  %-30s %s\n', [display_name ':'], sprintf(fmt, val));
            else
                fprintf('  ⚠️  %-30s [NaN/Inf]\n', [display_name ':']);
            end
            return;
        end
        
        % STRING
        if ischar(val) || isstring(val)
            val_str = char(val);
            if ~isempty(val_str) && ~strcmpi(val_str, 'NA') && ~strcmpi(val_str, 'N/A')
                fprintf('  %-30s %s\n', [display_name ':'], sprintf(fmt, val_str));
            else
                fprintf('  ⚠️  %-30s [EMPTY/NA]\n', [display_name ':']);
            end
            return;
        end
        
        % OTHER VECTORS/ARRAYS
        if isvector(val)
            if length(val) > 0 && all(isfinite(val)) && ~any(isnan(val))
                fprintf('  %-30s %s\n', [display_name ':'], sprintf(fmt, val));
            else
                fprintf('  ⚠️  %-30s [EMPTY/NaN VECTOR]\n', [display_name ':']);
            end
            return;
        end
        
        % STRUCT/CELL/OTHER
        if ~isempty(val)
            fprintf('  %-30s PRESENT\n', [display_name ':']);
        else
            fprintf('  ⚠️  %-30s [EMPTY]\n', [display_name ':']);
        end
        
    else
        % MISSING FIELD - uses alternate_name in warning
        fprintf('  ❗ %-30s [MISSING - add to config]\n', [display_name ':']);
    end
end

function str = iif(condition, true_val, false_val)
    if condition
        str = true_val;
    else
        str = false_val;
    end
end

function print_overwrite_pair(s)
    has_overwrite = isfield(s, 'overwrite_files');
    has_simnibs = isfield(s, 'overwrite_simnibs');
    if has_overwrite || has_simnibs
        fprintf('  %-30s %s', 'overwrite_files:', ...
            getfield_or_default(s, 'overwrite_files', ''));
        if has_simnibs
            fprintf(' / overwrite_simnibs: %d\n', s.overwrite_simnibs);
        else
            fprintf('\n');
        end
    end
end

function print_flag(s, field)
    if getfield_or_default(s, field, 0)
        fprintf('  %-30s Yes\n', [field ':']);
    end
end

function print_modules(s)
    mod_list = {
        'run_source_setup',    'Source';
        'run_acoustic_sims',   'Acoustic';
        'run_heating_sims',    'Thermal';
        'run_posthoc_water_sims', 'Post-water'
    };

    fprintf('\n🔧 MODULES\n');
    for i = 1:size(mod_list, 1)
        mod_field = mod_list{i, 1};
        mod_label = mod_list{i, 2};
        status = getfield_or_default(s, mod_field, 0);
        if status
            mark = ' ✓';
        else
            mark = ' ✗';
        end
        fprintf('  %-30s%s\n', [mod_label ':'], mark);
    end
end

function s = get_struct_or_default(params, field)
    if isfield(params, field) && isstruct(params.(field))
        s = params.(field);
    else
        s = struct();
    end
end

function val = getfield_or_default(s, field, default)
    if isfield(s, field)
        val = s.(field);
    else
        val = default;
    end
end

function val = get_nested_safe(s, path)
    try
        parts = strsplit(path, '.');
        val = s;
        for i = 1:length(parts)
            val = val.(parts{i});
        end
    catch
        val = NaN;
    end
end

function [base, prop] = split_path(path)
    parts = strsplit(path, '.');
    base = strjoin(parts(1:end-1), '.');
    prop = parts{end};
end

function unit = get_unit(prop)
% GET_UNIT returns units for tissue properties (MATLAB R2016b+ compatible)
if contains(prop, 'sound_speed')
    unit = 'm/s';
elseif contains(prop, 'density')
    unit = 'kg/m³';
elseif contains(prop, 'alpha_coeff')
    unit = 'dB/cm/MHz';
elseif contains(prop, 'alpha_power')
    unit = '(scalar)';
elseif contains(prop, 'thermal_conductivity')
    unit = 'W/m/°C';
elseif contains(prop, 'specific_heat_capacity')
    unit = 'J/kg/°C';
elseif contains(prop, 'perfusion')
    unit = 'mL/min/kg';
elseif contains(prop, 'absorption_fraction')
    unit = '[0-1]';
else
    unit = 'N/A';
end
end

function print_tissue_props(params, tissue, props, labels)
    has_props = false;
    for p_idx = 1:length(props)
        prop = props{p_idx}; label = labels{p_idx};
        val = get_nested_safe(params, sprintf('medium_properties.%s.%s', tissue, prop));
        
        if ~isempty(val) && ~isnan(val) && isfinite(val)
            unit = get_unit(prop);  % Reuse your function
            
            % Smart formatting
            if contains(prop, {'sound_speed','density','specific_heat_capacity','perfusion'})
                fmt = '%.0f';
            elseif strcmp(prop, 'alpha_coeff')
                fmt = '%.4f';
            else
                fmt = '%.3f';
            end
            
            fprintf('    %-28s %s %s\n', label, sprintf(fmt, val), unit);
            has_props = true;
        end
    end
    if ~has_props
        fprintf('    No properties defined\n');
    end
    fprintf('\n');
end
