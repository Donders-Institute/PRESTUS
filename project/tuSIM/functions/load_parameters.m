function parameters = load_parameters(varargin)
    parameters = yaml.loadFile('configs/default_config.yaml', "ConvertToArray", true);
    
    if nargin == 1
        extra_config_file = varargin{1};
        extra_parameters = yaml.loadFile(fullfile('configs',extra_config_file), "ConvertToArray", true);
        parameters = mergeStructure(parameters, extra_parameters);
    end
    if ~isfield(parameters.transducer,'source_phase_rad')
        assert(isfield(parameters.transducer,'source_phase_deg'), 'Source phase should be set in transducer parameters as source_phase_rad or source_phase_deg')
        parameters.transducer.source_phase_rad = parameters.transducer.source_phase_deg/180*pi;
    end
    if length(parameters.transducer.source_amp)==1 && parameters.transducer.n_elements > 1
        parameters.transducer.source_amp = repmat(parameters.transducer.source_amp, [1 parameters.transducer.n_elements]);
    end
    parameters.grid_step_m = parameters.grid_step_mm/1e3; % [meters]
    parameters.default_grid_dims = repmat(parameters.default_grid_size, [1 3]);
    if iscell(parameters.transducer.source_phase_rad)
        for i = 1:length(parameters.transducer.source_phase_rad)
            if ~isnumeric(parameters.transducer.source_phase_rad{i})
                parameters.transducer.source_phase_rad{i}= eval(parameters.transducer.source_phase_rad{i});
            end
        end
        parameters.transducer.source_phase_rad = cell2mat(parameters.transducer.source_phase_rad);
    end
end