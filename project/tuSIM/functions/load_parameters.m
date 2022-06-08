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
    if ~isfield(parameters.transducer,'dist_to_plane_mm')
        parameters.transducer.dist_to_plane_mm = sqrt(parameters.transducer.curv_radius_mm^2-(max(parameters.transducer.Elements_OD_mm)/2)^2);
        fprintf('Distance to transducer plane is not provided, calculated based on curvature radius and parameters.transducer diameter as %.2f mm\n', parameters.transducer.dist_to_plane_mm)
    end

    if length(parameters.transducer.source_amp)==1 && parameters.transducer.n_elements > 1
        parameters.transducer.source_amp = repmat(parameters.transducer.source_amp, [1 parameters.transducer.n_elements]);
    end
    
    parameters.grid_step_m = parameters.grid_step_mm/1e3; % [meters]
    if ~isfield(parameters, 'n_sim_dims')
        if isfield(parameters,'default_grid_dims')
            parameters.n_sim_dims = length(parameters.default_grid_dims);
        else
            parameters.n_sim_dims = 3;
        end            
    end
    
    if ~isfield(parameters,'default_grid_dims')
        parameters.default_grid_dims = repmat(parameters.default_grid_size, [1 parameters.n_sim_dims]);
    end
    if iscell(parameters.transducer.source_phase_rad)
        for i = 1:length(parameters.transducer.source_phase_rad)
            if ~isnumeric(parameters.transducer.source_phase_rad{i})
                parameters.transducer.source_phase_rad{i}= eval(parameters.transducer.source_phase_rad{i});
            end
        end
        parameters.transducer.source_phase_rad = cell2mat(parameters.transducer.source_phase_rad);
    end
end