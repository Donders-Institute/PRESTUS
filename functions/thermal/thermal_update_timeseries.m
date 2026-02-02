function timeseries = thermal_update_timeseries(parameters, medium_masks, timeseries, curT, curCEM43)
    % THERMAL_UPDATE_TIMESERIES Track max T/CEM43 per tissue layer
    % 
    % timeseries = thermal_update_timeseries(params, medium_masks, curT, curCEM43)
    % timeseries = thermal_update_timeseries(params, medium_masks, curT, curCEM43, timeseries)
    
    % Get masks
    mask = tissuemask_binary(parameters, medium_masks);
    layers = fieldnames(mask);

    % Initialize timeseries structure if missing
    if nargin < 3 || isempty(timeseries)
        timeseries.T = struct();
        timeseries.CEM43 = struct();
        for i_layer = 1:numel(layers)
            layer_mask = mask.(layers{i_layer});
            if nnz(layer_mask) > 0
                timeseries.T.(layers{i_layer}) = [];
                timeseries.Tdiff.(layers{i_layer}) = [];
                timeseries.CEM43.(layers{i_layer}) = [];
            end
        end
    else    
        % Loop layers
        for i_layer = 1:numel(layers)
            layer_mask = mask.(layers{i_layer});
            if nnz(layer_mask) > 0
                % Temperature max
                timeseries.T.(layers{i_layer})(end+1) = max(curT(layer_mask));
                % Temperature change
                timeseries.Tdiff.(layers{i_layer})(end+1) = max(curT(layer_mask))-...
                    parameters.thermal.temp_0.(layers{i_layer});
                % CEM43 max  
                timeseries.CEM43.(layers{i_layer})(end+1) = max(curCEM43(layer_mask));
            end
        end
    end
end

