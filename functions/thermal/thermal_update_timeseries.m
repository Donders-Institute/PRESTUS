function timeseries = thermal_update_timeseries(parameters, medium_masks, timeseries, curT, curCEM43, curCEM43_iso)
    % THERMAL_UPDATE_TIMESERIES Track max T/CEM43/CEM43_iso per tissue layer
    %
    % timeseries = thermal_update_timeseries(params, medium_masks)            — init
    % timeseries = thermal_update_timeseries(params, medium_masks, ts, T, CEM43)
    % timeseries = thermal_update_timeseries(params, medium_masks, ts, T, CEM43, CEM43_iso)

    % Get masks
    mask = tissuemask_binary(parameters, medium_masks);
    layers = fieldnames(mask);

    track_iso = nargin >= 6 && ~isempty(curCEM43_iso);

    % Initialize timeseries structure if missing
    if nargin < 3 || isempty(timeseries)
        timeseries.T      = struct();
        timeseries.Tdiff  = struct();
        timeseries.CEM43     = struct();
        timeseries.CEM43_iso = struct();
        for i_layer = 1:numel(layers)
            layer_mask = mask.(layers{i_layer});
            if nnz(layer_mask) > 0
                timeseries.T.(layers{i_layer})        = [];
                timeseries.Tdiff.(layers{i_layer})    = [];
                timeseries.CEM43.(layers{i_layer})     = [];
                timeseries.CEM43_iso.(layers{i_layer}) = [];
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
                % CEM43 max (kWave)
                timeseries.CEM43.(layers{i_layer})(end+1) = max(curCEM43(layer_mask));
                % CEM43 max (ISO)
                if track_iso
                    timeseries.CEM43_iso.(layers{i_layer})(end+1) = max(curCEM43_iso(layer_mask));
                end
            end
        end
    end
end

