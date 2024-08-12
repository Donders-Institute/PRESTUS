close all; clear;

load("C:\Users\marge\OneDrive - Radboud Universiteit\Documenten\GitHub\currentPRESTUS\new_3D_layer_check\visualize_layers.mat")

segmentation_check(parameters, t1_img_rr, segmented_img_rr, 1)

function segmentation_check(parameters, t1_img_rr, segmented_img_rr, down_sampling_factor)
    % prepare mask data
    disp("Prepare mask data...")
    label_names = ["skull_cortical" "skull_trabecular" "skin" "brain" "water"];
    colors = {[1 0 0], [1 1 0], [0 1 1], [0 1 0], [0 0 1]};  % RGB colors as cell array
    masks = cell(length(label_names), 1);
    for i  = 1:length(label_names)
        l_nam = label_names(i);
        idx = [];
        labels = fieldnames(parameters.layer_labels);
        idx_label = find(contains(labels, l_nam));
    
        if ~isempty(idx_label)
            for numidx = 1:numel(idx_label)
                idx = [idx, parameters.layer_labels.(labels{idx_label(numidx)})];
            end
        end
    
        % unsmoothed masks
        mask_unsmoothed = ismember(segmented_img_rr, idx);

        % Downsampling by a factor of n to decrease computation time, a
        % factor of 1 means no downsampling
        mask_unsmoothed = mask_unsmoothed(1:down_sampling_factor:end, 1:down_sampling_factor:end, 1:down_sampling_factor:end);
    
        masks{i} = mask_unsmoothed;
    
    end
    
    disp("Plotting...")
    % Downsampling by a factor of n to decrease computation time, a
    % factor of 1 means no downsampling
    t1_img_rr = t1_img_rr(1:down_sampling_factor:end, 1:down_sampling_factor:end, 1:down_sampling_factor:end);

    % Normalize T1 image to range [0, 1]
    t1_img_rr = mat2gray(t1_img_rr);
    
    % Initialize the 4D volume for sliceViewer
    combined_volume = zeros(size(t1_img_rr, 1), size(t1_img_rr, 2), size(t1_img_rr, 3), 3);
    
    % Add the T1 image to the combined volume as grayscale
    combined_volume(:,:,:,1) = t1_img_rr;
    combined_volume(:,:,:,2) = t1_img_rr;
    combined_volume(:,:,:,3) = t1_img_rr;

    default_volume = combined_volume;
    
    % Initialize figure
    fig = uifigure("HandleVisibility","on", "Name","Mask checker");
    fig.Position(3:4) = [580, 480];
    ViewPnl = uipanel(fig, "Title", "Slices","Position", [160, 20, 400, 400]);
    
    % Checkbox callback function
    function updateVolume(~, ~)
        % Remove old figure
        delete(ViewPnl);
        ViewPnl = uipanel(fig, "Title", "Slices","Position", [160, 20, 400, 400]);

        % Define alpha for blending
        alpha = ef.Value;

        % Reset combined volume
        combined_volume = default_volume;
        for i_lb = 1:length(label_names)
            if checkboxes(i_lb).Value == 1
                for j = 1:3
                    combined_volume(:,:,:,j) = combined_volume(:,:,:,j) .* (1 - alpha * masks{i_lb}) + ...
                        colors{i_lb}(j) * masks{i_lb} * alpha;
                end
            end
        end
        combined_volume = min(combined_volume, 1);
        sliceViewer(combined_volume, "Parent", ViewPnl);
        updateAxis();
    end
    
    
    % Create checkboxes
    checkboxes = gobjects(length(label_names), 1);
    for i_ch = 1:length(label_names)
        checkboxes(i_ch) = uicheckbox(fig, "Text", label_names{i_ch},...
            "Value", 1, "ValueChangedFcn", @updateVolume, "Position",...
            [20, 400 - (i_ch-1)*30, 120, 20], "FontColor", colors{i_ch}*0.6,...
            "FontWeight", "bold");
    end

    % Dropdown callback function
    function updateAxis(~, ~)
        % Remove old figure
        delete(ViewPnl);
        ViewPnl = uipanel(fig, "Title", "Slices","Position", [160, 20, 400, 400]);
        
        plane = dd.Value;
        s = sliceViewer(combined_volume, "Parent", ViewPnl);
        if plane == "Transverse plane"
            s.SliceDirection = [1 0 0];
        elseif plane == "Midsagittal plane"
            s.SliceDirection = [0 0 1];
        elseif plane == "Coronal plane"
            s.SliceDirection = [0 1 0];
        end
    end

    % Create axis dropdown
    dd = uidropdown(fig, "Items", ["Transverse plane", "Midsagittal plane", "Coronal plane"],...
        "Position", [160, 440, 200, 20], "ValueChangedFcn", @updateAxis);

    % Create alpha field
    uilabel(fig, "Text", "Opaqueness", "Position",[375, 440, 100, 20]);
    ef = uieditfield(fig, "numeric", "Value", 0.2, "Limits",[0 1],...
        "Position",[460, 440, 100, 20], "ValueChangedFcn", @updateVolume);

    % Initial display
    updateVolume();
    
end