function [mask, info] = air_cavity_mask(parameters, segmentation, pseudoCT, charm_air)
% AIR_CAVITY_MASK  Internal air cavities (sinuses, mastoid air cells) on the simulation grid
%
% Builds the mask that ACOUSTIC_ANALYSIS uses for the air-cavity beam-path
% check (headmodel.air_cavities). The mask is a read-out only: the medium is
% built from parameters.layers as usual, so cavity voxels keep the
% properties of their layer (water unless layers assigns them).
%
%   source 'simnibs'      - the SimNIBS charm 'Air-Internal' mask (AIR_CAVITY_CHARM),
%                           passed in as charm_air
%   source 'segmentation' - voxels whose segmentation label is in
%                           headmodel.air_cavities.labels
%   source 'hu'           - CT/pseudo-CT voxels in headmodel.air_cavities.hu_window,
%                           as BabelBrain's air regions: majority (median)
%                           filter over the 3x3(x3) neighbourhood, connected
%                           components (26-connectivity), the largest
%                           component (air around the head) removed and,
%                           unless keep_boundary, every component touching
%                           the grid edge
%   source 'none'         - no mask (returns [])
%
% Use as:
%   [mask, info] = air_cavity_mask(parameters, segmentation, pseudoCT, charm_air)
%
% Input:
%   parameters   - PRESTUS config (headmodel.air_cavities)
%   segmentation - segmentation label volume on the simulation grid
%   pseudoCT     - CT/pseudo-CT HU volume on the same grid, or [] (source 'hu' needs it)
%   charm_air    - logical charm air mask on the same grid, or [] (source 'simnibs' needs it)
%
% Output:
%   mask - logical volume the size of segmentation, or [] when off or unavailable
%   info - struct: source, n_vox, n_components, note (why the mask is empty or absent)
%
% See also: AIR_CAVITY_CHARM, AIR_CAVITY_METRICS, ACOUSTIC_ANALYSIS, BWCONNCOMP_3D

arguments
    parameters   (1,1) struct
    segmentation {mustBeNumericOrLogical}
    pseudoCT     {mustBeNumeric} = []
    charm_air    {mustBeNumericOrLogical} = []
end

    cfg = get_field(parameters, {'headmodel', 'air_cavities'}, struct());
    source = char(get_field(cfg, {'source'}, 'none'));
    info = struct('source', source, 'n_vox', 0, 'n_components', 0, 'note', '');
    mask = [];

    if strcmp(source, 'none')
        return
    end
    if ~strcmp(get_field(parameters, {'simulation', 'medium'}, ''), 'layered')
        info.note = 'air-cavity check needs simulation.medium: layered';     % water, phantom: silent
        return
    end
    if isempty(segmentation)
        info.note = 'no segmentation on the grid';
        fprintf('[air_cavity_mask] Skipped: %s.\n', info.note);
        return
    end

    switch source
        case 'simnibs'
            if isempty(charm_air) || ~isequal(size(charm_air), size(segmentation))
                info.note = 'no SimNIBS charm air mask on the simulation grid (see the air_cavity_charm log line)';
                fprintf('[air_cavity_mask] Skipped: %s.\n', info.note);
                return
            end
            mask = logical(charm_air);
            if any(mask(:))
                cc = conncomp(mask);
                info.n_components = cc.NumObjects;
            end

        case 'segmentation'
            labels = double(get_field(cfg, {'labels'}, []));
            mask = ismember(segmentation, labels);
            if any(mask(:))
                cc = conncomp(mask);
                info.n_components = cc.NumObjects;
            end

        case 'hu'
            if isempty(pseudoCT) || ~isequal(size(pseudoCT), size(segmentation))
                info.note = 'no CT/pseudo-CT on the simulation grid (pct.enabled, and a backend that keeps it)';
                fprintf('[air_cavity_mask] Skipped: %s.\n', info.note);
                return
            end
            win = double(get_field(cfg, {'hu_window'}, [-1200, -400]));
            air = pseudoCT >= win(1) & pseudoCT <= win(2);      % NaN padding compares false
            % Binary median over the 3x3(x3) neighbourhood = majority vote
            if ndims(air) == 3
                k = ones(3, 3, 3);
            else
                k = ones(3, 3);
            end
            air = convn(double(air), k, 'same') > numel(k) / 2;
            cc = conncomp(air);
            n_raw = cc.NumObjects;
            if n_raw > 0
                sizes = cellfun(@numel, cc.PixelIdxList);
                [~, i_largest] = max(sizes);
                drop = false(1, n_raw);
                drop(i_largest) = true;
                if ~logical(get_field(cfg, {'keep_boundary'}, false))
                    edge = false(size(air));
                    edge(1, :, :) = true; edge(end, :, :) = true;
                    edge(:, 1, :) = true; edge(:, end, :) = true;
                    if ndims(air) == 3
                        edge(:, :, 1) = true; edge(:, :, end) = true;
                    end
                    drop = drop | cellfun(@(ix) any(edge(ix)), cc.PixelIdxList);
                end
                air(vertcat(cc.PixelIdxList{drop})) = false;
                info.n_components = nnz(~drop);
            end
            mask = air;
            fprintf('[air_cavity_mask] HU window [%g, %g]: %d component(s), %d kept as internal air.\n', ...
                win(1), win(2), n_raw, info.n_components);

        otherwise
            error('air_cavity_mask:badSource', ...
                'headmodel.air_cavities.source must be ''none'', ''simnibs'', ''segmentation'' or ''hu''; got ''%s''.', source);
    end

    info.n_vox = nnz(mask);
    vox_ml = prod(parameters.grid.resolution_mm * ones(1, ndims(mask))) / 1e3;
    if info.n_vox == 0
        info.note = 'no cavity voxels on the simulation grid';
        fprintf('[air_cavity_mask] No air-cavity voxels on the simulation grid (source: %s).\n', source);
    else
        fprintf('[air_cavity_mask] %d air-cavity voxels (%.2f mL; source: %s).\n', ...
            info.n_vox, info.n_vox * vox_ml, source);
    end
end

%% Helpers ---------------------------------------------------------------------

function cc = conncomp(bw)
% 26-connected (2D: 8-connected) components: Image Processing Toolbox when licensed (fast on the
% large air volume around the head), else BWCONNCOMP_3D
    use_ipt = license('test', 'image_toolbox') && ~isempty(which('bwconncomp'));
    if use_ipt
        try
            use_ipt = logical(license('checkout', 'image_toolbox'));
        catch
            use_ipt = false;
        end
    end
    if use_ipt
        conn = 26;
        if ismatrix(bw), conn = 8; end             % the 2D equivalent of 26
        cc = bwconncomp(bw, conn);
    else
        cc = bwconncomp_3d(bw, 26);
    end
end
