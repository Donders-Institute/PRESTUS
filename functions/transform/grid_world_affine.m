function A = grid_world_affine(src)
% GRID_WORLD_AFFINE  Simulation-grid voxel -> world (RAS mm) affine of a run
%
% The full simulation grid of a run maps to world space through its planning
% image: sim_voxel -> world_RAS = planimg.inv_transf * planimg.t1_header.Transform.T
% (row-vector convention: [world 1] = [voxel 1] * A). The planning image comes
% either directly or from the run's planimg sidecar, which GRID_STAGE writes
% for every preprocessing run (sub-XXX_<medium>_planimg<affix>.mat, affix =
% io.preproc_affix, else io.output_affix).
%
% Use as:
%   A = grid_world_affine(planimg)
%   A = grid_world_affine(acoustic_info)   % .parameters of an acoustic cache
%
% Input:
%   src - planimg struct (.inv_transf, .t1_header.Transform), or a struct
%         with .parameters (io.dir_cache, subject_id, simulation.medium,
%         io.preproc_affix or io.output_affix) whose planimg sidecar is read
%
% Output:
%   A   - 4x4 affine, or [] when it cannot be recovered (no sidecar, or no
%         T1 header, e.g. a phantom or water-only run)
%
% See also: GRID_REMAP_AFFINE, REMAP_ACOUSTIC_FIELD, COMBINE_ASYNC_INTENSITY

    A = [];
    if ~isstruct(src)
        return;
    end
    if isfield(src, 'inv_transf')
        planimg = src;
    else
        planimg = load_sidecar(src);
        if isempty(planimg)
            return;
        end
    end
    if ~isfield(planimg, 'inv_transf') || isempty(planimg.inv_transf) || ...
            ~isfield(planimg, 't1_header') || isempty(planimg.t1_header) || ...
            ~isfield(planimg.t1_header, 'Transform')
        return;
    end
    A = planimg.inv_transf * planimg.t1_header.Transform.T;
end

function planimg = load_sidecar(acoustic_info)
    planimg = [];
    if ~isfield(acoustic_info, 'parameters')
        return;
    end
    p = acoustic_info.parameters;
    if ~isfield(p, 'io') || ~isfield(p, 'subject_id') || ~isfield(p, 'simulation') || ...
            ~isfield(p.io, 'dir_cache')
        return;
    end
    if isfield(p.io, 'preproc_affix') && ~isempty(p.io.preproc_affix)
        affix = p.io.preproc_affix;
    elseif isfield(p.io, 'output_affix')
        affix = p.io.output_affix;
    else
        return;
    end
    sidecar = fullfile(p.io.dir_cache, ...
        sprintf('sub-%03d_%s_planimg%s.mat', p.subject_id, p.simulation.medium, affix));
    if ~isfile(sidecar)
        return;
    end
    s = load(sidecar, 'planimg');
    if isfield(s, 'planimg')
        planimg = s.planimg;
    end
end
