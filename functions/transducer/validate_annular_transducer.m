function tr = validate_annular_transducer(tr, t_i)
% VALIDATE_ANNULAR_TRANSDUCER  Validate configuration fields of an annular transducer
%
% Checks that all required fields for an annular transducer definition are
% present and consistent in the input transducer struct.
%
% Use as:
%   tr = validate_annular_transducer(tr, t_i)
%
% Input:
%   tr  - (1,1) transducer definition struct
%   t_i - transducer index used in error messages
%
% Output:
%   tr - validated transducer struct (unchanged if all checks pass)
%
% See also: VALIDATE_MATRIX_TRANSDUCER, TRANSDUCER_SETUP

    assert(isfield(tr, 'annular'), ...
        'Transducer %i; Appropriate configuration for annular transducer is missing.', t_i);

    annular_tr = tr.annular;
    
    % Validate required geometric parameters
    assert(isfield(annular_tr, 'elem_id_mm'), ...
        'Transducer %i; Missing elem_id_mm for annular transducer. This defines inner diameters of elements.', t_i);
    
    assert(isfield(annular_tr, 'elem_od_mm'), ...
        'Transducer %i; Missing elem_od_mm for annular transducer. This defines outer diameters of elements.', t_i);
    
    assert(isfield(annular_tr, 'elem_n'), ...
        'Transducer %i; Missing elem_n parameter for annular transducer. This defines element count.', t_i);
    
    assert(numel(annular_tr.elem_id_mm) == annular_tr.elem_n, ...
        'Transducer %i; elem_id_mm length must match elem_n.', t_i);
            
    assert(numel(annular_tr.elem_od_mm) == annular_tr.elem_n, ...
        'Transducer %i; elem_od_mm length must match elem_n.', t_i);

    % Validate inner/outer diameter ordering
    assert(all(annular_tr.elem_od_mm > annular_tr.elem_id_mm), ...
        'Transducer %i; Outer diameter must be larger than inner diameter for all elements.', t_i);
    
    % Validate curvature definition
    assert(isfield(annular_tr, 'curv_radius_mm'), ...
        'Transducer %i; Missing curv_radius_mm field for annular transducer. Please specify radius of curvature.', t_i);

    % 3D steering unavailable for annular arrays → align with focus
    tr.align_to_focus = true;

    % Flat elements (curv_radius_mm = Inf): no geometric focus, the face is the
    % exit plane (see annular_bowl_depth_mm)
    if isinf(annular_tr.curv_radius_mm)
        if ~isfield(annular_tr, 'dist_geom_ep_mm') || isempty(annular_tr.dist_geom_ep_mm)
            annular_tr.dist_geom_ep_mm = Inf;
        end

    % Calculate distance to transducer plane if not provided
    elseif ~isfield(annular_tr, 'dist_geom_ep_mm') || isempty(annular_tr.dist_geom_ep_mm)
        assert(annular_tr.curv_radius_mm > max(annular_tr.elem_od_mm)/2, ...
            'Transducer %i; curv_radius_mm must exceed aperture radius.', t_i);

        annular_tr.dist_geom_ep_mm = sqrt(annular_tr.curv_radius_mm^2 - ...
                                (max(annular_tr.elem_od_mm) / 2)^2);

        fprintf('Transducer %i; Distance to transducer plane is not provided, calculated as %.2f mm\n', ...
                t_i, annular_tr.dist_geom_ep_mm);
    end

    % Evaluate source phase expressions if stored as cell arrays
    if isfield(annular_tr, 'elem_phase_rad') && iscell(annular_tr.elem_phase_rad)
        for p_i = 1:numel(annular_tr.elem_phase_rad)
            if ~isnumeric(annular_tr.elem_phase_rad{p_i})
                annular_tr.elem_phase_rad{p_i} = eval_phase_expression(annular_tr.elem_phase_rad{p_i});
            end
        end
        annular_tr.elem_phase_rad = cell2mat(annular_tr.elem_phase_rad);
    end 

    if isfield(annular_tr, 'elem_phase_deg') && iscell(annular_tr.elem_phase_deg)
        for p_i = 1:numel(annular_tr.elem_phase_deg)
            if ~isnumeric(annular_tr.elem_phase_deg{p_i})
                annular_tr.elem_phase_deg{p_i} = eval_phase_expression(annular_tr.elem_phase_deg{p_i});
            end
        end
        annular_tr.elem_phase_deg = cell2mat(annular_tr.elem_phase_deg);
    end 

    % Validate that source phase is set in radians or degrees
    if ~isfield(annular_tr, 'elem_phase_rad') && isfield(annular_tr, 'elem_phase_deg')
        annular_tr.elem_phase_rad = deg2rad(annular_tr.elem_phase_deg);
    elseif ~isfield(annular_tr, 'elem_phase_deg') && isfield(annular_tr, 'elem_phase_rad')
        annular_tr.elem_phase_deg = rad2deg(annular_tr.elem_phase_rad);
    elseif ~isfield(annular_tr, 'elem_phase_rad') && ~isfield(annular_tr, 'elem_phase_deg')
        error('Transducer %i; Phase must be specified as elem_phase_rad or elem_phase_deg.', t_i);
    elseif ~isempty(annular_tr.elem_phase_deg) && isnumeric(annular_tr.elem_phase_deg) && ...
            isnumeric(annular_tr.elem_phase_rad) && ...
            (numel(annular_tr.elem_phase_rad) ~= numel(annular_tr.elem_phase_deg) || ...
             any(abs(angle(exp(1i * (annular_tr.elem_phase_rad(:) - deg2rad(annular_tr.elem_phase_deg(:)))))) > 1e-6))
        % Both set and different (e.g. elem_phase_deg changed after loading,
        % leaving the derived elem_phase_rad stale): the source uses
        % elem_phase_rad, so make it follow elem_phase_deg
        warning('prestus:phaseMismatch', ...
            ['Transducer %i; elem_phase_rad does not match elem_phase_deg; using elem_phase_deg. ' ...
             'Set one of them, or both consistently.'], t_i);
        annular_tr.elem_phase_rad = deg2rad(annular_tr.elem_phase_deg);
    end 

    % Backward-compatibility shim: rename deprecated elem_amp → elem_p_amp
    if isfield(annular_tr, 'elem_amp') && (~isfield(annular_tr, 'elem_p_amp') || isempty(annular_tr.elem_p_amp))
        warning('prestus:deprecatedField', ...
            'annular.elem_amp is deprecated; rename to elem_p_amp in your config.');
        annular_tr.elem_p_amp = annular_tr.elem_amp;
        annular_tr = rmfield(annular_tr, 'elem_amp');
    end

    % Ensure source amplitude matches number of transducer elements
    if numel(annular_tr.elem_p_amp) == 1 && annular_tr.elem_n > 1
        annular_tr.elem_p_amp = repmat(annular_tr.elem_p_amp, [1, annular_tr.elem_n]);
    end

    % Default depth_mm for visualization (if not set in config).
    % 16 mm is a representative housing depth used only for 2-D overlay plots
    % (how far behind the exit plane the transducer body extends); it has no
    % effect on the simulation itself.  Override via annular.depth_mm in your config.
    if ~isfield(annular_tr, 'depth_mm') || isempty(annular_tr.depth_mm)
        annular_tr.depth_mm = 16;
    end

    % Encode updated annular transducer field
    tr.annular = annular_tr;
   
end

function v = eval_phase_expression(expr)
% Evaluate a config phase expression such as '274*pi/180'. Only numbers, pi
% and arithmetic are allowed, so a config value can never run other code.
    expr = char(string(expr));
    if isempty(regexp(regexprep(expr, '\<pi\>', ''), '^[0-9.eE+\-*/^()\[\],; \t]*$', 'once'))
        error('validate_annular_transducer:badPhaseExpression', ...
            'Phase "%s" must be a number or an arithmetic expression of numbers and pi.', expr);
    end
    v = eval(expr);
end
