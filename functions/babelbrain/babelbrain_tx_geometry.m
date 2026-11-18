function geom = babelbrain_tx_geometry(tr, tx_system, bb_root)
% BABELBRAIN_TX_GEOMETRY  Transducer geometry BabelBrain simulates for a PRESTUS transducer
%
% With an empty tx_system, BabelBrain simulates the PRESTUS transducer itself:
% its ring geometry (transducer.annular) on BabelBrain's generic focused
% annular template (tx_system 'PRESTUS' in run_babelbrain_case.py). The drive
% from the PRESTUS calibration library, fitted on the same geometry, then maps
% ring for ring. Single-element bowls are one ring. A named BabelBrain preset
% (e.g. 'CTX_500', or a phased array such as 'H317') uses BabelBrain's own
% device definition from its registry instead.
%
% Use as:
%   geom = babelbrain_tx_geometry(tr, tx_system, bb_root)
%
% Input:
%   tr        - PRESTUS transducer struct (parameters.transducer(1))
%   tx_system - parameters.babelbrain.tx_system ('' = PRESTUS geometry)
%   bb_root   - BabelBrain source root (contains TranscranialModeling/)
%
% Output:
%   geom - struct with fields
%       .tx_system      - value passed to run_babelbrain_case.py --tx_system
%       .type           - BabelBrain geometry family (e.g. 'focused_annular_array')
%       .yaml           - preset default.yaml path ('' for PRESTUS geometry)
%       .focal_length_m - focal length (radius of curvature) [m]; [] if unknown
%       .fl_src         - where focal_length_m came from (for the log)
%       .flags          - extra command-line arguments for run_babelbrain_case.py
%       .cfg            - preset default.yaml contents (empty struct for PRESTUS geometry)
%
% See also: BABELBRAIN_WRAPPER, EXPORT_BABELBRAIN_INPUTS, BABELBRAIN_PRESET

    arguments
        tr        (1,1) struct
        tx_system       = ''
        bb_root   (1,:) char = ''
    end
    tx_system = char(tx_system);

    has_rings = isfield(tr, 'annular') && ...
        all(isfield(tr.annular, {'elem_id_mm', 'elem_od_mm', 'curv_radius_mm'})) && ...
        ~isempty(tr.annular.elem_od_mm) && ~isempty(tr.annular.curv_radius_mm);

    if isempty(tx_system) || strcmp(tx_system, 'PRESTUS')
        if ~has_rings
            error('PRESTUS:babelbrain:noGeometry', ...
                ['BabelBrain simulates the PRESTUS transducer from transducer.annular ' ...
                 '(elem_id_mm, elem_od_mm, curv_radius_mm), which this transducer lacks. ' ...
                 'Set babelbrain.tx_system to a BabelBrain preset (e.g. a phased array such as H317).']);
        end
        id_m = double(tr.annular.elem_id_mm(:)') / 1000;
        od_m = double(tr.annular.elem_od_mm(:)') / 1000;
        rc_m = double(tr.annular.curv_radius_mm) / 1000;
        if numel(id_m) ~= numel(od_m)
            error('PRESTUS:babelbrain:ringMismatch', ...
                'transducer.annular.elem_id_mm and elem_od_mm need one value per ring (%d vs %d).', ...
                numel(id_m), numel(od_m));
        end
        aperture_m = max(od_m);   % radiating aperture; sizes BabelBrain's source face in the domain
        list = @(v) strjoin(arrayfun(@(x) sprintf('%.7f', x), v, 'UniformOutput', false), ',');
        geom = struct( ...
            'tx_system', 'PRESTUS', ...
            'type', 'focused_annular_array', ...
            'yaml', '', ...
            'focal_length_m', rc_m, ...
            'fl_src', 'transducer.annular.curv_radius_mm', ...
            'flags', {{'--aperture', sprintf('%.7f', aperture_m), ...
                       '--focal_length', sprintf('%.7f', rc_m), ...
                       '--in_diameters', list(id_m), ...
                       '--out_diameters', list(od_m)}}, ...
            'cfg', struct());
        fprintf(['[babelbrain_tx_geometry] PRESTUS transducer geometry: %d ring(s), ' ...
            'R_curv %.2f mm, OD %.2f mm\n'], numel(od_m), rc_m*1e3, aperture_m*1e3);
        return
    end

    % Named BabelBrain preset
    preset = babelbrain_preset(bb_root, tx_system);
    geom = struct('tx_system', preset.name, 'type', preset.type, 'yaml', preset.yaml, ...
        'focal_length_m', [], 'fl_src', 'BabelBrain preset yaml', 'flags', {{}}, 'cfg', preset.cfg);
    if isfield(preset.cfg, 'FocalLength')
        geom.focal_length_m = preset.cfg.FocalLength;
    end

    % Single-element presets (Single, BSonix) have no top-level geometry in
    % default.yaml: BabelBrain's GUI asks for aperture and focal length, and the
    % simple_focused RUN_SIM.RunCases reads both from its kwargs. Take them from
    % the PRESTUS transducer.
    if strcmp(preset.type, 'simple_focused')
        if ~has_rings
            error('PRESTUS:babelbrain:singleGeometry', ...
                ['BabelBrain preset ''%s'' needs the aperture and focal length from PRESTUS: ' ...
                 'set transducer.annular.curv_radius_mm and transducer.annular.elem_od_mm.'], preset.name);
        end
        if numel(tr.annular.elem_od_mm) > 1
            warning('PRESTUS:babelbrain:singleMultiElem', ...
                ['BabelBrain preset ''%s'' models one solid bowl; the %d-element transducer is ' ...
                 'simulated as a single element of outer diameter %.1f mm.'], ...
                preset.name, numel(tr.annular.elem_od_mm), max(tr.annular.elem_od_mm));
        end
        if any(tr.annular.elem_id_mm > 0)
            warning('PRESTUS:babelbrain:singleHole', ...
                'BabelBrain preset ''%s'' has no central opening; elem_id_mm %.1f mm is ignored.', ...
                preset.name, min(tr.annular.elem_id_mm));
        end
        aperture_m = max(tr.annular.elem_od_mm) / 1000;
        if isempty(geom.focal_length_m)
            geom.focal_length_m = tr.annular.curv_radius_mm / 1000;
            geom.fl_src = 'transducer.annular.curv_radius_mm';
        end
        geom.flags = {'--aperture', sprintf('%.6f', aperture_m), ...
                      '--focal_length', sprintf('%.6f', geom.focal_length_m)};
    end
end
