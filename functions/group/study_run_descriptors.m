function d = study_run_descriptors(parameters)
% STUDY_RUN_DESCRIPTORS  Settings that tell one run of a study apart from the others
%
% Collects the per-run settings that typically differ between subjects of a
% study (placement mode, target, transducer, expected focal depth, target
% intensity), for the study plan printed by PRESTUS_STUDY_START and for the
% per-run targeting sidecar (TARGETING_ANALYSIS) read by the group report.
%
% Use as:
%   d = study_run_descriptors(parameters)
%
% Input:
%   parameters - loaded PRESTUS config of one run
%
% Output:
%   d - struct with fields
%         placement_mode        placement.mode ('' when absent)
%         target_label          target name of the placement mode: localite
%                               target, heuristic/plantus target_name; '' otherwise
%         transducer_name       transducer(1).name (or serial)
%         n_transducers         number of transducers
%         focal_distance_ep_mm  expected focal distance from the exit plane, per transducer
%         target_isppa_wcm2     target free-water Isppa of transducer(1) (may be a vector)
%         condition             parameters.study.condition ('' outside a study)
%
% See also: PRESTUS_STUDY_START, TARGETING_ANALYSIS

    d = struct('placement_mode', "", 'target_label', "", 'transducer_name', "", ...
        'n_transducers', 0, 'focal_distance_ep_mm', [], 'target_isppa_wcm2', [], 'condition', "");

    mode = '';
    if isfield(parameters, 'placement') && isfield(parameters.placement, 'mode')
        mode = text_of(parameters.placement.mode);
    end
    d.placement_mode = string(mode);

    pl = struct();
    if isfield(parameters, 'placement'); pl = parameters.placement; end
    switch mode
        case 'localite'
            d.target_label = string(field_text(pl, {'localite', 'target'}));
            if d.target_label == ""
                pos = field_value(pl, {'localite', 'position'});
                if isnumeric(pos) && ~isempty(pos)
                    d.target_label = "series " + string(pos(1));
                end
            end
        case 'heuristic'
            d.target_label = string(field_text(pl, {'heuristic', 'target_name'}));
        case 'plantus'
            d.target_label = string(field_text(pl, {'plantus', 'target_name'}));
        case 'mni'
            v = field_value(pl, {'mni', 'focus_pos_mm'});
            if isnumeric(v) && numel(v) == 3
                d.target_label = string(sprintf('MNI [%g %g %g]', v));
            end
    end

    if isfield(parameters, 'transducer') && ~isempty(parameters.transducer)
        tr = parameters.transducer;
        d.n_transducers = numel(tr);
        name = '';
        if isfield(tr(1), 'name'); name = text_of(tr(1).name); end
        if isempty(name) && isfield(tr(1), 'serial'); name = text_of(tr(1).serial); end
        d.transducer_name = string(name);
        fd = nan(1, numel(tr));
        for k = 1:numel(tr)
            if isfield(tr(k), 'focal_distance_ep') && isnumeric(tr(k).focal_distance_ep) && ...
                    ~isempty(tr(k).focal_distance_ep)
                fd(k) = tr(k).focal_distance_ep(1);
            end
        end
        d.focal_distance_ep_mm = fd;
        if isfield(tr(1), 'target_isppa_wcm2') && isnumeric(tr(1).target_isppa_wcm2)
            d.target_isppa_wcm2 = tr(1).target_isppa_wcm2(:)';
        end
    end

    if isfield(parameters, 'study') && isfield(parameters.study, 'condition')
        d.condition = string(text_of(parameters.study.condition));
    end
end

% -------------------------------------------------------------------------
function v = field_value(s, path)
    v = [];
    for p = 1:numel(path)
        if ~isstruct(s) || ~isfield(s, path{p}); return; end
        s = s.(path{p});
    end
    v = s;
end

function t = field_text(s, path)
    t = text_of(field_value(s, path));
end

function t = text_of(v)
    if ischar(v)
        t = v;
    elseif isstring(v) && isscalar(v)
        t = char(v);
    else
        t = '';
    end
end
