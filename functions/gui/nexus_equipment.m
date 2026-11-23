function out = nexus_equipment(args)
% NEXUS_EQUIPMENT  Equipment-editor checks for the PRESTUS Nexus (through NEXUS_BRIDGE)
%
% The Equipment sub-tab (Transducer tab, advanced) writes a new or edited
% transducer entry, and any new driving systems, to a scratch folder and asks
% for CHECK_EQUIPMENT_CONFIG's report on them before saving, as CALIBRATION_GUI
% does for its 'New…' form (validate_equipment_entry). The check warns: it also
% flags some shipped entries (touching sub-rings, max_foc beyond the curvature
% radius), so the Nexus offers to save anyway.
%
% Use as:
%   out = nexus_equipment(args)
%
% Input:
%   args - struct with .action:
%     'check'  .folder: folder holding the entry files (tran-*_transducer.yaml,
%              ds-*_driver.yaml)
%     'layout' .matrix: a transducer.matrix block (define_here grid)
%
% Output:
%   out - 'check': .issues (cellstr, 'serial: issue' per line; empty when valid)
%         'layout': .x_mm, .y_mm, .z_mm (element centres, apex at the origin, beam
%                   along +z; from CONVERT_TO_ELEMENT_POS, as the simulation places
%                   them and as CALIBRATION_GUI's preview), .n (element count),
%                   .error ('' or why the layout failed)
%
% See also: CHECK_EQUIPMENT_CONFIG, CALIBRATION_GUI, NEXUS_BRIDGE

arguments
    args (1,1) struct
end

    switch args.action
        case 'check'
            [~, rep] = evalc('check_equipment_config(char(args.folder))');
            issues = {};
            for k = 1:numel(rep)
                for i = 1:numel(rep(k).issues)
                    issues{end+1} = sprintf('%s: %s', rep(k).serial, rep(k).issues{i}); %#ok<AGROW>
                end
            end
            out = struct('issues', {issues});
        case 'layout'
            out = struct('x_mm', [], 'y_mm', [], 'z_mm', [], 'n', 0, 'error', '');
            try
                m = args.matrix;
                if isfield(m, 'is_curved'); m.is_curved = logical(m.is_curved); end
                [~, p] = evalc('convert_to_element_pos(struct(), struct(''matrix'', m), [0; 0; 0], [0; 0; 0])');
                pos = p' * 1e3;
                out.x_mm = pos(:, 1)'; out.y_mm = pos(:, 2)'; out.z_mm = pos(:, 3)';
                out.n = size(pos, 1);
            catch ME
                out.error = ME.message;
            end
        otherwise
            error('nexus_equipment:action', 'Unknown action ''%s''.', args.action);
    end
end
