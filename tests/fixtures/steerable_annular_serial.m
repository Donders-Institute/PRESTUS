function [serial, tran] = steerable_annular_serial(equip)
% STEERABLE_ANNULAR_SERIAL  First steerable multi-element annular device of the equipment library
%
% Tests that fit, steer or refine per element pick their device by these
% properties instead of by library order, which starts with model entries
% such as the single-element, fixed-focus BSonix models. The device also has
% no shipped calibration entry of its own (config/calibration, serial only,
% no driving system), so a calibration session on it starts from the
% geometric model.
%
% Use as:
%   [serial, tran] = steerable_annular_serial()
%   [serial, tran] = steerable_annular_serial(equip)    % from LOAD_EQUIPMENT_CONFIG
%
% Output:
%   serial - transducer serial
%   tran   - its equipment entry

    if nargin < 1 || isempty(equip)
        equip = load_equipment_config();
    end
    shipped = fullfile(get_prestus_path(), 'config', 'calibration');
    serials = fieldnames(equip.trans);
    for i = 1:numel(serials)
        t = equip.trans.(serials{i});
        if ~isfield(t, 'transducer') || ~isfield(t.transducer, 'annular') || ...
                t.transducer.annular.elem_n < 2 || ...
                ~isfield(t, 'min_foc') || ~isfield(t, 'max_foc') || ~(t.max_foc > t.min_foc)
            continue
        end
        [~, has_entry] = library_entry_file(shipped, serials{i});
        if has_entry
            continue
        end
        serial = serials{i};
        tran = t;
        return
    end
    error('steerable_annular_serial:none', ...
        'No steerable multi-element annular transducer without a shipped entry in the equipment library.');
end
