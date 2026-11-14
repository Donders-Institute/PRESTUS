function fp = equipment_geometry_fingerprint(tran)
% EQUIPMENT_GEOMETRY_FINGERPRINT  Short hash of the transducer geometry a calibration was fitted for
%
% Calibration entries record it in meta.geometry_fingerprint, and
% LOAD_TRANSDUCER_FROM_LIBRARY warns when an entry is used with a transducer
% whose geometry differs (as BabelBrain checks the TxConfig stored in a
% calibration file). The hash covers the frequency and, for annular arrays,
% the ring diameters, radius of curvature and exit-plane distance, rounded
% to 1 µm (0.1 Hz for the frequency).
%
% Use as:
%   fp = equipment_geometry_fingerprint(tran)
%
% Input:
%   tran - equipment entry (load_equipment_config) or its .transducer block
%
% Output:
%   fp   - 12-character hex string; '' when the geometry is not annular
%
% See also: LOAD_TRANSDUCER_FROM_LIBRARY, LOAD_EQUIPMENT_CONFIG

    if isfield(tran, 'transducer'); tran = tran.transducer; end
    fp = '';
    if ~isfield(tran, 'annular') || ~isstruct(tran.annular)
        return;
    end
    a = tran.annular;
    need = {'elem_id_mm', 'elem_od_mm', 'curv_radius_mm', 'dist_geom_ep_mm'};
    if ~all(isfield(a, need)) || ~isfield(tran, 'freq_hz')
        return;
    end
    parts = {sprintf('f=%.1f', double(tran.freq_hz)), ...
             ['id=' sprintf('%.3f,', double(a.elem_id_mm(:)'))], ...
             ['od=' sprintf('%.3f,', double(a.elem_od_mm(:)'))], ...
             sprintf('roc=%.3f', double(a.curv_radius_mm)), ...
             sprintf('ep=%.3f', double(a.dist_geom_ep_mm))};
    md = java.security.MessageDigest.getInstance('SHA-256');
    h = typecast(md.digest(uint8(strjoin(parts, ';'))), 'uint8');
    fp = lower(reshape(dec2hex(h(1:6))', 1, []));
end
