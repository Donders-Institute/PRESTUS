function make_synthetic_localite_export(raw_root, sub, ses, poses, shift, names)
% MAKE_SYNTHETIC_LOCALITE_EXPORT  Write a synthetic raw Localite export for one session
%
% Writes <raw_root>/<sub>/<ses>/loc/Session_20240101120000000/ with a
% TriggerMarkers_Coil0 file holding one 80-pulse series per pose and an
% EntryTarget file with one named target per pose, 60 mm along the beam
% axis (matrix column 1) and 1 mm off it.
%
% Use as:
%   make_synthetic_localite_export(raw_root, sub, ses, poses, shift, names)
%
% Input:
%   raw_root - path.localite_raw
%   sub, ses - e.g. 'sub-001', 'ses-01'
%   poses    - cell of 4x4 Localite matrices, one series each
%   shift    - [3 x 1] mm added to every position (so sessions differ)
%   names    - cell of planned-target labels, one per pose
%
% See also: TEST_NEURONAV_BATCH, NEURONAV_DETECT_SERIES

    base = fullfile(raw_root, sub, ses, 'loc', 'Session_20240101120000000');
    mkdir(fullfile(base, 'TMSTrigger'));
    mkdir(fullfile(base, 'EntryTarget'));
    fid = fopen(fullfile(base, 'TMSTrigger', 'TriggerMarkers_Coil0_20240101123000000.xml'), 'w');
    fprintf(fid, '<?xml version="1.0" encoding="UTF-8"?>\n<TriggerMarkerList coordinateSpace="RAS">\n');
    t = 0;
    for s = 1:numel(poses)
        M = poses{s}; M(1:3, 4) = M(1:3, 4) + shift;
        for k = 1:80
            t = t + 1000;
            fprintf(fid, '  <TriggerMarker description="" recordingTime="%d">\n    <Matrix4D', t);
            for r = 0:3
                for c = 0:3
                    fprintf(fid, ' data%d%d="%.10g"', r, c, M(r+1, c+1));
                end
            end
            fprintf(fid, '/>\n  </TriggerMarker>\n');
        end
        t = t + 60000;
    end
    fprintf(fid, '</TriggerMarkerList>\n');
    fclose(fid);

    fid = fopen(fullfile(base, 'EntryTarget', 'EntryTarget20240101123100000.xml'), 'w');
    fprintf(fid, '<?xml version="1.0" encoding="UTF-8"?>\n<EntryTargetList>\n');
    for s = 1:numel(poses)
        M = poses{s};
        xyz = M(1:3, 4) + shift + 60 * M(1:3, 1) + [0; 0; 1];   % 1 mm off axis
        fprintf(fid, ['  <Entry><Marker description="" set="false"><ColVec3D data0="0" data1="0" data2="0"/></Marker></Entry>\n' ...
            '  <Target index="%d"><Marker description="%s" set="true"><ColVec3D data0="%.6f" data1="%.6f" data2="%.6f"/></Marker></Target>\n'], ...
            s - 1, names{s}, xyz);
    end
    fprintf(fid, '</EntryTargetList>\n');
    fclose(fid);
end
