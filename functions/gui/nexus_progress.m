function nexus_progress(msg)
% NEXUS_PROGRESS  Report progress of the current NEXUS_BRIDGE request to the Nexus
%
% Writes msg to progress/<id>.txt of the running request, which the Nexus
% shows while it waits. Outside the bridge it does nothing, so façade
% functions can call it unconditionally.
%
% Use as:
%   nexus_progress('Depth 45.0 mm (2/5)')
%
% Input:
%   msg - one line of text
%
% See also: NEXUS_BRIDGE

    f = getappdata(0, 'nexus_bridge_progress');
    if isempty(f); return; end
    fid = fopen([f '.tmp'], 'w', 'n', 'UTF-8');
    if fid < 0; return; end
    fwrite(fid, char(msg), 'char');
    fclose(fid);
    movefile([f '.tmp'], f, 'f');
end
