function v = telemetry_notice_version()
% TELEMETRY_NOTICE_VERSION  Version of the telemetry notice that consent is given to
%
% The consent record (~/.prestus/telemetry.json) stores the notice version it
% was given to (field 'notice'; records without it are version 1). Fields a
% later notice adds are sent only with consent to that notice or later:
%   1 - the original notice
%   2 - (2026-11-17) adds the medians of the modelled skull sound speed,
%       density, and attenuation, in all skull voxels and on the beam path
%   3 - (2026-11-23) consent per level (TELEMETRY_LEVEL_INFO); adds the
%       'outcomes' level: rounded brain intensity and transmission, brain MI,
%       targeting error, brain and skull temperature rise, duty cycle and
%       sonication time
% The Nexus writes the same version (prestus_nexus/telemetry.py, NOTICE_VERSION).
%
% Use as:
%   v = telemetry_notice_version()
%
% Output:
%   v - current notice version (integer)
%
% See also: TELEMETRY_SETUP, TELEMETRY_SET_CONSENT, TRACK_USAGE
v = 3;
end
