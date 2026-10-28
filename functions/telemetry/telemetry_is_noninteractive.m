function tf = telemetry_is_noninteractive()
% TELEMETRY_IS_NONINTERACTIVE  True when no telemetry prompt may be shown.
%
% Such sessions get the printed notice only: no console prompt, no dialog,
% nothing recorded, nothing sent unless consent is already recorded (see
% TELEMETRY_SETUP). True for:
%   - deployed apps and SLURM/PBS jobs (TELEMETRY_IS_HEADLESS_CLUSTER);
%   - `matlab -batch` sessions (batchStartupOptionUsed, R2019a+), whether
%     or not a display is attached: an unattended run (cron, CI, a script
%     started from a terminal) must never wait for an answer;
%   - sessions without the desktop and without a display (e.g.
%     `matlab -nodisplay -r ...` from cron), where stdin may not be a
%     terminal. An interactive `matlab -nodisplay` terminal session also
%     lands here; record a decision there with TELEMETRY_SET_CONSENT.
%
% On releases before R2019a, batchStartupOptionUsed does not exist; the
% call is guarded and the desktop/display check still applies.
%
% Use as:
%   tf = telemetry_is_noninteractive()
%
% Output:
%   tf - logical
%
% See also: TELEMETRY_SETUP, TELEMETRY_IS_HEADLESS_CLUSTER, TELEMETRY_HAS_DISPLAY

    tf = true;
    try
        tf = telemetry_is_headless_cluster() || is_batch_session() || ...
            (~usejava('desktop') && ~telemetry_has_display());
    catch
        tf = true;
    end
end

% -------------------------------------------------------------------------
function tf = is_batch_session()
    try
        tf = logical(batchStartupOptionUsed);
    catch
        tf = false;   % pre-R2019a: function absent
    end
end
