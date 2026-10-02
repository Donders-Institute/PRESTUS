function tf = telemetry_has_display()
% TELEMETRY_HAS_DISPLAY  True when MATLAB can render a figure/dialog window.
%
% False for compiled/deployed applications and for any session without a
% display attached (e.g. SLURM/PBS compute nodes, matlab -nodisplay). True
% for a normal desktop session and for `matlab -batch` run locally with a
% display attached — such sessions have no interactive stdin for input(),
% but can still show a modal dialog via uiconfirm.
%
% Use as:
%   tf = telemetry_has_display()
%
% Output:
%   tf - logical
%
% See also: TELEMETRY_SETUP

    tf = false;
    try
        if isdeployed
            return
        end
        tf = usejava('jvm') && feature('ShowFigureWindows');
    catch
        tf = false;
    end
end
