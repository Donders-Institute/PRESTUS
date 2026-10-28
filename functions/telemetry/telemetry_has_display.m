function tf = telemetry_has_display()
% TELEMETRY_HAS_DISPLAY  True when MATLAB can render a figure/dialog window.
%
% False for compiled/deployed applications and for any session without a
% display attached (e.g. SLURM/PBS compute nodes, matlab -nodisplay). True
% for a normal desktop session and for `matlab -batch` run locally with a
% display attached. A display alone does not make a session interactive:
% TELEMETRY_SETUP never prompts in a -batch session (see
% TELEMETRY_IS_NONINTERACTIVE).
%
% Use as:
%   tf = telemetry_has_display()
%
% Output:
%   tf - logical
%
% See also: TELEMETRY_SETUP, TELEMETRY_IS_NONINTERACTIVE

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
