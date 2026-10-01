function varargout = kwave_check_refinement(varargin)
% KWAVE_CHECK_REFINEMENT  Deprecated name of KWAVE_CALIBRATE_DRIVE
%
% The function calibrates the drive of per-depth refinements in free-water
% k-Wave (and reports the agreement with the measured profiles); it was
% renamed to say so. This wrapper forwards all arguments.
%
% See also: KWAVE_CALIBRATE_DRIVE

    warning('kwave_check_refinement:deprecated', ...
        'kwave_check_refinement is deprecated; use kwave_calibrate_drive (same arguments).');
    [varargout{1:nargout}] = kwave_calibrate_drive(varargin{:});
end
