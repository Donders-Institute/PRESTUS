function d = prefdir_prestus()
% PREFDIR_PRESTUS  Return (and create) the PRESTUS per-user preferences directory.
%
% Returns ~/.prestus/ on all platforms. Creates the directory if absent.
%
% Use as:
%   d = prefdir_prestus()

    d = fullfile(char(java.lang.System.getProperty('user.home')), '.prestus');
    if ~isfolder(d)
        mkdir(d);
    end
end
