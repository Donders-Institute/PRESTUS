function ld = simnibs_ld_library_path(parameters)
% SIMNIBS_LD_LIBRARY_PATH  LD_LIBRARY_PATH prefix for SimNIBS tools (charm, mni2subject, ...)
%
% hpc.ld_library_path when set; otherwise the libraries of the SimNIBS install
% in startup.simnibs_bin_path (<simnibs package>/external/lib/linux, the folder
% SimNIBS puts on its own search path), so the library path follows the selected
% SimNIBS environment. '' when neither exists (e.g. on macOS or Windows).
%
% Use as:
%   ld = simnibs_ld_library_path(parameters)
%
% Input:
%   parameters - PRESTUS config (hpc.ld_library_path, startup.simnibs_bin_path)
%
% Output:
%   ld - folder to prepend to LD_LIBRARY_PATH, or ''
%
% See also: SIMNIBS_INSTALL_INFO, SEGMENTATION_RUN, CONVERT_FINAL_TO_MNI_SIMNIBS

    ld = '';
    if isfield(parameters, 'hpc') && isfield(parameters.hpc, 'ld_library_path') ...
            && ~isempty(parameters.hpc.ld_library_path)
        ld = char(parameters.hpc.ld_library_path);
        return
    end
    if ispc || ~isfield(parameters, 'startup') || ~isfield(parameters.startup, 'simnibs_bin_path') ...
            || isempty(parameters.startup.simnibs_bin_path)
        return
    end
    info = simnibs_install_info(char(parameters.startup.simnibs_bin_path));
    ld = info.ld_library_path;
end
