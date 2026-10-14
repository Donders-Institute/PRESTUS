function [hu, rho, source] = medium_pct_density_calibration(parameters)
% MEDIUM_PCT_DENSITY_CALIBRATION  CT->density calibration table for pct.mapping_density 'k-plan'
%
% HU depend on the CT scanner, tube voltage (kVp) and reconstruction kernel,
% so an HU->density table belongs to the CT a pseudo-CT's UTE->HU line was
% fitted to. pct.density_calibration selects the table:
%   ''  (default)     k-Plan's default table (config/ct_calibration/kplan_default.csv)
%   name              a table shipped in config/ct_calibration/ (file name without .csv),
%                     e.g. 'ucl_lowdose_80kvp' (pairs with skull_mapping 'treeby'),
%                     'cph2025' (pairs with a custom mapping of petra2density's line)
%   path              any CSV: two columns HU, density [kg/m^3]; lines starting
%                     with # are comments. Relative paths resolve against the
%                     current folder.
% The table is used piecewise linearly with linear extrapolation beyond its
% end points (FIT_PAIRWISELINEAR), as for the k-Plan default.
%
% Use as:
%   [hu, rho, source] = medium_pct_density_calibration(parameters)
%
% Input:
%   parameters - PRESTUS parameter struct (pct.density_calibration, optional)
%
% Output:
%   hu, rho - calibration points, sorted by HU (row vectors) [HU], [kg/m^3]
%   source  - resolved file path
%
% See also: MEDIUM_PCT_DENSITY, FIT_PAIRWISELINEAR

arguments
    parameters (1,1) struct
end

    sel = '';
    if isfield(parameters, 'pct') && isfield(parameters.pct, 'density_calibration')
        sel = char(parameters.pct.density_calibration);
    end
    cal_dir = fullfile(get_prestus_path(), 'config', 'ct_calibration');
    if isempty(sel)
        source = fullfile(cal_dir, 'kplan_default.csv');
    elseif isfile(fullfile(cal_dir, [sel '.csv']))
        source = fullfile(cal_dir, [sel '.csv']);
    elseif isfile(sel)
        source = sel;
    else
        shipped = dir(fullfile(cal_dir, '*.csv'));
        error('prestus:pct:densityCalibration', ['pct.density_calibration ''%s'' is neither a table in %s ' ...
            '(%s) nor an existing file.'], sel, cal_dir, strjoin(erase({shipped.name}, '.csv'), ', '));
    end

    T = readmatrix(source, 'FileType', 'text', 'CommentStyle', '#', 'Delimiter', ',');
    assert(size(T, 2) >= 2 && size(T, 1) >= 2 && all(isfinite(T(:, 1:2)), 'all'), 'prestus:pct:densityCalibration', ...
        'Density calibration %s: need at least two rows of HU, density [kg/m^3].', source);
    [hu, o] = sort(T(:, 1)'); rho = T(o, 2)';
    assert(all(diff(hu) > 0) && all(diff(rho) > 0), 'prestus:pct:densityCalibration', ...
        'Density calibration %s: HU and density must both increase strictly.', source);
end
