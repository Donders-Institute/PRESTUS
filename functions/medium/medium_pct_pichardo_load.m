function [F_sos, F_att] = medium_pct_pichardo_load(parameters)
% MEDIUM_PCT_PICHARDO_LOAD  Load and cache Pichardo 2D interpolant from BabelBrain
%
% Reads MapPichardo.h5 from the BabelBrain submodule and returns two
% griddedInterpolant objects for ρ→SoS and ρ→Att at a given frequency.
% Results are cached in a persistent variable to avoid reloading on every
% call within the same MATLAB session.
%
% The Pichardo map encodes frequency-dependent longitudinal SoS [m/s] and
% attenuation [Np/m] as a function of tissue density [kg/m³] and frequency
% [MHz]. Density and frequency axes are the row/column indices of the 500×500
% grid stored in the H5 file.
%
% Use as:
%   [F_sos, F_att] = medium_pct_pichardo_load(parameters)
%
% Input:
%   parameters - PRESTUS config struct; uses parameters.babelbrain.root (optional)
%                to locate BabelBrain, otherwise falls back to the standard
%                external/BabelBrain submodule path.
%
% Output:
%   F_sos - griddedInterpolant(rho_kg_m3, freq_MHz) → SoS [m/s]
%   F_att - griddedInterpolant(rho_kg_m3, freq_MHz) → Att [Np/m]
%
% See also: MEDIUM_PCT_SOUNDSPEED, MEDIUM_PCT_ATTENUATION

persistent cached_F_sos cached_F_att cached_bb_root

% Resolve BabelBrain root
bb_root = '';
if isfield(parameters, 'babelbrain') && isfield(parameters.babelbrain, 'root')
    bb_root = parameters.babelbrain.root;
end
if isempty(bb_root)
    prestus_root = get_prestus_path();
    bb_root = fullfile(prestus_root, 'external', 'BabelBrain');
end
h5_path = fullfile(bb_root, 'TranscranialModeling', 'MapPichardo.h5');

% Return cached result if the root hasn't changed
if ~isempty(cached_F_sos) && isequal(cached_bb_root, bb_root)
    F_sos = cached_F_sos;
    F_att = cached_F_att;
    return
end

if ~isfile(h5_path)
    error(['medium_pct_pichardo_load: MapPichardo.h5 not found at:\n  %s\n' ...
           'Ensure the BabelBrain submodule is initialised:\n' ...
           '  git submodule update --init external/BabelBrain'], h5_path);
end

rho      = double(h5read(h5_path, '/rho'));    % [kg/m³], length N
freq     = double(h5read(h5_path, '/freq'));   % [MHz],   length M
MapSoS   = double(h5read(h5_path, '/MapSoS')); % [m/s],   N×M
MapAtt   = double(h5read(h5_path, '/MapAtt')); % [Np/m],  N×M

% griddedInterpolant expects monotonically increasing grid vectors.
% BabelBrain stores the map as MapSoS(rho_idx, freq_idx) via RectBivariateSpline
% which takes (rho, freq) as separate 1-D axes — same indexing as griddedInterpolant.
if ~issorted(rho),  [rho,  ri] = sort(rho);  MapSoS = MapSoS(ri,:); MapAtt = MapAtt(ri,:); end
if ~issorted(freq), [freq, fi] = sort(freq); MapSoS = MapSoS(:,fi); MapAtt = MapAtt(:,fi); end

F_sos = griddedInterpolant({rho, freq}, MapSoS, 'linear', 'nearest');
F_att = griddedInterpolant({rho, freq}, MapAtt, 'linear', 'nearest');

cached_F_sos  = F_sos;
cached_F_att  = F_att;
cached_bb_root = bb_root;
