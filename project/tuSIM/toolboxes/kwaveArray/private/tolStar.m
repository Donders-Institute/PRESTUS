function [lin_ind, varargout] = tolStar(tolerance, kgrid, point)
%TOLSTAR Compute spatial extent of BLI above given tolerance.
%
% DESCRIPTION:
%     tolStar computes a set of subscripts corresponding to locations where
%     the magnitude of a multidimensional sinc function has not yet decayed
%     below a specified tolerance. These subscripts are relative to the 
%     nearest grid node to a given Cartesian point. Note, a sinc
%     approximation is used for the band-limited interpolant (BLI).
%
% USAGE:
%     [lin_ind, x_ind] = tolStar(tolerance, kgrid, point)
%     [lin_ind, x_ind, y_ind] = tolStar(tolerance, kgrid, point)
%     [lin_ind, x_ind, y_ind, z_ind] = tolStar(tolerance, kgrid, point)
%
% INPUTS:
%     tolerance             - Scalar value controlling where the spatial
%                             extent of the BLI at each point is trunctated
%                             as a  portion of the maximum value.
%     kgrid                 - Object of the kWaveGrid class defining the
%                             Cartesian and k-space grid fields.
%     point                 - Cartesian coordinates defining location of
%                             the BLI.
% 
% OUTPUTS:
%     lin_ind               - Linear indices following MATLAB's column-wise
%                             linear matrix index ordering.
%     x_ind, y_ind, z_ind   - Grid indices (y and z only returned in 2D and
%                             3D).
%     
% ABOUT:
%     author                - Elliott Wise
%     date                  - 16th March 2017
%     last update           - 21st May 2018
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2017-2018 Elliott Wise

% This file is part of k-Wave. k-Wave is free software: you can
% redistribute it and/or modify it under the terms of the GNU Lesser
% General Public License as published by the Free Software Foundation,
% either version 3 of the License, or (at your option) any later version.
% 
% k-Wave is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
% more details. 
% 
% You should have received a copy of the GNU Lesser General Public License
% along with k-Wave. If not, see <http://www.gnu.org/licenses/>.

% store a canonical star for use whenever tol remains unchanged
persistent tol subs0

% compute a new canonical star if the tolerance changes, where the star
% gives indices relative to [0 0 0]
if isempty(tol) || (tolerance ~= tol) || (length(subs0) ~= kgrid.dim)
% (no subs yet)    (new tolerance)       (new dimensionality)
    
    % assign tolerance value
    tol = tolerance;
    
    % the on-axis decay of the BLI is given by dx/(pi*x), thus the grid 
    % point at which the BLI has decayed to the tolerance value can be 
    % calculated by 1/(pi*tolerance)
    decay_subs = ceil(1/(pi*tol));
    
    % compute grid indices along axes
    is0 = -decay_subs:decay_subs;          
    
    % replicate grid vectors
    switch kgrid.dim
        case 2
            [is0, js0] = ndgrid(is0, is0);
        case 3
            [is0, js0, ks0] = ndgrid(is0, is0, is0);
    end
    
    % only keep grid indices where the BLI is above the tolerance value
    % (i.e., within the star)
    switch kgrid.dim
        case 1
            subs0 = {is0};
        case 2
            instar = logical(abs(is0 .* js0) <= decay_subs);
            is0 = is0(instar);
            js0 = js0(instar);
            subs0 = {is0, js0};
        case 3
            instar = logical(abs(is0.*js0.*ks0) <= decay_subs);
            is0 = is0(instar);
            js0 = js0(instar);
            ks0 = ks0(instar);
            subs0 = {is0, js0, ks0};
    end
    
end

% get nearest subs to given point and shift canonical subs by this amount
i = round((point(1) + kgrid.x_size/2) * kgrid.Nx/kgrid.x_size + 1);
is = subs0{1} + i;
if kgrid.dim > 1
    j = round((point(2) + kgrid.y_size/2) * kgrid.Ny/kgrid.y_size + 1);
    js = subs0{2} + j;
end
if kgrid.dim > 2
    k = round((point(3) + kgrid.z_size/2) * kgrid.Nz/kgrid.z_size + 1);
    ks = subs0{3} + k;
end

% remove those beyond the bounds of the simulation
inbounds = (is >= 1) & (is <= kgrid.Nx);
if kgrid.dim > 1
    inbounds = inbounds & (js >= 1) & (js <= kgrid.Ny);
end
if kgrid.dim > 2
    inbounds = inbounds & (ks >= 1) & (ks <= kgrid.Nz);
end
is = is(inbounds);
if kgrid.dim > 1
    js = js(inbounds); 
end
if kgrid.dim > 2
    ks = ks(inbounds); 
end

% convert to indices
switch kgrid.dim
    case 1
        lin_ind = is;
    case 2
        lin_ind = sub2ind([kgrid.Nx, kgrid.Ny], is, js);
    case 3
        lin_ind = sub2ind([kgrid.Nx, kgrid.Ny, kgrid.Nz], is, js, ks);
end

% generate variable output
varargout{1} = is;
if (kgrid.dim > 1) && (nargout > 1)
    varargout{2} = js; 
end
if (kgrid.dim > 2) && (nargout > 2)
    varargout{3} = ks; 
end