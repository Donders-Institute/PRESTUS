classdef CWSourceMatrix
% CWSourceMatrix  Lazy rank-2 factored CW source signal matrix for k-Wave.
%
%   For a continuous-wave (CW) source every voxel's time signal is:
%
%       source.p(i, t) = Re(C_i)*cos(omega*t) - Im(C_i)*sin(omega*t)
%
%   where C_i = W_src(i,:) * C_elem is the complex amplitude at voxel i.
%   This is a rank-2 matrix: 2*num_src + 2*Nt values suffice rather than
%   the full num_src*Nt dense matrix.
%
%   CWSourceMatrix stores Re_C, Im_C (num_src x 1 single) and computes
%   columns on demand in subsref, so k-Wave's time loop accesses
%   source.p(:, t) without ever materialising the full matrix.
%
%   single() and double() are overloaded to return `this` unchanged, so
%   k-Wave's kspaceFirstOrder_dataCast CPU-cast step passes through
%   harmlessly. gpuArray() is overloaded to upload Re_C/Im_C to the GPU,
%   after which subsref returns gpuArray slices suitable for k-Wave's GPU
%   time loop. isnumeric() returns true so k-Wave's input validation passes.
%
%   Valid for matlab_cpu (DataCast='single'/'double') and matlab_gpu
%   (DataCast='gpuArray-single'/'gpuArray-double') backends.
%   NOT suitable for cpp backends — those write source.p to HDF5 and
%   require a dense matrix; use the full distributed_source_signal there.
%
% Construction:
%   obj = CWSourceMatrix(Re_C, Im_C, omega, dt, Nt)
%
%   Re_C   - num_src x 1  real part of per-voxel complex amplitude [Pa]
%   Im_C   - num_src x 1  imaginary part                           [Pa]
%   omega  - scalar       angular frequency [rad/s]
%   dt     - scalar       time step [s]
%   Nt     - scalar       number of time steps

    properties (Access = private)
        Re_C    % num_src x 1 single
        Im_C    % num_src x 1 single
        omega   % scalar double  [rad/s]
        dt      % scalar double  [s]
        n_src   % scalar int32
        Nt_     % scalar int32  (private; public Nt clashes with kWave vars)
    end

    methods

        function obj = CWSourceMatrix(Re_C, Im_C, omega, dt, Nt)
            obj.Re_C  = single(Re_C(:));
            obj.Im_C  = single(Im_C(:));
            obj.omega = double(omega);
            obj.dt    = double(dt);
            obj.n_src = int32(numel(Re_C));
            obj.Nt_   = int32(Nt);
        end

        % ----------------------------------------------------------------
        % Indexing: support source.p(:, t_index) used by k-Wave time loop
        % and source.p(row, :) used by k-Wave's input checking to infer Nt
        % via flags.source_p = length(source.p(1, :)).
        % Dot-access (properties, methods) is routed to builtin.
        % ----------------------------------------------------------------
        function varargout = subsref(obj, s)
            if strcmp(s(1).type, '.')
                [varargout{1:nargout}] = builtin('subsref', obj, s);
            elseif strcmp(s(1).type, '()') && numel(s) == 1 && numel(s.subs) == 2
                row_sub = s.subs{1};
                col_sub = s.subs{2};
                if ischar(col_sub) || isstring(col_sub)
                    % source.p(rows, ':') — return full time series for
                    % selected rows. k-Wave calls source.p(1,:) in input
                    % checking to determine Nt via length(); returning the
                    % correct 1×Nt vector ensures flags.source_p = Nt so
                    % the source is applied at every time step.
                    t_vec = (0 : double(obj.Nt_)-1) * obj.dt;   % 1 × Nt
                    cs    = cos(obj.omega * t_vec);
                    sn    = sin(obj.omega * t_vec);
                    if isequal(row_sub, ':')
                        varargout{1} = obj.Re_C * cs - obj.Im_C * sn;   % num_src × Nt
                    else
                        varargout{1} = obj.Re_C(row_sub) * cs - obj.Im_C(row_sub) * sn;
                    end
                else
                    % source.p(:, t_index) or source.p(rows, t_index) —
                    % evaluate the CW signal at a specific time index.
                    % t = 0 at index 1, matching k-Wave's t_array = (0:Nt-1)*dt
                    col_sub = double(col_sub);
                    t  = (col_sub - 1) * obj.dt;
                    cs = cos(obj.omega * t);
                    sn = sin(obj.omega * t);
                    if isequal(row_sub, ':')
                        varargout{1} = obj.Re_C * cs - obj.Im_C * sn;
                    else
                        varargout{1} = obj.Re_C(row_sub) * cs - obj.Im_C(row_sub) * sn;
                    end
                end
            else
                error('CWSourceMatrix:unsupportedIndex', ...
                    'CWSourceMatrix only supports (rows, col) indexing; got %s.', ...
                    strjoin(cellfun(@(s_) s_.type, num2cell(s), 'UniformOutput', false)));
            end
        end

        % ----------------------------------------------------------------
        % Shape queries
        % ----------------------------------------------------------------
        function varargout = size(obj, varargin)
            sz = [double(obj.n_src), double(obj.Nt_)];
            if isempty(varargin)
                if nargout <= 1
                    varargout{1} = sz;
                else
                    varargout{1} = sz(1);
                    varargout{2} = sz(2);
                end
            else
                varargout{1} = sz(varargin{1});
            end
        end

        function n = numel(obj)
            n = double(obj.n_src) * double(obj.Nt_);
        end

        function n = length(obj)
            n = max(double(obj.n_src), double(obj.Nt_));
        end

        function e = end(obj, k, ~)
            sz = [double(obj.n_src), double(obj.Nt_)];
            e  = sz(k);
        end

        % ----------------------------------------------------------------
        % Type predicates (k-Wave validation uses these)
        % ----------------------------------------------------------------
        function tf = isnumeric(~),  tf = true;  end
        function tf = isreal(~),     tf = true;  end
        function tf = isfloat(~),    tf = true;  end
        function tf = isempty(obj),  tf = obj.n_src == 0 || obj.Nt_ == 0; end

        % ----------------------------------------------------------------
        % Cast pass-throughs: kspaceFirstOrder_dataCast calls
        %   source.p = single(source.p)  or  source.p = double(source.p)
        % We return `this` unchanged so the lazy object survives dataCast.
        % ----------------------------------------------------------------
        function obj = single(obj), end
        function obj = double(obj), end

        % ----------------------------------------------------------------
        % GPU support: kspaceFirstOrder_dataCast calls
        %   source.p = gpuArray(single(source.p))
        % for matlab_gpu. Uploading Re_C/Im_C to GPU means subsequent
        % subsref calls return gpuArray slices that k-Wave's GPU time
        % loop can use directly without a CPU–GPU round-trip per step.
        % ----------------------------------------------------------------
        function obj = gpuArray(obj)
            obj.Re_C = gpuArray(obj.Re_C);
            obj.Im_C = gpuArray(obj.Im_C);
        end

        function obj = gather(obj)
            obj.Re_C = gather(obj.Re_C);
            obj.Im_C = gather(obj.Im_C);
        end

        % ----------------------------------------------------------------
        % Materialise (for debugging or cpp fallback)
        % ----------------------------------------------------------------
        function M = full(obj)
            t_vec = (0 : double(obj.Nt_)-1) * obj.dt;          % 1 x Nt
            Re = gather(obj.Re_C);
            Im = gather(obj.Im_C);
            M = Re * cos(obj.omega * t_vec) ...
              - Im * sin(obj.omega * t_vec);              % num_src x Nt
        end

    end
end
