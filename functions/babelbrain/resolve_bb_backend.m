function [backend, device] = resolve_bb_backend(backend, device)
% RESOLVE_BB_BACKEND  Resolve 'auto' BabelBrain computing_backend / device_name
%
% Maps parameters.babelbrain.computing_backend ('auto', a backend name such as
% 'cuda', or 1-4) to a concrete numeric backend (1=CUDA, 2=OpenCL, 3=Metal,
% 4=CPU) and device name.
% Priority: Metal (macOS) > CUDA (nvidia-smi) > CPU. A numeric backend is
% returned unchanged; only an 'auto' device is then filled in.
%
% Use as:
%   [backend, device] = resolve_bb_backend(backend, device)
%
% See also: BABELBRAIN_WRAPPER, EXPORT_BABELBRAIN_INPUTS

    if isstring(backend); backend = char(backend); end
    if isstring(device);  device  = char(device);  end

    % Named backends ('cuda' | 'opencl' | 'metal' | 'cpu') -> BabelBrain integer
    if ischar(backend) && ~strcmpi(backend, 'auto')
        names = {'cuda', 'opencl', 'metal', 'cpu'};
        idx = find(strcmpi(backend, names), 1);
        if isempty(idx)
            error('resolve_bb_backend:unknownBackend', ...
                'Unknown babelbrain.computing_backend ''%s''; use auto, cuda, opencl, metal, cpu or 1-4.', backend);
        end
        backend = idx;
    end

    if ~(ischar(backend) && strcmpi(backend, 'auto'))
        % Already a numeric backend — only resolve device if still 'auto'
        if ischar(device) && strcmp(device, 'auto')
            device = default_device_for_backend(backend);
        end
        return
    end

    if ismac
        backend = 3;   % Metal
        device  = 'Apple';
        fprintf('[resolve_bb_backend] Auto-detected backend: Metal (macOS)\n');
        return
    end

    % Linux/Windows: probe for CUDA
    [status, out] = run_command({'nvidia-smi', '--query-gpu=name', '--format=csv,noheader'});
    if status == 0
        lines = strtrim(strsplit(strtrim(out), newline));
        lines = lines(~cellfun(@isempty, lines));
        if ~isempty(lines)
            backend = 1;   % CUDA
            device  = lines{1};
            fprintf('[resolve_bb_backend] Auto-detected backend: CUDA, device="%s"\n', device);
            return
        end
    end

    % Fallback: CPU — warn so users on GPU machines notice driver/detection issues
    backend = 4;
    device  = 'CPU';
    warning('babelbrain_wrapper:cpuFallback', ...
        '[resolve_bb_backend] No GPU detected (nvidia-smi failed or returned no devices); falling back to CPU. Simulation may be slow.');
end

function device = default_device_for_backend(backend)
% Return a safe default device string for a known numeric backend. 'CPU' is
% the sentinel run_babelbrain_case.py auto-detects from for CUDA/OpenCL; for
% CUDA probe nvidia-smi so standalone RunCases scripts get a real name too.
    switch backend
        case 3;  device = 'Apple';
        case 1
            device = 'CPU';
            [status, out] = run_command({'nvidia-smi', '--query-gpu=name', '--format=csv,noheader'});
            if status == 0
                lines = strtrim(strsplit(strtrim(out), newline));
                lines = lines(~cellfun(@isempty, lines));
                if ~isempty(lines); device = lines{1}; end
            end
        otherwise; device = 'CPU';
    end
end
