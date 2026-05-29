function [trans_pos_ras, focus_pos_ras] = launch_web_alignment(parameters, nii_path, trans_pos_ras, focus_pos_ras, options)

% LAUNCH_WEB_ALIGNMENT  Interactive browser-based transducer placement editor
%
% Exports the NIfTI volume and current transducer/focus positions to a
% temporary working directory, starts a local Python HTTP server, and opens
% the NiiVue-based viewer in the default browser. Blocks until the user
% clicks "Save & close" in the viewer, then reads back the updated RAS+
% positions and returns them.
%
% The viewer allows the user to:
%   - Navigate the volume in three orthogonal views
%   - Drag the transducer (red) and focus (blue) markers to new positions
%   - Inspect the beam axis in real time
%   - Save the result or cancel to keep the original positions
%
% Use as:
%   [trans_pos_ras, focus_pos_ras] = ...
%       launch_web_alignment(parameters, nii_path, trans_pos_ras, focus_pos_ras)
%   [trans_pos_ras, focus_pos_ras] = ...
%       launch_web_alignment(..., Name=Value)
%
% Input:
%   parameters    - (1,1) simulation parameters struct
%   nii_path      - char | string path to the NIfTI volume to display
%   trans_pos_ras - [1x3] initial transducer position in RAS+ mm
%   focus_pos_ras - [1x3] initial focus position in RAS+ mm
%
% Name-Value options:
%   Port     - (1,1) integer HTTP port (default: 8742)
%   Timeout  - (1,1) seconds to wait before giving up (default: 600)
%   Browser  - char command to open browser ('auto' uses system default)
%
% Output:
%   trans_pos_ras - [1x3] updated transducer position in RAS+ mm
%                  (unchanged if user cancels or times out)
%   focus_pos_ras - [1x3] updated focus position in RAS+ mm
%
% Requirements:
%   Python 3.x must be on PATH (used for the embedded HTTP server).
%
% See also: POSITION_TRANSDUCER_LOCALITE, TRANSFORM_COORDINATES

    arguments
        parameters    (1,1) struct
        nii_path      (1,:) char
        trans_pos_ras (1,3) double
        focus_pos_ras (1,3) double
        options.Port    (1,1) double  = 8742
        options.Timeout (1,1) double  = 600
        options.Browser (1,:) char    = 'auto'
    end

    port    = options.Port;
    timeout = options.Timeout;

    %% Resolve paths
    this_dir   = fileparts(mfilename('fullpath'));
    viewer_dir = fullfile(this_dir, 'viewer');
    server_py  = fullfile(this_dir, 'server.py');
    work_dir   = fullfile(tempdir, sprintf('prestus_web_alignment_%d', port));
    if ~exist(work_dir, 'dir'); mkdir(work_dir); end

    result_file = fullfile(work_dir, 'placement_result.json');
    if exist(result_file, 'file'); delete(result_file); end

    %% Build transducer geometry summary for the viewer
    tr = parameters.transducer(1);
    tr_type = tr.type;
    curv_radius_mm = tr.(tr_type).curv_radius_mm;
    if strcmp(tr_type, 'annular')
        diameter_mm = tr.annular.elem_od_mm(end);
    else
        diameter_mm = tr.(tr_type).outer_diameter_mm;
    end

    %% Write initial placement JSON
    init_json = struct();
    init_json.trans_pos_ras = trans_pos_ras;
    init_json.focus_pos_ras = focus_pos_ras;
    depth_mm = tr.(tr_type).depth_mm;
    if isempty(depth_mm); depth_mm = 0; end

    init_json.transducer = struct( ...
        'curv_radius_mm', curv_radius_mm, ...
        'diameter_mm',    diameter_mm, ...
        'depth_mm',       depth_mm);
    % Preserve the original extension so NiiVue receives the correct file type
    % (copying a plain .nii as .nii.gz causes NiiVue to fail gzip decompression)
    if endsWith(nii_path, '.nii.gz')
        nii_dest_name = 'volume.nii.gz';
    else
        nii_dest_name = 'volume.nii';
    end
    init_json.nii_filename = nii_dest_name;

    init_file = fullfile(work_dir, 'placement_init.json');
    fid = fopen(init_file, 'w');
    fprintf(fid, '%s', jsonencode(init_json, 'PrettyPrint', true));
    fclose(fid);

    %% Copy NIfTI to work dir (server only serves from work_dir)
    dest_nii = fullfile(work_dir, nii_dest_name);
    if ~strcmp(nii_path, dest_nii)
        copyfile(nii_path, dest_nii);
    end

    %% Copy viewer static files to work_dir so the server can find them
    viewer_files = {'index.html', 'app.js', 'style.css'};
    for i = 1:numel(viewer_files)
        src = fullfile(viewer_dir, viewer_files{i});
        dst = fullfile(work_dir,   viewer_files{i});
        if ~exist(dst, 'file') || file_newer_than(src, dst)
            copyfile(src, dst);
        end
    end

    %% Start Python HTTP server
    cmd_server = sprintf('python3 "%s" --port %d --workdir "%s" &', ...
        server_py, port, work_dir);
    [status, ~] = system(cmd_server);
    if status ~= 0
        error('PRESTUS:webAlignment:serverFailed', ...
            'Failed to start web alignment server. Ensure python3 is on PATH.');
    end

    % Give the server a moment to bind
    pause(1.0);

    %% Open browser
    url = sprintf('http://localhost:%d/index.html', port);
    if strcmp(options.Browser, 'auto')
        if ismac
            system(sprintf('open "%s" &', url));
        elseif isunix
            system(sprintf('xdg-open "%s" &', url));
        else
            system(sprintf('start "" "%s"', url));
        end
    else
        system(sprintf('%s "%s" &', options.Browser, url));
    end

    fprintf('Web alignment viewer opened at %s\n', url);
    fprintf('Adjust transducer and focus positions, then click "Save & close".\n');
    fprintf('Waiting up to %d seconds ...\n', timeout);

    %% Poll for result
    t_start = tic;
    result_received = false;
    while toc(t_start) < timeout
        if exist(result_file, 'file')
            result_received = true;
            break;
        end
        pause(1.0);
    end

    %% Kill the server
    kill_cmd = sprintf('lsof -ti tcp:%d | xargs kill -9 2>/dev/null; true', port);
    system(kill_cmd);

    %% Parse result
    if ~result_received
        warn('Web alignment timed out after %d s — original positions retained.', timeout);
        return;
    end

    fid = fopen(result_file, 'r');
    raw = fread(fid, '*char')';
    fclose(fid);

    result = jsondecode(raw);

    if isfield(result, 'cancelled') && result.cancelled
        fprintf('Web alignment cancelled — original positions retained.\n');
        return;
    end

    trans_pos_ras = result.trans_pos_ras(:)';
    focus_pos_ras = result.focus_pos_ras(:)';

    fprintf('Web alignment saved:\n');
    fprintf('  trans_pos_ras = [%.2f  %.2f  %.2f] mm\n', trans_pos_ras);
    fprintf('  focus_pos_ras = [%.2f  %.2f  %.2f] mm\n', focus_pos_ras);
end

% ---------------------------------------------------------------------------
function newer = file_newer_than(a, b)
    da = dir(a);
    db = dir(b);
    newer = ~isempty(da) && (isempty(db) || da.datenum > db.datenum);
end
