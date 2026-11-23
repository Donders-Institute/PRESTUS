function rid = generate_run_id()
% GENERATE_RUN_ID  A random UUID v4 string for telemetry (run ID and install ID).
%
% Drawn from java.util.UUID (seeded by the operating system), not from
% MATLAB's rand: every MATLAB session starts rand at the same default seed,
% so the first draw of each new session (every `matlab -batch` job) would be
% the same ID, and drawing would also move the user's own random stream.
% Without Java, a private stream seeded from the clock is used instead.
%
% Use as:
%   rid = generate_run_id()
%
% Output:
%   rid - char, e.g. '3f2b9c1e-8d4a-4b6f-9e21-7c5d0a1b2c3d'
%
% See also: TRACK_USAGE

    if usejava('jvm')
        rid = char(java.util.UUID.randomUUID().toString());
        return
    end
    s = RandStream('mt19937ar', 'Seed', 'shuffle');        % never the global stream
    bytes = uint8(floor(rand(s, 1, 16) * 256));
    bytes(7) = bitor(bitand(bytes(7), uint8(15)), uint8(64));    % version 4
    bytes(9) = bitor(bitand(bytes(9), uint8(63)), uint8(128));   % variant
    rid = sprintf('%02x%02x%02x%02x-%02x%02x-%02x%02x-%02x%02x-%02x%02x%02x%02x%02x%02x', bytes);
end
