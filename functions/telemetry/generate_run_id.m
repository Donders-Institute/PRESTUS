function rid = generate_run_id()
% GENERATE_RUN_ID  Generate a random UUID v4 string for telemetry.
    bytes = uint8(floor(rand(1,16) * 256));
    bytes(7) = bitor(bitand(bytes(7), uint8(15)), uint8(64));
    bytes(9) = bitor(bitand(bytes(9), uint8(63)), uint8(128));
    rid = sprintf('%02x%02x%02x%02x-%02x%02x-%02x%02x-%02x%02x-%02x%02x%02x%02x%02x%02x', bytes);
end
