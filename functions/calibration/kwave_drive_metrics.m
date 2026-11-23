function m = kwave_drive_metrics(z_meas, I_meas, z_sim, I_sim, isppa_ref, options)
% KWAVE_DRIVE_METRICS  Absolute agreement of a simulated free-water profile with its target
%
% Scores a simulated on-axis intensity profile at the library's analytical
% drive against what the drive is meant to produce: the hydrophone profile
% scaled to the reference ISPPA (the library stores drives for isppa_ref;
% the measured shape is the target). Unlike the peak-normalised R² of the
% fit, these metrics keep the absolute level, so they show how far the
% analytical drive misses the target ISPPA, and how much of the remaining
% mismatch is shape once the drive is corrected
% (drive × sqrt(isppa_ref / isppa_sim), i.e. intensity × isppa_ref / isppa_sim).
%
% Use as:
%   m = kwave_drive_metrics(z_meas, I_meas, z_sim, I_sim, isppa_ref)
%   m = kwave_drive_metrics(..., options)
%
% Input:
%   z_meas, I_meas - hydrophone profile: [N x 1] mm from the exit plane, intensity (any scale)
%   z_sim, I_sim   - simulated profile: [M x 1] mm from the exit plane, W/cm²
%   isppa_ref      - reference ISPPA the drive is meant to give [W/cm²]
%   options        - struct, optional:
%     .NearFieldSkipMm - ignore z < this [mm] (default 10)
%     .I_analytical    - [N x 1] analytical (Rayleigh) profile at the same drive on z_meas, W/cm²
%
% Output (m):
%   .isppa_sim_wcm2      - simulated peak beyond the near field
%   .isppa_ratio         - isppa_sim / isppa_ref (1 = drive exactly right)
%   .isppa_err_pct       - 100 * (isppa_ratio - 1)
%   .i_at_target_ratio   - simulated intensity at the measured FLHM centre / isppa_ref
%   .nrmse_before        - rms(I_sim - I_target) / isppa_ref, analytical drive
%   .nrmse_after         - same after the drive correction
%   .r2_abs_before, .r2_abs_after - R² of the absolute profiles (not peak-normalised)
%   .isppa_analytical_wcm2, .kwave_over_analytical - with options.I_analytical
%
% See also: KWAVE_CALIBRATE_DRIVE, COMPUTE_FLHM_CENTER

arguments
    z_meas {mustBeNumeric}
    I_meas {mustBeNumeric}
    z_sim {mustBeNumeric}
    I_sim {mustBeNumeric}
    isppa_ref {mustBeNumeric}
    options (1,1) struct = struct()
end

    if nargin < 6; options = struct(); end
    skip = opt(options, 'NearFieldSkipMm', 10);
    z_meas = z_meas(:); I_meas = I_meas(:); z_sim = z_sim(:); I_sim = I_sim(:);

    sel = z_meas >= skip & z_meas >= min(z_sim) & z_meas <= max(z_sim) & isfinite(I_meas);
    zt  = z_meas(sel);
    It  = I_meas(sel) / max(I_meas(sel)) * isppa_ref;          % target: measured shape at isppa_ref
    Is  = interp1(z_sim, I_sim, zt, 'linear');

    m.isppa_sim_wcm2 = max(I_sim(z_sim >= skip));
    m.isppa_ratio    = m.isppa_sim_wcm2 / isppa_ref;
    m.isppa_err_pct  = 100 * (m.isppa_ratio - 1);
    zc = compute_flhm_center(It, zt);
    m.i_at_target_ratio = interp1(z_sim, I_sim, zc, 'linear') / isppa_ref;

    Ia = Is / m.isppa_ratio;                                    % after the drive correction
    m.nrmse_before  = sqrt(mean((Is - It).^2)) / isppa_ref;
    m.nrmse_after   = sqrt(mean((Ia - It).^2)) / isppa_ref;
    ss = sum((It - mean(It)).^2);
    m.r2_abs_before = 1 - sum((Is - It).^2) / ss;
    m.r2_abs_after  = 1 - sum((Ia - It).^2) / ss;

    m.isppa_analytical_wcm2 = NaN; m.kwave_over_analytical = NaN;
    if isfield(options, 'I_analytical') && ~isempty(options.I_analytical)
        Ian = options.I_analytical(:);
        m.isppa_analytical_wcm2 = max(Ian(z_meas >= skip));
        m.kwave_over_analytical = m.isppa_sim_wcm2 / m.isppa_analytical_wcm2;
    end
end

function v = opt(s, f, d)
    if isfield(s, f) && ~isempty(s.(f)); v = s.(f); else; v = d; end
end
