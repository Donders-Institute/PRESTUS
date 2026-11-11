function np_m = db_cm_mhz_to_np_m(coeff_db, freq_hz, power)
% DB_CM_MHZ_TO_NP_M  Convert an attenuation or absorption coefficient to nepers per metre
%
% Tissue tables give attenuation and absorption as a coefficient in
% dB/(cm MHz^y): the loss at 1 MHz, scaled to other frequencies by the power
% law f^y (y = alpha_power). Heat deposition needs the loss in nepers per
% metre at the source frequency:
%   np_m = coeff_db * (f / 1 MHz)^y * 100 / (20 log10(e))
% (100 converts per cm to per m; 20 log10(e) = 8.686 dB per neper).
%
% Use as:
%   np_m = db_cm_mhz_to_np_m(0.59, 500e3, 1.2)    % brain attenuation at 500 kHz
%   np_m = db_cm_mhz_to_np_m(1, freq_hz, alpha_power)  % conversion factor per voxel
%
% Input:
%   coeff_db - coefficient [dB/(cm MHz^y)]; scalar or array
%   freq_hz  - source frequency [Hz]
%   power    - power-law exponent y [-]; scalar or array the size of coeff_db
%
% Output:
%   np_m - coefficient at freq_hz [Np/m]
%
% See also: MEDIUM_SETUP, NIFTI_MEDIUM_PROPERTIES, DB2NEPER

    np_m = coeff_db .* (freq_hz / 1e6) .^ power * 100 / (20 * log10(exp(1)));
end
