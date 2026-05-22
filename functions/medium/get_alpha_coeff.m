function alpha_db_cm = get_alpha_coeff(medium, freq)
% GET_ALPHA_COEFF  Compute amplitude attenuation coefficient for a given medium and frequency
%
% Looks up tissue-specific alpha_coeff and alpha_power from the PRESTUS
% parameter set (via load_parameters) and evaluates the power-law model
% alpha = alpha_coeff * (freq/1e6)^alpha_power. Tissue parameters are
% based on ITRUSST benchmark values.
%
% Use as:
%   alpha_db_cm = get_alpha_coeff(medium, freq)
%
% Input:
%   medium      - tissue type name matching a field in load_parameters().medium
%                 (e.g., 'brain', 'skull')
%   freq        - source frequency [Hz]
%
% Output:
%   alpha_db_cm - amplitude attenuation coefficient [dB/(MHz cm)]
%
% See also: FITPOWERLAWPARAMSMULTI, MEDIUM_SETUP

arguments
    medium (1,:) char
    freq   (1,1) {mustBeNumeric, mustBePositive}
end

    % Load tissue-specific parameters for the specified medium
    medium = load_parameters().medium.(medium);

    % Compute the attenuation coefficient using tissue-specific parameters
    % `alpha_coeff` and `alpha_power` define the scaling relationship
    alpha_db_cm = medium.alpha_coeff * ((freq / 1e6) .^ medium.alpha_power); % [db / MHz cm]

end