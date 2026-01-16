function alpha_db_cm = get_alpha_coeff(medium, freq)

% GET_ALPHA_COEFF Computes the amplitude attenuation coefficient for a given medium and frequency.
%
% This function calculates the amplitude attenuation coefficient (`alpha_db_cm`) 
% in decibels per centimeter for a specified medium and frequency. The calculation 
% assumes linear scaling based on tissue-specific parameters from ITRUSST benchmarks.
%
% Input:
%   medium - String specifying the tissue type (e.g., 'brain', 'skull').
%   freq   - Scalar specifying the frequency in Hz.
%
% Output:
%   alpha_db_cm - Amplitude attenuation coefficient in decibels per MHz per centimeter.

    % Load tissue-specific parameters for the specified medium
    medium = load_parameters().medium.(medium);

    % Compute the attenuation coefficient using tissue-specific parameters
    % `alpha_0_true` and `alpha_power_true` define the scaling relationship
    alpha_db_cm = medium.alpha_0_true * ((freq / 1e6) .^ medium.alpha_power_true); % [db / MHz cm]

end