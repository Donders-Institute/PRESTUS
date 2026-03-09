function [alpha_coeff, alpha_power] = medium_pct_attenuation(parameters, medium, alpha_coeff, alpha_power, pseudoCT, skull_idx, algorithm)

switch algorithm
    case 'k-plan'

        kPlan_alpha = 13.3; % https://dispatch.k-plan.io/static/docs/simulation-pipeline.html
        kPlan_alpha_power = 1;
        % Note that we allow different values to be specified in the config.
        % If replication of k-Wave is the goal, the above values should be specified.
        % Throw a warning in the case of deviations.
        if medium.skull.alpha_coeff ~= kPlan_alpha || ...
            medium.skull.alpha_power ~= kPlan_alpha_power
            warning('Specified attenuation varies from k-Plan setup.')
        end
        alpha_coeff(skull_idx) = medium.skull.alpha_coeff;
        alpha_power(skull_idx) = medium.skull.alpha_power;

    case 'mueller'

        alpha_min     = 4;        % cortical bone at 500 kHz [dB/cm] [Aubry et al., 2022] 
        alpha_max     = 8.7;      % bone at 500 kHz [dB/cm] [Fry 1978]

        % Finds maximum and minimum values
        HU_min = min(pseudoCT(skull_idx));
        HU_max = max(pseudoCT(skull_idx));

        % estimate attenuation based on (pseudo-)HU
        alpha_pseudoCT(skull_idx) = alpha_min + (alpha_max - alpha_min) * ...
            (1 - (pseudoCT(skull_idx) - HU_min) / (HU_max - HU_min)).^0.5;
        alpha_power(skull_idx) = medium.skull.alpha_power;
        % convert alpha at 500 kHz into prefactor alpha0 (dB/MHz/cm) according to specified alpha_power
        % (definition of lower and upper attenuation bounds is derived from 500kHz)
        alpha_coeff(skull_idx) = alpha_pseudoCT(skull_idx)./(0.5^medium.skull.alpha_power);

    case 'aubry'

        alpha_min     = 0.2;    % Aubry et al., 2003
        alpha_max     = 8;      % Aubry et al., 2003

        phi(skull_idx) = 1-(pseudoCT(skull_idx)/max(pseudoCT(skull_idx)));
        alpha_coeff(skull_idx) = alpha_min + (alpha_max - alpha_min) * (phi(skull_idx).^0.5);
        % regularize sound speed to a minimum of water
        alpha_power(skull_idx) = medium.skull.alpha_power;

    case 'none'

        alpha_coeff(skull_idx) = medium.skull.alpha_coeff;
        alpha_power(skull_idx) = medium.skull.alpha_power;

    otherwise
        error("Specified pCT attenuation mapping is not supported.")
end