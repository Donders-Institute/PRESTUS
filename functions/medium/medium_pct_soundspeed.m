function [sound_speed] = medium_pct_soundspeed(parameters, sound_speed, density, pseudoCT, skull_idx, algorithm)

switch algorithm
    case 'k-plan'

        sound_speed(skull_idx) = 1.33.*density(skull_idx) + 167;

    case 'marsac'

        c_water       = parameters.medium.water.sound_speed;     % sound speed [m/s]
        c_skull       = 3360;     % max. speed of sound in skull [m/s]
        rho_water     = parameters.medium.water.density;      % density [kg/m^3]
        rho_bone      = 2100;     % max. skull density [kg/m3]

        sound_speed(skull_idx) = c_water + (c_skull - c_water) * ...
            (density(skull_idx) - rho_water) / (rho_bone - rho_water);
        
    case 'aubry'

        c_water = parameters.medium.water.sound_speed;
        c_bone = parameters.medium.skull.sound_speed;

        phi(skull_idx) = 1-(pseudoCT(skull_idx)/max(pseudoCT(skull_idx)));
        sound_speed(skull_idx) = c_water * phi(skull_idx) + ...
            c_bone * (1-phi(skull_idx));
        % regularize sound speed to a minimum of water
        sound_speed(skull_idx) = max(sound_speed(skull_idx),c_water);

    case 'none'

        sound_speed(skull_idx) = parameters.medium.skull.sound_speed;

    otherwise
        error("Specified CT sound speed mapping is not supported.")
end
