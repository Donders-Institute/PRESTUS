function alpha_db_cm = get_alpha_coeff(medium, freq)
% medium should be brain or skull
% frequency is provided in Hz

if strcmp(medium, 'skull')
    % From Connor et al (2002 Phys. Med. Biol. 47, https://iopscience.iop.org/article/10.1088/0031-9155/47/22/302/pdf):
    % "While a precise consensus on the amplitude attenuation coefficient is not
    % available from the literature, a value of 167 × f/1e6 [Np m] for
    % cortical bone and 300 × f/1e6 [Np m] for trabecular bone is consistent
    % (Fry and Barger 1978, Theismann and Pfander 1949, Martin and McElhaney
    % 1971)".
%     alpha_np_m = freq*(0.5*(167+300))/1e6 % [Np / MHz m]
%     alpha_db_m = 20 * log10(exp(1)) * alpha_np_m; % [db / MHz m]
%     alpha_db_cm = alpha_db_m  / 100; % [db / MHz cm]
%     
    % assuming linear scaling from ITRUSST parameters
    alpha_db_cm = freq/0.5e6*load_parameters().medium.skull.alpha_coef_true;
    
elseif strcmp(medium, 'brain')
    alpha_db_cm = freq/0.5e6*load_parameters().medium.brain.alpha_coef_true;   
end
end