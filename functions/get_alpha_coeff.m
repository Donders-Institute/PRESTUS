function alpha_db_cm = get_alpha_coeff(medium, freq)
% medium should be brain or skull
% frequency is provided in Hz
%     From Connor et al (2002 Phys. Med. Biol. 47, https://iopscience.iop.org/article/10.1088/0031-9155/47/22/302/pdf):
%     "While a precise consensus on the amplitude attenuation coefficient is not
%     available from the literature, a value of 167 xf/1e6 [Np m] for
%     cortical bone and 300 xf/1e6 [Np m] for trabecular bone is consistent
%     (Fry and Barger 1978, Theismann and Pfander 1949, Martin and McElhaney
%     1971)".
%     alpha_np_m = freq*(0.5*(167+300))/1e6 % [Np / MHz m]
%     alpha_db_m = 20 * log10(exp(1)) * alpha_np_m; % [db / MHz m]
%     alpha_db_cm = alpha_db_m  / 100 % [db / MHz cm]
%     
%     assuming linear scaling from ITRUSST parameters
    medium = load_parameters().medium.(medium);
    alpha_db_cm = medium.alpha_0_true*((freq/1e6).^medium.alpha_power_true);   
end
