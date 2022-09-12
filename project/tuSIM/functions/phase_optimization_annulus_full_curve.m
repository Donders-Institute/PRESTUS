function error = phase_optimization_annulus(phase, parameters, velocity, axial_position, desired_intensity_curve, plot_results, opt_limits, weights)
  arguments
      phase
      parameters
      velocity
      axial_position
      desired_intensity_curve
      plot_results = 0
      opt_limits = [1, length(axial_position)]
      weights = 0
  end
  
  limit_ind = (axial_position>=opt_limits(1)&axial_position<=opt_limits(2));
  
  p_axial_oneil = focusedAnnulusONeil(parameters.transducer.curv_radius_mm/1e3, ...
    [parameters.transducer.Elements_ID_mm; parameters.transducer.Elements_OD_mm]/1e3, repmat(velocity,1,4), ...
    [0 phase], parameters.transducer.source_freq_hz, parameters.medium.water.sound_speed, ...
    parameters.medium.water.density, (axial_position-0.5)*1e-3);
  
  i_axial_oneil = p_axial_oneil.^2/(2*parameters.medium.water.sound_speed*parameters.medium.water.density) .* 1e-4;
  [~, max_pos] = max(desired_intensity_curve);
  %actual_focal_dist_mm = axial_position(p_axial_oneil==max(p_axial_oneil));
  if weights == 0
      weights = normpdf(axial_position, axial_position(max_pos)+0.5, 10);
  end
  weights = weights/sum(weights);
  if plot_results
  figure
  hold on
  plot(axial_position,  i_axial_oneil )
  plot(axial_position, desired_intensity_curve)
  yyaxis right
  plot(axial_position, weights);
  hold off
  end
  error_v = (i_axial_oneil - desired_intensity_curve).^2.*weights;
  error_v = error_v(limit_ind);
  error = mean(error_v);
end