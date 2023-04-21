function [error, ax1, ax2, h] = phase_optimization_annulus_full_curve(phase, parameters, velocity, axial_position, desired_intensity_curve, plot_results, opt_limits, weights)
  arguments
      phase
      parameters
      velocity
      axial_position % 0 is the transducer surface
      desired_intensity_curve
      plot_results = 0
      opt_limits = [1, max(axial_position)]
      weights = 0
  end
  
  limit_ind = (axial_position>=opt_limits(1)&axial_position<=opt_limits(2));
  
  p_axial_oneil = focusedAnnulusONeil(parameters.transducer.curv_radius_mm/1e3, ...
    [parameters.transducer.Elements_ID_mm; parameters.transducer.Elements_OD_mm]/1e3, repmat(velocity,1,length(phase)+1), ...
    [0 phase], parameters.transducer.source_freq_hz, parameters.medium.water.sound_speed, ...
    parameters.medium.water.density, (axial_position-0.5)*1e-3);
  
  i_axial_oneil = p_axial_oneil.^2/(2*parameters.medium.water.sound_speed*parameters.medium.water.density) .* 1e-4;
  [~, max_pos] = max(desired_intensity_curve);
  %actual_focal_dist_mm = axial_position(p_axial_oneil==max(p_axial_oneil));
  if weights == 0
      [flhm_center, flhm_center_index] = get_flhm_center_position(axial_position, desired_intensity_curve);
      weights = normpdf(axial_position, axial_position(flhm_center_index)+0.5, axial_position(flhm_center_index)/3);
  end
  
  weights = weights/sum(weights);
  error_v = (i_axial_oneil - desired_intensity_curve).^2.*weights;
  error_v = error_v(limit_ind);
  error = mean(error_v);

  if plot_results
      h = figure;
      ax1 = subplot(1,2,1);
      hold on
      plot(axial_position,  i_axial_oneil )
      plot(axial_position, desired_intensity_curve)
      yyaxis right
      plot(axial_position, weights);
      hold off
      legend(["fitted profile","real profile","cost function","error"])
      ax2 = subplot(1,2,2);
      plot(axial_position(limit_ind), error_v,'-o');
      legend(["error"])
  end
end