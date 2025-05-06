function [sensor_data] = convert_axisymmetric_to_2d(sensor_data)

% Convert the sensor data from symmetric about x=0 to a mirrored 2d setup.
sensor_data.p_final = cat(1, sensor_data.p_final, flipud(sensor_data.p_final));
sensor_data.p_max_all = cat(1, sensor_data.p_max_all, flipud(sensor_data.p_max_all));

end