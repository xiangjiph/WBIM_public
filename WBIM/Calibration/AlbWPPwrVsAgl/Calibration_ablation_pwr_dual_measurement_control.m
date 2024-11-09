%% Connect to power meter
pm102 = WBIMPowerMeter(WBIMConfig.ABLATION_THERMAL_PM_NAME);
pm100d = WBIMPowerMeter(WBIMConfig.ABLATION_PHOTODIODE_PM_NAME);
wavelength_nm = 800;
pm102.wavelength_nm = wavelength_nm;
pm100d.wavelength_nm = wavelength_nm;

pm100d.meter.setPowerAutoRange(1);
%% Measurement Parameters
date_str = datestr(now, 'YYYYmmDD');
data = struct;
data.hwp_angle_list_deg = [-25 : 1 : 40];
data.num_angle = numel(data.hwp_angle_list_deg);
data.equilibrium_wait_time_s = 5;
data.acq_avg_time_s = 0.1;
data.num_measurement = 100;
[data.power_W, data.oth_power_W] = deal(zeros(data.num_measurement, data.num_angle));
calibration_folder = fullfile(DataManager.SCRIPT_PATH, 'Calibration', 'AlbWPPwrVsAgl');
data.filepath = fullfile(calibration_folder, ...
    sprintf('Calibration_AlbWPPwrVsAgl_dual_%s.mat', date_str));
%% Iterate through all wavelengths
pm102.avg_time_s = data.acq_avg_time_s;
pm100d.avg_time_s = data.acq_avg_time_s;

t_tic = tic;
wbim.ablation.enable_laser_output();
for i = 1 : data.num_angle
    tmp_ang = data.hwp_angle_list_deg(i);
    wbim.ablation.h_hwp.moveTo(tmp_ang);
    pause(data.equilibrium_wait_time_s);
    % Do multiple measurements
    data.power_W(:, i) = pm102.measure_power(data.num_measurement, 0);
    data.oth_power_W(:, i) = pm100d.measure_power(data.num_measurement, 0);
    fprintf('Fiish measuring power for angle %.2f. Elapsed time is %.2f seconds. \n', ...
        tmp_ang, toc(t_tic));
end
wbim.ablation.disable_laser_output();
fprintf('Finish measurements. Total elapsed time is %.2f seconds.\n', toc(t_tic));
% Disconnect
pm102.delete();
pm100d.delete();
%%
save(data.filepath, '-struct', 'data');
%% 
data.power_W_mean = mean(data.power_W, 1).';
data.oth_power_W_mean = mean(data.oth_power_W, 1).';
%%
result_s = WBIMCalibration.FitPowerVsAngle(data.hwp_angle_list_deg.', ...
    data.power_W_mean);
result_s.data_filepath = data.filepath;
result_s.filepath = fullfile(calibration_folder, sprintf('AlbPowVsAng_%s.mat', date_str));
save(result_s.filepath, '-struct', 'result_s');
%% 
result_p = WBIMCalibration.FitPowerVsAngle(data.hwp_angle_list_deg, data.oth_power_W_mean);
exp_p_power = result_s.fit_coeff(1) * (0.5 + 0.5 * cosd(4 * (result_p.angle + result_s.fit_coeff(2) + 45)));
figure;plot(result_p.angle, exp_p_power);

corr_coeff = result_p.power ./ exp_p_power;

result_p.data_filepath = data.filepath;
result_p.filepath = fullfile(calibration_folder, sprintf('AlbOthPowVsAng_%s.mat', date_str));
result_p.angle2power = griddedInterpolant(result_p.angle, result_p.power, ...
    'linear', 'none');
save(result_p.filepath, '-struct', 'result_p');

[fig_p, ax_p] = WBIMCalibration.PlotPowerVsAngleFit(result_p);
% fig_fn = fullfile(calibration_folder, 'Orthogonal_path_power_vs_angle_fitting.png');
% fun_print_image_in_several_formats(fig_p, fig_fn);