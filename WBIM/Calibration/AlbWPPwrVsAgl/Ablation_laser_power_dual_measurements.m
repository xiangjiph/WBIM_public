date_str = '20230526';
DataManager = WBIMFileManager;
calibration_folder = fullfile('D:\My Drive\Microscope\Calibration\AblationLaserPowerVsAngle', ...
    date_str);
data_fp = fullfile(calibration_folder, sprintf('AblationPowerAngleCalibration_%s.xlsx', ...
    date_str));
beam_power = readtable(data_fp);
% beam_power = beam_power(1:end-1, :);
%%
result_s = WBIMCalibration.FitPowerVsAngle(beam_power.Angle, ...
    beam_power.Ablation_path_W_);
result_s.data_filepath = data_fp;
result_s.filepath = fullfile(DataManager.SCRIPT_PATH, 'Calibration/AlbWPPwrVsAgl', ...
    sprintf('AlbPowVsAng_%s.mat', date_str));
save(result_s.filepath, '-struct', 'result_s');
%%
[fig_s, ax_s] = WBIMCalibration.PlotPowerVsAngleFit(result_s);
ax_s.YLabel.String = 'S-beam power (W)';
ax_s.XLabel.String = 'Angle \theta (^\circ)';
fig_fn = fullfile(calibration_folder, 'Ablation_path_power_vs_angle_fitting.png');
fun_print_image_in_several_formats(fig_s, fig_fn);
%%
result_p = WBIMCalibration.FitPowerVsAngle(beam_power.Angle, beam_power.Othogonal_path_mW_);
exp_p_power = result_s.fit_coeff(1) * (0.5 + 0.5 * cosd(4 * (result_p.angle + result_s.fit_coeff(2) + 45)));
figure;plot(result_p.angle, exp_p_power);

corr_coeff = result_p.power ./ exp_p_power;

result_p.data_filepath = data_fp;
result_p.filepath = fullfile(DataManager.SCRIPT_PATH, 'Calibration/AlbWPPwrVsAgl', ...
    sprintf('AlbOthPowVsAng_%s.mat', date_str));
result_p.angle2power = griddedInterpolant(result_p.angle, result_p.power, ...
    'linear', 'none');
save(result_p.filepath, '-struct', 'result_p');

[fig_p, ax_p] = WBIMCalibration.PlotPowerVsAngleFit(result_p);
fig_fn = fullfile(calibration_folder, 'Orthogonal_path_power_vs_angle_fitting.png');
fun_print_image_in_several_formats(fig_p, fig_fn);
