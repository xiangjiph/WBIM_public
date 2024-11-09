DataManager = WBIMFileManager;
exp_folder = 'G:\My Drive\Microscope\Calibration\ImagingLaserPowerVsAngle\20211116';
measurement_name = 'ImgPowVsAng_20211116_2';
exp_data_fp = fullfile(exp_folder, sprintf('%s.csv', measurement_name));
exp_data = readtable(exp_data_fp);
%%
result = struct;
result.data_filepath = exp_data_fp;
result.data_x = exp_data.Angle_deg;
result.data_y = exp_data.Power_mW;
result.fit_fun = @(c, theta) c(1) .* cosd(2 * (theta + c(2))).^2;

result.fit_coeff = lsqcurvefit(result.fit_fun, [150, 0], result.data_x, result.data_y);
result.fit_y = result.fit_fun(result.fit_coeff, result.data_x);
result.fit_r = corr(result.data_y, result.fit_y);
result.filepath = fullfile(DataManager.SCRIPT_PATH, 'Calibration', 'ImgWPPwrVsAgl', ...
    sprintf('%s.mat', measurement_name));
save(result.filepath, '-struct', 'result');
xml_write(strrep(result.filepath, '.mat', '.xml'), result);
%% Visualization
fig_hdl = figure;
ax_hdl = axes(fig_hdl);
scatter(ax_hdl, result.data_x, result.data_y, 'ro');
% Fit with cos
hold(ax_hdl, 'on');
plot_fit_hdl = plot(ax_hdl, result.data_x, result.fit_y, 'k', 'LineWidth', 1);
ax_hdl.XLabel.String = 'Angle (degrees)';
ax_hdl.YLabel.String = 'Power (mW)';
if result.fit_coeff(2) > 0
    ax_hdl.Title.String = sprintf('%.1f Cos^2[2(\\theta + %.1f)]', result.fit_coeff);
else
    ax_hdl.Title.String = sprintf('%.1f Cos^2[2(\\theta - %.1f)]', abs(result.fit_coeff));
end
fig_fp = fullfile(exp_folder, sprintf('%s_fit.png', measurement_name));
% fun_print_image_in_several_formats(fig_hdl, fig_fp);
%%
