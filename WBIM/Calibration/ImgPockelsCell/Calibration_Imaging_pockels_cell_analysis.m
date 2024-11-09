DataManager = WBIMFileManager;
calibration_folder = 'G:\My Drive\Microscope\Calibration\ImagingLaserPockelsCell\20210802';
data_fp = fullfile(calibration_folder, 'Pockels_cell_calibration.csv');
beam_power = readtable(data_fp);
% beam_power = beam_power(1:end-1, :);
%%
result = struct;
result.data_x = beam_power.Voltage_V_;
result.data_y = beam_power.Power_mW_;
result.fit_fun = @(c, v) c(1) .* cosd(c(2) .* v + c(3)).^2;
result.fit_coeff = lsqcurvefit(result.fit_fun, [max(data_y), 0.25, 90], result.data_x, result.data_y);
result.fit_y = result.fit_fun(fit_coeff, result.data_x);
result.fit_rho = corr(result.data_y, result.fit_y);
result.data_filepath = data_fp;
result.filepath = fullfile(DataManager.SCRIPT_PATH, 'Calibration/ImgPockelsCell', ...
    'ImgPockelPwrVsVol.mat');
save(result.filepath, '-struct', 'result');
%%
fig_hdl = figure;
ax_hdl = axes(fig_hdl);
scatter(ax_hdl, data_x, data_y, 'ro');
hold(ax_hdl, 'on');


plot_fit_hdl = plot(ax_hdl, result.data_x, result.fit_y, 'k', 'LineWidth', 1);

ax_hdl.XLabel.String = 'Voltage (V)';
ax_hdl.YLabel.String = 'Power (mW)';
legend(ax_hdl, plot_fit_hdl, sprintf('\\rho = %.3f', result.fit_rho), 'Location', 'best');
ax_hdl.Title.String = sprintf('%.3e Cos^2[%.3e V + %.3e]', result.fit_coeff);

fig_fp = fullfile(calibration_folder, sprintf('Imaging_pockels_cell_pwr_vs_vol_20210802.png'));
fun_print_image_in_several_formats(fig_hdl, fig_fp);


% result = WBIMCalibration.FitPowerVsAngle(beam_power.Voltage_V_, beam_power.Power_mW_);
% result.data_filepath = data_fp;
% result.filepath = fullfile(DataManager.SCRIPT_PATH, 'Calibration/ImgPockelsCell', ...
%     'ImgPockelsCell_20210802.mat');
% save(result.filepath, '-struct', 'result_s');
% 
% [fig_s, ax_s] = WBIMCalibration.PlotPowerVsAngleFit(result);
% ax_s.YLabel.String = 'Power (mW)';
% ax_s.XLabel.String = 'Voltage (V)';
