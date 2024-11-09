rs = dabs.resources.ResourceStore();
pockel_ao = rs.filterByName('/vDAQ0/AO5');
hSI.hShutters.hShutters{1}.open();
hSI.hShutters.hShutters{1}.close();
pockel_ao.setValue(0.05);
%% TODO
% 1. Figure out how control Thorlabs powermeter from MATLAB
% 2. Write a for loop to measure the imaging laser output power versus
% control voltage to the Pockels Cell. 
% 3. Repeat the measurements for 810 nm, 900 nm, 950 nm, 1029 nm. 

%%
DataManager = WBIMFileManager;
calibration_folder = 'D:\My Drive\Microscope\Calibration\ImagingLaserPockelsCell\20220607';
data_fp = fullfile(calibration_folder, 'WBIMCalibration_20220607_PockelCell_2.xlsx');
beam_power = readtable(data_fp);
% beam_power = beam_power(1:end-1, :);
[~, max_ind] = max(beam_power.Fractional_power);
%% Overwrite SI LUT
hSI.hBeams.hBeams{1}.powerFraction2ModulationVoltLut = cat(2, ...
    beam_power.Fractional_power(1:max_ind), beam_power.Control_voltage_V_(1:max_ind));
hSI.hBeams.hBeams{1}.powerFraction2PowerWattLut = [0, beam_power.Power_W_(1); ...
    1, beam_power.Power_W_(max_ind)];

%%
result = struct;
result.data_x = beam_power.Control_voltage_V_;
result.data_y = beam_power.Power_W_;
result.fit_fun = @(c, v) c(1) .* cosd(c(2) .* v + c(3)).^2;
result.fit_coeff = lsqcurvefit(result.fit_fun, [max(result.data_y), 0.25, 90], result.data_x, result.data_y);
result.fit_y = result.fit_fun(result.fit_coeff, result.data_x);
result.fit_rho = corr(result.data_y, result.fit_y);
result.data_filepath = data_fp;
result.filepath = fullfile(DataManager.SCRIPT_PATH, 'Calibration/ImgPockelsCell', ...
    sprintf('ImgPockelPwrVsVol_%s.mat', datestr(now, 'yyyymmdd')));
save(result.filepath, '-struct', 'result');
%%
fig_hdl = figure;
ax_hdl = axes(fig_hdl);
scatter(ax_hdl, result.data_x, result.data_y, 'ro');
hold(ax_hdl, 'on');
plot_fit_hdl = plot(ax_hdl, result.data_x, result.fit_y, 'k', 'LineWidth', 1);
ax_hdl.XLabel.String = 'Voltage (V)';
ax_hdl.YLabel.String = 'Power (W)';
legend(ax_hdl, plot_fit_hdl, sprintf('\\rho = %.3f', result.fit_rho), 'Location', 'best');
ax_hdl.Title.String = sprintf('%.3e Cos^2[%.3e V + %.3e]', result.fit_coeff);

fig_fp = fullfile(calibration_folder, sprintf('Imaging_pockels_cell_pwr_vs_vol_%s.png', ...
    datestr(now, 'yyyymmdd')));
fun_print_image_in_several_formats(fig_hdl, fig_fp);
