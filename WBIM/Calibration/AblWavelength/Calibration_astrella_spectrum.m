data_folder = 'D:\My Drive\Microscope\Calibration\AblationLaserSpectrum';
data_fn = 'Astrella_spectrum.xlsx';
data_fp = fullfile(data_folder, data_fn);
spectrum_data = readtable(data_fp);
%%
is_selected_Q = spectrum_data.Wavelength_nm_ > 700 & ...
    spectrum_data.Wavelength_nm_ < 900;
wavelength_nm = spectrum_data.Wavelength_nm_(is_selected_Q);
intensity = spectrum_data.Intensity(is_selected_Q);

fig_hdl = figure;
ax_hdl_1 = subplot(2, 1, 1);
plt_hdl_s = plot(ax_hdl_1, wavelength_nm, intensity, 'LineWidth', 2);
ax_hdl_1.XLabel.String = 'Wavelength (nm)';
ax_hdl_1.YLabel.String = 'Intensity';
grid(ax_hdl_1, 'on')
%% Fit with Gaussian
f_guass = @(x, xdata)x(1) .* exp(- (xdata - x(2)).^2 ./ (2 .* x(3).^2)) + x(4);
x0 = [4000, 800, 10, 160];
[x, resnorm, ~, exitflag, foutput] = lsqcurvefit(f_guass, x0, ...
    wavelength_nm, intensity);
hold(ax_hdl_1, 'on');
y_list = f_guass(x, wavelength_nm);
plt_hdl_f = plot(ax_hdl_1, wavelength_nm, y_list, 'LineWidth', 1.5);
legend(ax_hdl_1, [plt_hdl_s, plt_hdl_f], 'Spectrum',...
    sprintf('Gaussian fit\n\\lambda_0 = %.2f nm\n\\sigma = %.2f nm', x(2), x(3)));
%% Enclosed energy


%%
fig_fp = fullfile(data_folder, 'Astrella_spectrum_w_gauss_fit.png');
% fun_print_image_in_several_formats(fig_hdl, fig_fp);

corrcoef(intensity, y_list)
%%
grating_density = 830;
lambda_0_nm = 797.21;
incident_angle_rad = lambda_0_nm * 1e-3 / (1e3 / grating_density);
lambda_fwhm_nm = x(3) * 2 * sqrt(2 * log(2));
diff_angle_rad = (wavelength_nm - lambda_0_nm) * 1e-3 / (1e3 / grating_density);
diff_angle_deg = diff_angle_rad * 180 / pi;

dist_to_grating_mm = 150 + 30 * 25.4;
diff_dist_to_ctr_mm = diff_angle_rad * dist_to_grating_mm;
beam_sigma_mm = 1 / cos(incident_angle_rad);

dist_to_ctrl_mm = linspace(diff_dist_to_ctr_mm(1) - 4 * beam_sigma_mm, ...
    diff_dist_to_ctr_mm(end) + 4 * beam_sigma_mm, 1000).';
monochromatic_gauss = exp(- pdist2(diff_dist_to_ctr_mm, dist_to_ctrl_mm) .^ 2 /...
    (2 * beam_sigma_mm ^ 2)) .* (intensity - x(4));
overall_gauss = sum(monochromatic_gauss, 1);

fig_hdl = figure;
ax_hdl = axes(fig_hdl);
plt_hdl = plot(ax_hdl, dist_to_ctrl_mm, overall_gauss);


