calibration_folder = 'D:\My Drive\Microscope\Calibration\AblationLaserDivergence';
beam_power = readtable(fullfile(calibration_folder, 'Astrella_beam_size.xlsx'));
%%
power_n = beam_power.Position2;
iris_diameter_mm = beam_power.IrisDiamter_mm_;
%%
power_n = power_n ./ power_n(end);
fit_fun = @(s, d)1 - exp(-d.^2 ./ (8*s^2));
fit_s = lsqcurvefit(fit_fun, 6, iris_diameter_mm, power_n);
fprintf('sigma: %.2f mm\n', fit_s);
fig_hdl = figure;
ax_hdl = axes(fig_hdl);
scatter(ax_hdl, iris_diameter_mm, power_n);
hold(ax_hdl, 'on');
plt_x = 0 :0.5: max(iris_diameter_mm);
plot(ax_hdl, plt_x, fit_fun(fit_s, plt_x));
ax_hdl.Title.String = sprintf('\\sigma = %.2f mm, 1/e^2 diameter = %.2f mm', fit_s, fit_s * 4);
ax_hdl.XLabel.String = 'Iris diameter (mm)';
ax_hdl.YLabel.String = 'Normalized intensity';