classdef WBIMCalibration < handle
    
   properties
       
       
   end
   %%
   methods
       function obj = WBIMCalibration()
           
       end
   end
   %%
   methods(Static)
       function result = FitPowerVsAngle(rot_angle_deg, pwr)
           assert(numel(rot_angle_deg) == numel(pwr));
           if isrow(rot_angle_deg)
               rot_angle_deg = rot_angle_deg.';
           end
           if isrow(pwr)
              pwr = pwr.'; 
           end
           result.angle = rot_angle_deg;
           result.power = pwr;
           result.fit_fun = @(c, theta) c(1) .* cosd(2 * (theta + c(2))).^2;
           result.fit_coeff = lsqcurvefit(result.fit_fun, [max(pwr), 0], result.angle, result.power);
           result.fit_y = result.fit_fun(result.fit_coeff, result.angle);
           result.fit_rho = corr(result.power, result.fit_y);           
           result.data_filepath = [];
           result.filepath = [];
       end
       
       function [fig_hdl, ax_hdl] = PlotPowerVsAngleFit(result)
           fig_hdl = figure;
           ax_hdl = axes(fig_hdl);
           scatter(ax_hdl, result.angle, result.power, 'ro');
           % Fit with cos
           hold(ax_hdl, 'on');
           plot_fit_hdl = plot(ax_hdl, result.angle, result.fit_y, 'k', 'LineWidth', 1);
           if result.fit_coeff(2) > 0
               ax_hdl.Title.String = sprintf('%.2f Cos(2(\\theta + %.2f))^2', result.fit_coeff);
           else
               ax_hdl.Title.String = sprintf('%.2f Cos(2(\\theta - %.2f))^2', abs(result.fit_coeff));
           end
           legend(ax_hdl, plot_fit_hdl, sprintf('\\rho = %.3f', result.fit_rho), 'Location', 'best');
           grid(ax_hdl, 'on');
       end
   end
end