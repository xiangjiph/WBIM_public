classdef WBIMAnalysisAblationTest < handle
    
    methods(Static) % Utilities
        
        function [fit_para] = fit_1d_int_profile_with_gaussian(x, int, opt)
            arguments
                x (:, 1) double
                int (:, 1) double
                opt.visQ (1, 1) logical = false;
                opt.normalized_Q (1, 1) logical = false;
                opt.vis_x_label = 'X';
                opt.vis_y_label = 'Intensity';
                opt.vis_title = [];
            end
            est_c = min(int);
            est_A = max(int) - est_c;
            
            [fit_para_0, ~] = gaussfitn(x, int, {est_c, est_A, [], []});
            fit_para_0 = cat(1, fit_para_0{:});
            fit_para = struct;
            fit_para.sigma = sqrt(fit_para_0(4));
            fit_para.mu = fit_para_0(3);
            fit_para.A = fit_para_0(2);
            fit_para.c = fit_para_0(1);
            fit_para.FWHM = 2 * sqrt(2 * log(2)) * fit_para.sigma;

            if opt.normalized_Q
                int = (int - fit_para.c) / fit_para.A;
                opt.vis_y_label = 'Normalized Intensity';
            end
            if opt.visQ

                fig_hdl = figure;
                ax_hdl = axes(fig_hdl);
                if opt.normalized_Q
                    data_hdl = plot(ax_hdl, x - fit_para_0(3), int, '-o');
                    hold(ax_hdl, 'on');
                    r_max = max(abs(x([1, end]) - fit_para_0(3)));
                    plt_x = linspace(-r_max, r_max, 1000);
                    fit_fun = @(para, x) para(1) + para(2) * exp(- (x).^2 / (2 * para(4)));
                    fit_hdl = plot(ax_hdl, plt_x, fit_fun([0, 1, 0, fit_para_0(4)], plt_x));
                else
                    data_hdl = plot(ax_hdl, x, int, '-o');
                    hold(ax_hdl, 'on');
                    plt_x = linspace(x(1), x(end), 1000);
                    fit_fun = @(para, x) para(1) + para(2) * exp(- (x - para(3)).^2 / (2 * para(4)));
                    fit_hdl = plot(ax_hdl, plt_x, fit_fun(fit_para_0, plt_x));
                end
                ax_hdl.XLabel.String = opt.vis_x_label;
                ax_hdl.YLabel.String = opt.vis_y_label;
                legend(ax_hdl, [data_hdl, fit_hdl], {'Data', sprintf('Gaussian FWHM = %.2f \\mum', ...
                    fit_para.FWHM)}, 'Location', 'best');
                if ~isempty(opt.vis_title)
                    ax_hdl.Title.String = opt.vis_title;
                end
                if opt.normalized_Q
                    ax_hdl.YLim = [0, 1];
                end
                fit_para.fig_hdl = fig_hdl;
            end
        end
    end
    
    methods(Static)
        function [result] = analyze_single_ablation_site(im_c, opt)
            arguments
                im_c
                opt.pixel_size_um (1, 3) double = [1, 1, 1]; % in [y, x, z] order
                opt.visQ (1, 1) logical = false;
                opt.save_fig_Q (1, 1) logical = false;
                opt.save_folder = [];
                opt.save_im_prefix = []; 
                opt.normalized_Q (1, 1) logical = false;
            end
            if isscalar(opt.pixel_size_um)
                opt.pixel_size_um = [1,1,1] * opt.pixel_size_um;
            end
            
            % Get the orthogonal view
            im_c_zxy = permute(im_c, [3, 2, 1]);
            
            % Find the center on zx plane by mean projection
            avg_int_zx = mean(im_c_zxy, 3);
            % Find the peak
            [~, max_ind] = max(avg_int_zx, [], 'all', 'linear');
            max_sub = fun_ind2sub(size(avg_int_zx), max_ind);
                        
            % Fit gaussian
            x_profile = avg_int_zx(max_sub(1), :);
            x_pos_um = (1 : numel(x_profile)) * opt.pixel_size_um(2);
            fit_para = WBIMAnalysisAblationTest.fit_1d_int_profile_with_gaussian(...
                x_pos_um, x_profile, 'normalized_Q', opt.normalized_Q,...
                'visQ', opt.visQ, 'vis_x_label', 'X (\mum)');

            
            z_profile = avg_int_zx(:, max_sub(2));
            z_pos_um = (1 : numel(z_profile)) * opt.pixel_size_um(3);
            fit_para_z = WBIMAnalysisAblationTest.fit_1d_int_profile_with_gaussian(...
                z_pos_um, z_profile, 'normalized_Q', opt.normalized_Q,...
                'visQ', opt.visQ, 'vis_x_label', 'Z (\mum)');
            
            % Average around the focal plane:
            z_avg_range = max_sub(1) + [-1, 1] * ceil(fit_para_z.FWHM / 2 / opt.pixel_size_um(3));
            z_avg_range = max(1, z_avg_range);
            x_avg_range = max_sub(2) + [-1, 1] * ceil(fit_para.FWHM / 2 / opt.pixel_size_um(2));
            avg_int_yx = mean(im_c(:, :, z_avg_range(1) : z_avg_range(2)), 3);
            y_profile = mean(avg_int_yx(:, x_avg_range(1) : x_avg_range(2)), 2);
            y_pos_um = (1 : numel(y_profile)) * opt.pixel_size_um(1);
            fit_para_y = WBIMAnalysisAblationTest.fit_1d_int_profile_with_gaussian(...
                y_pos_um, y_profile, 'normalized_Q', opt.normalized_Q,...
                'visQ', opt.visQ, 'vis_x_label', 'Y (\mum)');

            
            if opt.visQ
                fig_hdl = figure;
                ax_1 = subplot(1,2,1);
                imagesc(avg_int_zx);
                ax_1.DataAspectRatio = [1, opt.pixel_size_um(1) / opt.pixel_size_um(3), 1];
                ax_1.XTickLabel = arrayfun(@(x) num2str(x * opt.pixel_size_um(2)), ax_1.XTick, 'UniformOutput', false);
                ax_1.XLabel.String = 'X (\mum)';
                ax_1.YLabel.String = 'Z (\mum)';
                
                ax_2 = subplot(1,2,2);
                imagesc(avg_int_yx);
                ax_2.DataAspectRatio = [1,1,1];
                ax_2.XTickLabel = arrayfun(@(x) num2str(x * opt.pixel_size_um(2)), ax_2.XTick, 'UniformOutput', false);
                ax_2.YTickLabel = arrayfun(@(x) num2str(x * opt.pixel_size_um(1)), ax_2.YTick, 'UniformOutput', false);
                ax_2.XLabel.String = 'X (\mum)';
                ax_2.YLabel.String = 'Y (\mum)';

            end
            
            result = struct;
            result.x = fit_para;
            result.y = fit_para_y;
            result.z = fit_para_z;
            if opt.save_fig_Q
                if opt.normalized_Q
                    opt.save_im_prefix = sprintf('%s_normalized', opt.save_im_prefix);
                end

                fig_fp = fullfile(opt.save_folder, sprintf('%s_psfx.png', ...
                    opt.save_im_prefix));
                fun_print_image_in_several_formats(fit_para.fig_hdl, fig_fp);
                
                
                fig_fp = fullfile(opt.save_folder, sprintf('%s_psfz.png', ...
                    opt.save_im_prefix));
                fun_print_image_in_several_formats(fit_para_z.fig_hdl, fig_fp);
                
                fig_fp = fullfile(opt.save_folder, sprintf('%s_psfy.png', ...
                    opt.save_im_prefix));
                fun_print_image_in_several_formats(fit_para_y.fig_hdl, fig_fp);
                
                fig_fp = fullfile(opt.save_folder, sprintf('%s_mean_proj.png', ...
                    opt.save_im_prefix));
                fun_print_image_in_several_formats(fig_hdl, fig_fp);
            end            
        end
        
        
    end
    
    
    
    
end