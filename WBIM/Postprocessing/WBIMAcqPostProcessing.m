classdef WBIMAcqPostProcessing < handle
    %% To do list:
    % 1. *** Add chunk size 
    % 2. * Use low-level hdf5    
    properties(Abstract, Hidden, Constant)
        COMPONENT_NAME
    end
    
    properties(Abstract)
        
    end
    
    properties
%         tile_info WBIMTileMetadata
        log_filepath char
        err_log_filepath char
    end
    properties(Hidden, Dependent)
       error_logger_name 
    end
    
    properties(Hidden)
       h_logger mlog.Logger 
       h_err_logger mlog.Logger
    end
    %%
    methods
        function val = get.error_logger_name(obj)
           val = sprintf('%s_err', obj.COMPONENT_NAME);
        end
    end
    %%
    methods
        function obj = WBIMAcqPostProcessing()
            
        end
        
        function init_logger(obj)
            [log_folder, ~] = fileparts(obj.log_filepath);
            if ~isfolder(log_folder)
                mkdir(log_folder);
            end
            obj.h_logger = mlog.Logger(obj.COMPONENT_NAME, obj.log_filepath);
            obj.h_logger.FileThreshold = mlog.Level.DEBUG;
            obj.h_logger.CommandWindowThreshold = mlog.Level.WARNING;
            obj.h_logger.MessageReceivedEventThreshold = mlog.Level.NONE;
            
            obj.h_err_logger = mlog.Logger(obj.error_logger_name, obj.err_log_filepath);
            obj.h_err_logger.FileThreshold = mlog.Level.DEBUG;
            obj.h_err_logger.CommandWindowThreshold = mlog.Level.DEBUG;
            obj.h_err_logger.MessageReceivedEventThreshold = mlog.Level.NONE;
        end
        
    end
    %% 
    methods(Static)
        function stat_str = analyze_single_tile_mip(im)
            assert(ismatrix(im));
            im_max_val = single(intmax(class(im)));
            im = single(im);
            im_f = medfilt2(im);
            stat_str = struct;
            stat_str.size = size(im);
            stat_str.num_pxl = prod(stat_str.size);
            stat_str.mean_raw = mean(im, 'all');
            stat_str.std_raw = std(im, 1, 'all');
            stat_str.max_raw = max(im, [], 'all');
            stat_str.frac_above_half_max_raw = nnz(im > im_max_val/2) / stat_str.num_pxl;
            stat_str.frac_saturated_raw = nnz(im >= (im_max_val * 0.95)) / stat_str.num_pxl;
            
            stat_str.mean_mft = mean(im_f, 'all');
            stat_str.std_mft = std(im_f, 1, 'all');
            stat_str.max_mft = max(im_f, [], 'all');
            stat_str.frac_above_half_max_mft = nnz(im_f > im_max_val/2) / stat_str.num_pxl;
            stat_str.frac_saturated_mft = nnz(im_f >= (im_max_val * 0.95)) / stat_str.num_pxl;
            
            stat_str.mean_rf = stat_str.mean_mft / (eps + stat_str.mean_raw);
            stat_str.std_rf = stat_str.std_mft / (eps + stat_str.std_raw);
            stat_str.max_rf = stat_str.max_mft / (eps + stat_str.max_raw);
            stat_str.frac_above_half_max_rf = stat_str.frac_above_half_max_mft ...
                / (eps + stat_str.frac_above_half_max_raw);
            stat_str.frac_saturated_rf = stat_str.frac_saturated_mft...
                / (eps + stat_str.frac_saturated_raw);
            
            % Normalized median-filtered image intensity PDF
            hist_bin_edge = 0 : 0.1 : 1;
            stat_str.nmpdf = histcounts(rescale(im_f), hist_bin_edge, 'Normalization', 'probability');
            stat_str.mpdf = histcounts(im_f ./ im_max_val, hist_bin_edge, 'Normalization', 'probability');
        end
        
        function exit_code = write_si_metadata(fp, info)
            [folder, ~, ~] = fileparts(fp);
            if ~isfolder(folder)
                mkdir(folder);
            end
            f_hdl = fopen(fp, 'a+');
            fprintf(f_hdl, '%s', info);
            exit_code = fclose(f_hdl);
        end
    end    
    %% Mask
    methods(Static)
        %% Computation
            %% Utilities
        function [im_ds] = smooth_and_downsample_sections_in_stack(im, downsample_facotr)
            if all(downsample_facotr == 1)
                im_ds = im;
            else
                im_size = size(im, [1,2]);
                num_sec = size(im, 3);
                if isscalar(downsample_facotr)
                    downsample_facotr = [downsample_facotr, downsample_facotr];
                end
                im_size_ds = round(im_size ./ downsample_facotr);
                im_ds = cell(num_sec, 1);
                for i = 1 : num_sec
                    tmp_im = imgaussfilt(im(:, :, i), downsample_facotr ./2);
                    im_ds{i} = imresize(tmp_im, im_size_ds);
                end
                im_ds = cat(3, im_ds{:});
            end

        end
        
        function [mask_str] = construct_ablation_mask_yx_um_itp(mask, voxel_size_um, ...
                options)
            arguments
               mask
               voxel_size_um
               options.interpolation = 'nearest';
               options.extrapolation = 'nearest';
            end            
            mask_str = fun_initialized_structure_array_with_fieldname_list(...
                {'mask_size', 'mask_size_um'});
            if ~isempty(mask)
                array_size = size(mask);
                sub_1 = (0.5 : array_size(1) - 0.5) * voxel_size_um(1);
                sub_2 = (0.5 : array_size(2) - 0.5) * voxel_size_um(2);
                sub_3 = (0.5 : array_size(3) - 0.5) * voxel_size_um(3);
                [sub_1, sub_2, sub_3] = ndgrid(sub_1, sub_2, sub_3);
                mask_str = struct;
                mask_str.itp = griddedInterpolant(sub_1, sub_2, sub_3, ...
                    single(mask), options.interpolation, options.extrapolation);
                mask_str.mask_size = array_size;
                mask_str.mask_size_um = array_size .* voxel_size_um;
            end
        end
        
        function ds_smip_mask = construct_z_fwhm_mask(result)
            ds_smip_mask = false([result.num_layers, prod(result.im_size)]);
            for i = 1 : result.num_pxl
                tmp_hm_idx = (result.fitst_lt_hm_idx(i)) : (result.first_gt_hm_idx(i) - 1);
                ds_smip_mask(tmp_hm_idx, result.valid_ind(i)) = true;
            end
            ds_smip_mask = reshape(ds_smip_mask, result.num_layers, result.im_size(1), result.im_size(2));
            ds_smip_mask = permute(ds_smip_mask, [2,3,1]);
        end
        
        function [result] = compute_z_int_profile_fractional_mask(ds_smip_5um,...
                ds_mask, options)
            arguments
                ds_smip_5um {mustBeNumeric} % 3D image array
                ds_mask logical
                options.int_frac (1,1) double = 0.5;
                options.est_bg = []                
            end
            
            if isempty(options.est_bg)
                if ~all(ds_mask, 'all')
                    ds_mip = max(ds_smip_5um, [], 3);
                    options.est_bg = median(ds_mip(~ds_mask));
                else
                    options.est_bg = 0;
                end
            end
            result = struct;
            
            result.im_size = size(ds_smip_5um, [1,2]);
            result.num_layers = size(ds_smip_5um, 3);
            result.valid_ind = find(ds_mask);
            result.num_pxl = numel(result.valid_ind);
            
            ds_smip_5um = permute(ds_smip_5um, [3,1,2]);
            ds_smip_5um = reshape(ds_smip_5um, result.num_layers, []);
            valid_profiles = ds_smip_5um(:, result.valid_ind);
            int_max = max(valid_profiles, [], 1);
            est_half_int = (int_max - options.est_bg ) * options.int_frac + options.est_bg ;
            above_hm_Q = valid_profiles > est_half_int;
            % Find the first value lower than half intensity
            hm_Q_diff = diff(cat(1, false(1, result.num_pxl), above_hm_Q, ...
                false(1, result.num_pxl)));
            % Return the first maximum value position
            [max_Q_value, result.fitst_lt_hm_idx] = max(hm_Q_diff, [], 1); 
            % Return the first minimum value position
            [min_Q_value, result.first_gt_hm_idx] = min(hm_Q_diff, [], 1); 
            
            assert(all(max_Q_value == 1) & all(min_Q_value == -1));
            assert(all(result.fitst_lt_hm_idx < result.first_gt_hm_idx), ...
                'The first point with intensity less than half max should occurs before the first point with intensity higher than half max');
            result.mask = WBIMAcqPostProcessing.construct_z_fwhm_mask(result);
            result.depth_map = uint8(sum(result.mask, 3));
        end
          
        function margin_stat = estimate_image_background_on_the_edge(im, margin_num_pxls, ...
                void_value)
            if nargin < 3
                void_value = nan;
            end            
            im_dim = ndims(im);
            if im_dim == 3
                im = max(im, [], 3);
            end
            % Assuming the stitched image has nan 
            if isnan(void_value)
                is_void_Q = isnan(im);
            elseif isfinite(void_value)
                is_void_Q = (im == void_value);
            end
            if ~any(is_void_Q, 'all')
%                 warning('None of the pixel has value equal to the specificed voiud value');                
                margin_stat = fun_analysis_get_basic_statistics([]);
                margin_stat.th_5_sigma = [];
                margin_stat.th_3_sigma = [];
                margin_stat.th_3_ipt = [];
            else
                is_void_Q_p = padarray(is_void_Q, [1, 1], true);
                is_margin_Q = imdilate(is_void_Q_p, strel('disk', margin_num_pxls)) & ...
                    ~is_void_Q_p;
                is_margin_Q = is_margin_Q(2:end-1, 2:end-1);
                margin_values = im(is_margin_Q);
                assert(all(isfinite(margin_values)), 'All values in the margin band should be finite');
                margin_stat = fun_analysis_get_basic_statistics(margin_values);
                margin_stat.th_5_sigma = margin_stat.mean + margin_stat.std * 5;
                margin_stat.th_3_sigma = margin_stat.mean + margin_stat.std * 3;
                margin_stat.th_3_ipt = margin_stat.median + 3 * diff(margin_stat.prctile_val([7,9]));
            end
        end
        
        function bg_stat = estimate_image_background_with_prior_th(im, est_th)
            validateattributes(im, {'numeric'}, {});
            validateattributes(est_th, {'numeric'}, {'scalar'});
            est_bg_mask = im < est_th;
            est_bg_int = im(est_bg_mask);
            bg_stat = fun_analysis_get_basic_statistics(est_bg_int);
            bg_stat.th_5_sigma = bg_stat.mean + bg_stat.std * 5;
            bg_stat.th_3_sigma = bg_stat.mean + bg_stat.std * 3;
            bg_stat.itp = diff(bg_stat.prctile_val([7,9]));
            bg_stat.th_3_ipt = bg_stat.median + 3 * bg_stat.itp;
        end
        
        function [major_mask, varargout] = remove_small_cc_by_distance(im_mask, min_major_ccf, ...
                max_sm_cc_dist_dspxl)
            validateattributes(im_mask, 'logical', {'2d'});
            validateattributes(min_major_ccf, 'numeric', {'scalar'});
            validateattributes(max_sm_cc_dist_dspxl, 'numeric', {'scalar'});
            
            cc = bwconncomp(im_mask);
            cc_size = cellfun(@numel, cc.PixelIdxList);
            cc_frac = cc_size / sum(cc_size);
            is_major_cc_Q = cc_frac > min_major_ccf;
            if any(is_major_cc_Q)                
                major_mask = false(size(im_mask));
                cc.major_cc_ind = cc.PixelIdxList(is_major_cc_Q);
                major_mask(cat(1, cc.major_cc_ind{:})) = true;
                % Add the small components based on distance?
                mm_dist = bwdist(major_mask);
                
                cc.small_cc_ind = cc.PixelIdxList(~is_major_cc_Q);
                % Average distnace to the nearest major connected components
                small_cc_avg_dist = cellfun(@(x) min(mm_dist(x)), ...
                    cc.small_cc_ind);
                close_smcc_Q = (small_cc_avg_dist < max_sm_cc_dist_dspxl);
                cc.small_cc_ind = cc.small_cc_ind(close_smcc_Q);
                major_mask(cat(1, cc.small_cc_ind{:})) = true;
                cc.NumObjects = numel(cc.major_cc_ind) + numel(cc.small_cc_ind);
            else
                major_mask = im_mask;
                fprintf('No connected component account for > %.2f of the total cc area\n', min_major_ccf);
            end
            if nargout > 1
                varargout{1} = cc;
            end
        end
        
        function mask_sample_yx_um = compute_interpolated_ablation_mask(...
                abl_mask_itp, abl_z_r_um)
            arguments
                abl_mask_itp
                abl_z_r_um (1, :) double
            end            
            [sub1, sub2, sub3] = ndgrid(1 : abl_mask_itp.mask_size_um(1), ...
                1 : abl_mask_itp.mask_size_um(2), abl_z_r_um);
            mask_sample_yx_um = (abl_mask_itp.itp(sub1, sub2, sub3) > 0.5);
        end
        
        function result = compute_step_mip_fractional_max_mask(stitched_mip, ...
                frac_max, options)
            arguments
               stitched_mip
               frac_max (1,1) double = 0.5
               options.est_bg_th = []
               options.bg_est_margin_width_pxl = 5;
               options.visQ (1,1) logical = false;
               options.expand_disk_r_pxl = 2;
               options.clean_disk_r_pxl = 4;
               options.remove_distant_small_cc_Q (1,1) logical = true;
               options.rdscc_min_major_cc_frac (1,1) double = 0.1;
               options.rdscc_min_dist_pxl (1,1) double = 40;
            end
            result = struct;
            result = fun_initialized_structure_array_with_fieldname_list({'im_size', ...
                'num_layers', 'valid_ind', 'num_pxl', 'int_max', 'fitst_lt_hm_idx', ...
                'first_gt_hm_idx', 'mask', 'depth_map'});
            result.im_size = size(stitched_mip, [1,2]);
            result.num_layers = size(stitched_mip, 3);
            
            mip_0 = max(stitched_mip, [], 3);
            margin_stat = WBIMAcqPostProcessing.estimate_image_background_on_the_edge(...
                mip_0, options.bg_est_margin_width_pxl);
            if ~isempty(margin_stat.th_3_ipt) 
                if isempty(options.est_bg_th)
                    options.est_bg_th = margin_stat.th_3_ipt;
                else
                    options.est_bg_th = max(options.est_bg_th, margin_stat.th_3_ipt);
                end                
            else
                assert(~isempty(options.est_bg_th), 'Failed to estimated background intensity threshold. Please provide an estimation');
            end
            mip_mask = mip_0 > options.est_bg_th;
            if any(mip_mask, 'all')
                if options.remove_distant_small_cc_Q
                    [mip_mask, ~] = WBIMAcqPostProcessing.remove_small_cc_by_distance(...
                        mip_mask, options.rdscc_min_major_cc_frac, ...
                        options.rdscc_min_dist_pxl);
                end
                result.valid_ind = find(mip_mask);
                result.num_pxl = numel(result.valid_ind);
                
                stitched_mip = permute(stitched_mip, [3,1,2]);
                stitched_mip = reshape(stitched_mip, result.num_layers, []);
                valid_profiles = stitched_mip(:, result.valid_ind);
                result.int_max = max(valid_profiles, [], 1);
                
                if ~isempty(margin_stat.mean)
                    est_bg = margin_stat.mean;
                else
                    est_bg = median(mip_0(~mip_mask), 'omitnan');
                end
                % Analyze intensity profile along z direction
                est_half_int = (result.int_max - est_bg) * frac_max + est_bg ;
                above_hm_Q = valid_profiles > est_half_int;
                % Find the first value lower than half intensity
                hm_Q_diff = diff(cat(1, false(1, result.num_pxl), above_hm_Q, ...
                    false(1, result.num_pxl)));
                % Return the first maximum value position
                [max_Q_value, result.fitst_lt_hm_idx] = max(hm_Q_diff, [], 1);
                % Return the first minimum value position
                [min_Q_value, result.first_gt_hm_idx] = min(hm_Q_diff, [], 1);
                
                assert(all(max_Q_value == 1) & all(min_Q_value == -1));
                assert(all(result.fitst_lt_hm_idx < result.first_gt_hm_idx), ...
                    'The first point with intensity less than half max should occurs before the first point with intensity higher than half max');
                result.mask = WBIMAcqPostProcessing.construct_z_fwhm_mask(result);
                if options.expand_disk_r_pxl > 0
                    result.mask = imdilate(result.mask, strel('disk', options.expand_disk_r_pxl));
                end
                if options.clean_disk_r_pxl > 0
                    result.mask = imclose(result.mask, strel('disk', options.clean_disk_r_pxl));
                end
                result.depth_map = uint8(sum(result.mask, 3));
                result.fitst_lt_hm_idx = uint8(result.fitst_lt_hm_idx);
                result.first_gt_hm_idx = uint8(result.first_gt_hm_idx);
                if options.visQ
                    WBIMAcqPostProcessing.vis_abl_thickness_map(result);
                end                
            end
        end
        
        function [shg_in_vsl_mask, vsl_in_shg_mask] = compute_step_mip_spectrum_unmixing_mask(...
                mip_vsl, mip_shg, options)
            arguments
                mip_vsl
                mip_shg
                options.bg_offset (1,1) double = 500;
                options.visQ (1,1) logical = false;
                options.r_vsl2shg (1,1) double = 0.5; % maximum breakthrough intensity ratio from vsl to shg
                
                options.r_shg2vsl (1,1) double = 0.5; % maximum breakthrough intensity ratio from shg to vsl 
                options.min_int_shg (1,1) double = 5000;
                options.min_int_vsl (1,1) double = 5000;
            end
            num_sec = size(mip_vsl, 3);
            [shg_in_vsl_mask, vsl_in_shg_mask] = deal(false(size(mip_vsl)));
            for i = 1 : num_sec
                tmp_vsl = mip_vsl(:, :, i);
                tmp_shg = mip_shg(:, :, i);
                tmp_vsl_mask = tmp_vsl > options.min_int_vsl;
                tmp_shg_mask = tmp_shg > options.min_int_shg;
                
                sng_r_vsl2shg = (tmp_vsl + options.bg_offset) ./ (tmp_shg + options.bg_offset);
                
                shg_in_vsl_mask(:, :, i) = tmp_shg_mask & ...
                    (sng_r_vsl2shg < options.r_shg2vsl) & (tmp_vsl > 0);
                
                vsl_in_shg_mask(:, :, i) = tmp_vsl_mask & ...
                    (sng_r_vsl2shg > (1/options.r_vsl2shg)) & ...
                    (sng_r_vsl2shg > 1) & (tmp_shg > 0);
            end            
        end
        %% Surface detection
        % Not used. To be modified
        function [stat_str, varargout] = compute_spectrum_unmixing_factor(im1, im2, im2_th, ...
                fit_range, visQ)
            % Estimate the break-through ratio of signal in im2 to im1
            arguments
                im1
                im2
                im2_th (1,:) double
                fit_range (1,2) double
                visQ (1,1) logical = false;
            end
            switch numel(im2_th)
                case 1
                    im2_mask = im2 > im2_th;
                case 2
                    im2_mask = im2 > im2_th(1) & im2 < im2_th(2);
                otherwise
                    error('im2_th can only have 1 or 2 elements');
            end            
            im1_data = im1(im2_mask);
            im2_data = im2(im2_mask);
            im2_bin_edge = linspace(fit_range(1), fit_range(2), 50);
            [im2_ind, im2_bin_val] = fun_bin_data_to_idx_list_by_edges(im2_data, im2_bin_edge);
            valid_bin_Q = ~cellfun(@isempty, im2_ind);
            im2_ind = im2_ind(valid_bin_Q);
            im2_bin_val = im2_bin_val(valid_bin_Q);
            im1_binned_stat = fun_analysis_get_basic_statistics_in_bins(im1_data, im2_ind);
            dm_lf = fitlm(im2_bin_val, [im1_binned_stat.median]);
            stat_str.k = dm_lf.Coefficients.Estimate(2);
            stat_str.b = dm_lf.Coefficients.Estimate(1);
            stat_str.R2 = dm_lf.Rsquared.Ordinary;
            if visQ
                fig_hdl = figure;
                fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* 1.5;
                ax_0 = nexttile();
                imagesc(ax_0, im2);
                ax_0.Title.String = 'Im2';
                ax_1 = nexttile();
                imagesc(ax_1, im1);
                ax_1.Title.String = 'Im1';
                ax_2 = nexttile();
                imshowpair(im1, im2_mask);
                ax_2.Title.String = 'Im1 image overlaid with IM2 mask';
                ax_3 = nexttile();
                histogram2(ax_3, im2_data, im1_data, 'DisplayStyle', 'tile');
                ax_3.YLim(2) = 5e3;
                hold(ax_3, 'on');
                plot(ax_3, im2_bin_val, [im1_binned_stat.median], 'rx');                    
                line_hdl = plot(ax_3, im2_bin_val, im2_bin_val * stat_str.k + ...
                    stat_str.b, 'LineWidth', 2);
                legend(ax_3, line_hdl, sprintf('k = %.3f\nR^2 = %.3f\n', stat_str.k, ...
                    stat_str.R2));                
                ax_3.XLabel.String = 'Im2 intensity';
                ax_3.YLabel.String = 'Im1  intensity';
                ax_3.Title.String = 'Joint distribution of masked pixel intensity';
                ax_4 = nexttile();
%                 imagesc(ax_4, max(0, im1 - im2 * stat_str.k));
                imagesc(ax_4, max(0, im1 - im2 * 0.1));
                ax_4.Title.String = 'Im1 image after linear unmixing';
                [ax_0.DataAspectRatio, ax_1.DataAspectRatio, ax_2.DataAspectRatio, ax_4.DataAspectRatio] = deal([1,1,1]);                
            else
                fig_hdl = [];
            end
            if nargout > 1
                varargout{1} = fig_hdl;
            end
        end
        
        function [dIdz_n, varargout] = compute_z_normalized_derivative(smip_array)
            % This function computes the normalized derivative along the
            % 3rd dimension. dI(z)/dz / I(z).
            % Compute forward / backward difference at the ends, use
            % center difference for the middle points
            
            z_diff = cat(3, smip_array(:, :, 2) - smip_array(:, :, 1) + eps, ...
                (smip_array(:, :, 3:end) - smip_array(:, :, 1:end-2) + eps) ./2, ...
                smip_array(:, :, end) - smip_array(:, :, end-1) + eps);
            dIdz_n = z_diff ./ (smip_array + eps);
            if nargout > 1
                varargout{1} = z_diff;
            end
        end
        
        % Not used
        function [max_ls_um] = compute_z_int_scatter_length_from_norm_grad(smip_array, ...
                z_pxl_size_um)
            arguments
                smip_array
                z_pxl_size_um (1,1) double = 1
            end
            [dIdz_n, dIdz] = WBIMAcqPostProcessing.compute_z_normalized_derivative(...
                smip_array);
            max_ls_um = -z_pxl_size_um./dIdz_n;
            % Remove all the postive gradient 
            max_ls_um(max_ls_um<=0) = nan;            
        end
        
        % Not used
        function [fit_str, varargout] = compute_z_int_scatter_length_from_exp_fit(int_trace, z_trace_um, ...
                options)
            arguments
                int_trace (:, 1) 
                z_trace_um (:, 1) = 1 : numel(int_trace);
                options.saturation_int (1,1) = 31000
                options.background_int (1,1) = 3000;
                options.visQ (1,1) logical = false;
            end
            if ~isfloat(int_trace)
                int_trace = double(int_trace);
            end
            int_trace(int_trace >= options.saturation_int) = nan;
            [int_max, max_idx] = max(int_trace);
            % Find the first descending trace after the peak 
            selectedQ = true(size(int_trace));
            selectedQ(1:(max_idx-1)) = false;
            selectedQ = selectedQ & (int_trace >= options.background_int);
            % Fit with exponential function
            fit_x = z_trace_um(selectedQ);
            fit_y = log(int_trace(selectedQ));
            
            fit_hdl = fitlm(fit_x, fit_y);
            fit_str = struct;
            fit_str.k = fit_hdl.Coefficients.Estimate(2);
            fit_str.ls_um = -1/fit_str.k;
            fit_str.b = fit_hdl.Coefficients.Estimate(1);
            fit_str.R2 = fit_hdl.Rsquared.Ordinary;
            
            y_hat = exp(fit_str.b + fit_str.k .* fit_x);
            if options.visQ
                fig_hdl = figure;
                ax_hdl = axes(fig_hdl);
                plot(ax_hdl, z_trace_um, int_trace);
                ax_hdl.YScale = 'log';
                hold(ax_hdl, 'on');
                plot(ax_hdl, fit_x, y_hat);
                grid(ax_hdl, 'on');
                ax_hdl.XLabel.String = 'Depth (\mum)';
                ax_hdl.YLabel.String = 'Intensity';
                legend(ax_hdl, 'Data',  sprintf('l_s = %.2f \\mum\nR^2 = %.2f',...
                    fit_str.ls_um, fit_str.R2));
            else
                fig_hdl = [];
            end
            if nargout > 1
                varargout{1} = fig_hdl;
            end
        end
        
        function result = analyze_stack_scattering_length_in_skull(tile_info, opts)
            arguments
                tile_info (1,1) WBIMTileMetadata
                opts.vsl_th (1,1) double = 6e3;
                opts.min_dist_2_vsl (1,1) double = 10;
                opts.shg_th (1,1) double = 5e3;
                opts.est_bg_int (1,1) double = 3e3;
                opts.est_saturation_int (1,1) double = 3.1e4;
                opts.valid_ls_range_um (1,2) double = [10, 100]
                opts.visQ (1,1) logical = false;
                opts.save_fig_Q (1,1) logical = false;
            end
            persistent DM
            if isempty(DM)
                DM = WBIMFileManager;
            end
            test_stacks = tile_info.load_tile([WBIMChannelName.Vessel, WBIMChannelName.SHG]);
            test_stacks = cellfun(@medfilt3, test_stacks, 'UniformOutput', false);
            % Downsample the image stack
            test_stacks = cellfun(@(x) imresize3(x, round(tile_info.stack_size_um)),...
                test_stacks, 'UniformOutput', false);
            test_mips = cellfun(@(x) max(x, [], 3), test_stacks, 'UniformOutput', false);
%             skull_pixels = ~(imdilate(test_mips{1} > opts.vsl_th, strel('disk', opts.min_dist_2_vsl)))...
%                 & (test_mips{2} > opts.shg_th);
            skull_pixels = (test_mips{2} > opts.shg_th);
            skull_ind = find(skull_pixels);
            result = struct('ls_map', [], 'ind', [], 'R2', [], 'stat', fun_analysis_get_basic_statistics([]));
            if ~isempty(skull_ind)
                num_skull_ind = numel(skull_ind);
                shg_z_profile = reshape(permute(test_stacks{2}, [3,1,2]), tile_info.stack_size(3), []);
                shg_z_profile = shg_z_profile(:, skull_ind);
                
                
                [ls_um, ls_fit_R2] = deal(nan(num_skull_ind, 1));
                profile_z_um = (1 : tile_info.stack_size_um(3));
                est_bg_max = opts.est_bg_int;
                est_saturated_int = opts.est_saturation_int;
                t_tic = tic;
                parfor i = 1 : num_skull_ind
                    warning('off', 'stats:LinearModel:RankDefDesignMat');
                    tmp_int = shg_z_profile(:, i);
                    if any(est_bg_max < tmp_int)
                        tmp_str = WBIMAcqPostProcessing.compute_z_int_scatter_length_from_exp_fit(...
                            tmp_int, profile_z_um, 'visQ', false, 'background_int', est_bg_max, ...
                            'saturation_int', est_saturated_int);
                        ls_um(i) = tmp_str.ls_um;
                        ls_fit_R2(i) = tmp_str.R2;
                    end
                    if mod(i, 1000) == 0
                        fprintf('Finish computing %d/%d profiles. Elapsed time is %.2f seconds\n',...
                            i, num_skull_ind, toc(t_tic));
                    end
                    warning('on', 'stats:LinearModel:RankDefDesignMat');
                end
                is_valid_Q = ls_um > opts.valid_ls_range_um(1) & ...
                    ls_um < opts.valid_ls_range_um(2);
                ls_um = ls_um(is_valid_Q);
                skull_ind = skull_ind(is_valid_Q);
                ls_fit_R2 = ls_fit_R2(is_valid_Q);
                
                result.ls_map = nan(round(tile_info.stack_size_um(1:2)));
                result.ls_map(skull_ind) = ls_um;
                result.ind = skull_ind;
                result.R2 = ls_fit_R2;
                result.stat = fun_analysis_get_basic_statistics(ls_um);
                
                if opts.visQ || opts.save_fig_Q
                    fig_hdl = figure('Visible', 'off');
                    fig_hdl.Position(3) = fig_hdl.Position(3) * 2;
                    ax_1 = subplot(1,2,1);
                    histogram2(ax_1, ls_um, result.R2, 'DisplayStyle', 'tile');
                    ax_1.Title.String = sprintf('Median: %.1f \\mum', result.stat.median);
                    ax_1.XLabel.String = 'Scattering length (\mum)';
                    ax_1.YLabel.String = 'R^2';
                    cbar = colorbar(ax_1);
                    cbar.Label.String = 'Number of points';
                    ax_1.YLim(1) = 0.80;
                    ax_1.ColorScale = 'linear';
                    ax_2 = subplot(1,2,2);
                    imagesc(ax_2, result.ls_map);
                    ax_2.DataAspectRatio = [1,1,1];
                    cbar_hdl = colorbar(ax_2);
                    ax_2.CLim = [0, 100];
                    cbar_hdl.Label.String = 'l_s (\mum)';
                    if opts.save_fig_Q
                        fig_fp = fullfile(DM.fp_experiment(tile_info.experiment_group, tile_info.experiment),...
                            'visualization', 'Scattering_length', ...
                            sprintf('%s_%s_layer_%d_grid_%d_scatter_length_map.png', ...
                            tile_info.experiment_group, tile_info.experiment, tile_info.layer, tile_info.grid_ind));
                        fun_print_image_in_several_formats(fig_hdl, fig_fp);
                    end
                    if opts.visQ
                        fig_hdl.Visible = 'on';
                    else
                        delete(fig_hdl);
                    end
                end
            end            
        end
        
        
        
        
        % Not used
        function [shg_mask_surf, varargout] = mask_refinement_surface_frac_z_peak_int(...
                smip_shg, shg_mask, options)
            
            arguments
                smip_shg
                shg_mask
                options.sec_range (1,2) double = [1, size(smip_shg, 3)];
                options.max_peak_int (1,1) double = inf;
                options.peak_int_fraction (1,1) double = 0.5;
                % minimum of the local intensity maximum
                options.min_lm_int (1,1) double = 7e2
                % Local maximum detection
                options.local_max_wd_sz (1,1) double = 3
                options.local_max_min_frac_max_int (1,1) double = 0.2;
                % Smooth surface map
                options.sf_map_imc_disk_r_pxl (1,1) double = 2;
                % Clean up ablation mask 
                options.abl_mask_min_cc_pxl_2d (1,1) double = 0;
                options.visQ (1,1) logical = false;
            end
            % Estimate background level. Use median or mean?
            est_bg_med = median(smip_shg(~shg_mask), 'omitnan');
            smip_shg(isnan(smip_shg)) = 0;        
            % Reshape array to get the z intensity profiles
            mask_yx_size = size(smip_shg, [1,2]);
            num_sections = size(smip_shg, 3);
            sec_idx_list = options.sec_range(1) : options.sec_range(2);
            search_map = any(shg_mask(:, :, sec_idx_list), 3);
            mask_pxl_ind = find(search_map);            
            num_pxl = numel(mask_pxl_ind);
            % Get masked pixel z intensity profile 
            pxl_int_profile = reshape(permute(smip_shg, [3,1,2]), num_sections, []);
            pxl_int_profile = pxl_int_profile(:, mask_pxl_ind);
            % Search for the local maxima in the intensity profile 
            int_lm_Q = movmax(pxl_int_profile, options.local_max_wd_sz, 1);
            % Local maximum should be higher than the minimum masked pixel
            % value, and higher than 20% of the peak intenisty value
            % (along the axial direction)
            z_int_max = max(pxl_int_profile, [], 1);
            int_lm_Q = pxl_int_profile >= int_lm_Q & ...
                pxl_int_profile > max(options.min_lm_int, ...
                (z_int_max * options.local_max_min_frac_max_int));
            [~, int_lm_idx] = max(int_lm_Q, [], 1);
            local_max_int = min(pxl_int_profile(...
                sub2ind(size(pxl_int_profile), int_lm_idx, 1:num_pxl)), options.max_peak_int);
            tmp_exp_th = (local_max_int - est_bg_med) * options.peak_int_fraction + est_bg_med;
            % Fractional max intensity point:
            tmp_exp_check_int_ab_th_Q = (pxl_int_profile > tmp_exp_th);
            hm_Q_diff = diff(cat(1, false(1, num_pxl), tmp_exp_check_int_ab_th_Q, ...
                false(1, num_pxl)));
            % Return the first greater than fractional maximum value position
            [max_Q_value, first_gt_fm_idx] = max(hm_Q_diff, [], 1);
            assert(all(max_Q_value==1));
            
            lm_int_sec_map = zeros(mask_yx_size);
            lm_int_sec_map(mask_pxl_ind) = int_lm_idx;
            % z-index of the surface in the mask 
            surf_sec_map = zeros(mask_yx_size);
            surf_sec_map(mask_pxl_ind) = first_gt_fm_idx;
            % Smooth the map
            surf_sec_map = medfilt2(surf_sec_map);
            if options.sf_map_imc_disk_r_pxl > 0
                surf_sec_map = imclose(surf_sec_map, strel('disk', options.sf_map_imc_disk_r_pxl));
            end            
            tmp_surf_map = surf_sec_map;
            tmp_surf_map(tmp_surf_map==0) = inf;
            % Constructing the surface mask
            surf_mask = repmat(permute(1:num_sections, [3,1,2]), mask_yx_size);
            % @ge for including the surface; @gt for excluding the surface
            surf_mask = bsxfun(@ge, surf_mask, tmp_surf_map);
            shg_mask_surf = surf_mask & shg_mask;
            if options.abl_mask_min_cc_pxl_2d
                % Remove small component per layer
                shg_mask_surf = WBIMAcqPostProcessing.filter_sections(shg_mask_surf, ...
                    @bwareaopen, options.abl_mask_min_cc_pxl_2d);
            end
            surf_mask_in_range_mip = any(shg_mask_surf(:, :, sec_idx_list), 3);
            if options.visQ                
                fig_hdl = figure;
                num_sp = 4;
                ax_0 = subplot(2,2,1);
                imagesc(ax_0, surf_sec_map);
                ax_0.Title.String = 'Surface layer';
                cmap = colormap(ax_0);
                cmap(1, :) = [0,0,0];
                colormap(ax_0, cmap);
                colorbar(ax_0);
                ax_1 = subplot(2,2,2);
                imagesc(ax_1, lm_int_sec_map);
                colorbar(ax_1);
                ax_1.Title.String = 'Peak intensity layer';
                colormap(ax_1, cmap);
                ax_2 = subplot(2,2,3);
                imagesc(ax_2, surf_mask_in_range_mip);
                colorbar(ax_2);
                ax_2.Title.String = 'Ablation mask in range';
                colormap(ax_2, 'gray');
                ax_3 = subplot(2,2,4);
                th_map = zeros(size(surf_sec_map));
                th_map(mask_pxl_ind) = tmp_exp_th;
                imagesc(ax_3, th_map);
                ax_3.Title.String = 'Threshold map';
                colorbar(ax_3);
                [ax_0.DataAspectRatio, ax_1.DataAspectRatio, ...
                    ax_2.DataAspectRatio, ax_3.DataAspectRatio] = deal([1,1,1]);
                if nargout > 1
                    varargout{1} = fig_hdl;
                end
            end                        
        end
        
            %% Vessel
        function [smip_mask_yx_um_itp, varargout] = compute_ablation_mask_yx_um_itp_vsl(stitched_mip, ...
                local_tile_info, est_th)
            % Deprecated
            % Use FWHM to compute the mask for the vascular channel
            % ablation 
            if nargin < 3
                est_th = [];
            end
            voxel_size_um = local_tile_info.step_mip_pixel_yxz_um;
            % Parameters
            ds_pxl_size_um = 10;
            clean_disk_d_um = 40;
            bg_est_margin_width_pxl = 5;
            % Downsample the image stack            
            ds_fraction = ds_pxl_size_um ./ voxel_size_um(1:2);
            
            ds_smip_10um = WBIMAcqPostProcessing.smooth_and_downsample_sections_in_stack(...
                stitched_mip, ds_fraction);
            ds_mip = max(ds_smip_10um, [], 3);
            ds_margin_stat = WBIMAcqPostProcessing.estimate_image_background_on_the_edge(...
                ds_mip, bg_est_margin_width_pxl);
            % Mask for profile analysis 
            if isempty(est_th)
                est_th = ds_margin_stat.th_3_sigma;
            elseif ~isempty(ds_margin_stat.th_3_sigma)
                est_th = max(est_th, ds_margin_stat.th_3_sigma);                
            end
            ds_mask = ds_mip > est_th;
            % Remove small connected components: 
%             [ds_mask_c, ds_cc] = WBIMAcqPostProcessing.remove_small_cc_by_distance(...
%                 ds_mask, 0.1, 20);
            clean_disk_r_pxl = round(clean_disk_d_um ./ ds_pxl_size_um / 2);
            min_cc_size = round(pi * clean_disk_r_pxl^2 * 4);
            ds_mask_c = bwareaopen(ds_mask, min_cc_size);
            % Should not fill the hole here. Some pixels might
            % have value smaller than ds_margin_stat.mean - cannot find
            % zero-crossing points.             
            % Remove small connected componennts?         
%             ds_mask_c = imclose(ds_mask_c, strel('disk', clean_disk_r_pxl));
            [ds_smip_mask_yx, ds_smip_fwhm] = WBIMAcqPostProcessing.compute_z_int_profile_fractional_mask(ds_smip_10um, ...
                ds_mask_c, ds_margin_stat.mean);          
            ds_smip_mask_yx = imclose(ds_smip_mask_yx, strel('disk', clean_disk_r_pxl));
            % Construct interpolation object
            smip_mask_yx_um_itp = WBIMAcqPostProcessing.construct_ablation_mask_yx_um_itp(...
                ds_smip_mask_yx, [ds_pxl_size_um, ds_pxl_size_um, voxel_size_um(3)]);                      
            smip_mask_yx_um_itp.abs_fp_z_um_mm = local_tile_info.abs_fp_z_um_done;
            smip_mask_yx_um_itp.step_mip_abs_fp_z_um = local_tile_info.step_mip_abs_fp_z_um;

%             smip_mask_yx_um = imresize3(smip_mask_yx_um_itp, round(size(stitched_mip) .* resize_voxel_size),...
%                 'Method', 'nearest');
            if nargout > 1
                ds_smip_fwhm.pxiel_size_um = ds_pxl_size_um;
                varargout{1} = ds_smip_fwhm;
            end
        end
        
        function vsl_smip_mask_ht = compute_vsl_mask_by_thresholding(smip_vsl, smip_shg, options)
            % Assume smip_vsl has the same size as smip_shg
            arguments
                smip_vsl
                smip_shg = []
                options.vsl2shg_int_min_fraction (1,1) double = 1;
                options.threshold (1,2) double = [2e3, 1e4];
                options.min_cc_vxl (1,1) double = 8;
                options.dilation_radius_pxl (1,1) double = 6; % 5 um resolution -> 30 um spacing 
                options.cmaxQ (1,1) logical = true;
            end
            vsl_smip_mask_ht = hysterisis_thresholding(smip_vsl, options.threshold(1), ...
                options.threshold(2));
            if ~isempty(smip_shg)
                vsl_smip_mask_ht = vsl_smip_mask_ht & (smip_vsl > smip_shg .* options.vsl2shg_int_min_fraction);
            end
            if options.min_cc_vxl > 0
                vsl_smip_mask_ht = bwareaopen(vsl_smip_mask_ht, options.min_cc_vxl);
            end
            vsl_smip_mask_ht = imfill(vsl_smip_mask_ht, 'holes');
            % Assume the voxel size of both images to be isotropic -
            % for using imdilate below
            if options.dilation_radius_pxl > 0
                vsl_smip_mask_ht = imdilate(vsl_smip_mask_ht, strel('sphere', options.dilation_radius_pxl));
            end
            if options.cmaxQ
                vsl_smip_mask_ht = cummax(vsl_smip_mask_ht, 3);
            end            
        end
        
        function [result, varargout] = estimate_surface_by_thresholding(smip_vsl, smip_shg, options)
            arguments
                smip_vsl
                smip_shg = []
                options.sec_range (1,2) double = [1, size(smip_vsl, 3)];
                options.threshold (1,:) double = [1e3, 3e3]
                options.vsl2shg_int_min_fraction (1,1) double = 1;
                options.min_cc_vxl (1,1) double = 8;
                options.dilation_radius_pxl (1,1) double = 3; % 15 um 
                options.cmaxQ (1,1) logical = true;
                options.visQ (1,1) logical = false;
            end
            smip_vsl(isnan(smip_vsl)) = 0;
            smip_shg(isnan(smip_shg)) = 0;
            smip_mask = WBIMAcqPostProcessing.compute_vsl_mask_by_thresholding(...
                smip_vsl, smip_shg, 'threshold', options.threshold, ...
                'vsl2shg_int_min_fraction', options.vsl2shg_int_min_fraction, ...
                'min_cc_vxl', options.min_cc_vxl, 'dilation_radius_pxl', options.dilation_radius_pxl, ...
                'cmaxQ', options.cmaxQ);
            
            result = struct;
            result.mask = smip_mask;
            [tmp_max_val, result.surface_sec_map] = max(smip_mask, [], 3);
            result.surface_sec_map = result.surface_sec_map .* tmp_max_val;
            result.surface_in_range_map = result.surface_sec_map >= options.sec_range(1) & ...
                result.surface_sec_map <= options.sec_range(2);
            result.abl_mask = result.mask;
            if options.visQ
                fig_hdl = figure;
                num_sp = 2;
                ax_0 = subplot(1,num_sp,1);
                imagesc(ax_0, result.surface_sec_map);
                ax_0.Title.String = 'Surface layer';
                cmap = colormap(ax_0);
                cmap(1, :) = [0,0,0];
                colormap(ax_0, cmap);
                colorbar(ax_0);
%                 ax_1 = subplot(1,3,2);
%                 imagesc(ax_1, result.lm_sec_map);
%                 colorbar(ax_1);
%                 ax_1.Title.String = 'Peak intensity layer';
%                 colormap(ax_1, cmap);
                ax_2 = subplot(1,num_sp,2);
                imagesc(ax_2, result.surface_in_range_map);
                colorbar(ax_2);
                ax_2.Title.String = 'Surface in range';
                colormap(ax_2, 'gray');
                [ax_0.DataAspectRatio, ax_2.DataAspectRatio] = deal([1,1,1]);
                if nargout > 1
                    varargout{1} = fig_hdl;
                end
            end            
        end
        
            %% SHG
        function [smip_mask_yx_um_itp, varargout] = compute_ablation_mask_yx_um_itp_SHG(stitched_mip, ...
                local_tile_info, est_th, est_fwhm_max)
            arguments
               stitched_mip
               local_tile_info
               est_th = []
               est_fwhm_max = []
            end
            voxel_size_um = local_tile_info.step_mip_pixel_yxz_um;
            
            % Downsample the image stack
            ds_pxl_size_um = 5;
            ds_fraction = ds_pxl_size_um ./ voxel_size_um(1:2);
            clean_disk_d_um = 40;
            clean_disk_r_pxl = min(round(clean_disk_d_um ./ ds_pxl_size_um / 2));
            
            ds_smip_5um = WBIMAcqPostProcessing.smooth_and_downsample_sections_in_stack(...
                stitched_mip, ds_fraction);
            ds_mip = max(ds_smip_5um, [], 3);
            ds_margin_stat = WBIMAcqPostProcessing.estimate_image_background_on_the_edge(...
                ds_mip, 5);
            % Mask for profile analysis 
            if isempty(est_th)
                est_th = ds_margin_stat.th_3_sigma;
            elseif ~isempty(ds_margin_stat.th_3_sigma)
                est_th = max(min(est_th, ds_margin_stat.th_5_sigma), ...
                    ds_margin_stat.th_3_sigma);
            end
            ds_mask = ds_mip > est_th;
            [ds_mask_c, ds_cc] = WBIMAcqPostProcessing.remove_small_cc_by_distance(...
                ds_mask, 0.1, 20);
            % Should not fill the hole here. Some pixels might
            % have value smaller than ds_margin_stat.mean - cannot find
            % zero-crossing points.             
            if ~isempty(est_fwhm_max)
                % Cap the maximum intensity for FWHM detection
                validateattributes(est_fwhm_max, {'numeric'}, {'nonnegative', 'scalar'});
                ds_smip_5um = max(ds_smip_5um, est_fwhm_max);
            end
            ds_smip_fwhm = WBIMAcqPostProcessing.compute_z_int_profile_fractional_mask(ds_smip_5um, ...
                ds_mask_c, 'est_bg', ds_margin_stat.mean);            
            % Remove small connected componennts?             
            ds_smip_fwhm.mask = imclose(ds_smip_fwhm.mask, strel('disk', clean_disk_r_pxl));
            
            smip_mask_yx_um_itp = WBIMAcqPostProcessing.construct_ablation_mask_yx_um_itp(...
                ds_smip_fwhm.mask, [ds_pxl_size_um, ds_pxl_size_um, voxel_size_um(3)]);
            smip_mask_yx_um_itp.abs_fp_z_um_mm = local_tile_info.abs_fp_z_um_done;
            smip_mask_yx_um_itp.step_mip_abs_fp_z_um = local_tile_info.step_mip_abs_fp_z_um;
            % Resize the mask back to 1 um resolution 
%             resize_voxel_size = voxel_size_um;
%             resize_voxel_size(3) = 1;
%             smip_mask_yx_um = imresize3(smip_mask_yx_um, round(size(stitched_mip) .* resize_voxel_size),...
%                 'Method', 'nearest');
            if nargout > 1
                ds_smip_fwhm.pxiel_size_um = ds_pxl_size_um;
                varargout{1} = ds_smip_fwhm;
            end
        end        
        
        function [shg_mask] = compute_shg_mask_by_thresholding(smip_shg, smip_vsl, options)
            arguments
                smip_shg 
                smip_vsl = []
                options.threshold (1,:) double = [3e3, 1e4];
                options.min_cc_vxl (1,1) double = 8;
                options.max_2d_hole_pxl (1,1) double = 9;
                options.bwclose_sph_r_pxl (1,1) double = 2; % 25 um diameter
            end
            % Deal with cross talk later. Need to check with the new
            % filters
            if isscalar(options.threshold)
                shg_mask = smip_shg > options.threshold;
            elseif numel(options.threshold) == 2
                shg_mask = hysterisis_thresholding(smip_shg, options.threshold(1), ...
                    options.threshold(2));
            end
            if options.max_2d_hole_pxl > 0
               shg_mask = bw_fill_small_holes_in_sections(shg_mask, options.max_2d_hole_pxl); 
            end            
            if options.bwclose_sph_r_pxl > 0
                shg_mask = imclose(shg_mask, strel('sphere', options.bwclose_sph_r_pxl));
            end            
            if options.min_cc_vxl > 0
                shg_mask = bwareaopen(shg_mask, options.min_cc_vxl);
            end
        end
        
        function [result, varargout] = estimate_surface_by_frac_z_peak_int(smip_shg, options)
            arguments
                smip_shg
                options.sec_range (1,2) double = [1, size(smip_shg, 3)];
                options.max_peak_int (1,1) double = inf;
                options.peak_int_fraction (1,1) double = 0.5;
                % Threshold 
                options.threshold (1,2) double = [7e2, 3e3];
                options.smip_mask_min_cc_size_pxl (1,:) double = 27; % in 3D mask
                % Local maximum detection
                options.local_max_wd_sz (1,1) double = 3
                options.local_max_min_frac_max_int (1,1) double = 0.2;
                % Smooth surface map
                options.sf_map_imc_disk_r_pxl (1,1) double = 2;
                % Clean up ablation mask 
                options.abl_mask_min_cc_pxl_2d (1,1) double = 0;
                options.visQ (1,1) logical = false;
            end
            smip_shg(isnan(smip_shg)) = 0;
            smip_mask = WBIMAcqPostProcessing.compute_shg_mask_by_thresholding(...
                smip_shg, [], 'threshold', options.threshold, ...
                'min_cc_vxl', options.smip_mask_min_cc_size_pxl, 'bwclose_sph_r_pxl', 0);            
            % Reshape array to get the z intensity profiles
            mask_yx_size = size(smip_shg, [1,2]);
            num_sections = size(smip_shg, 3);
            sec_idx_list = options.sec_range(1) : options.sec_range(2);
            search_map = any(smip_mask(:, :, sec_idx_list), 3);
            mask_pxl_ind = find(search_map);            
            num_pxl = numel(mask_pxl_ind);
            % Get masked pixel z intensity profile 
            pxl_int_profile = reshape(permute(smip_shg, [3,1,2]), num_sections, []);
            pxl_int_profile = pxl_int_profile(:, mask_pxl_ind);
            % Search for the local maxima in the intensity profile 
            int_lm_Q = movmax(pxl_int_profile, options.local_max_wd_sz, 1);
            % Local maximum should be higher than the minimum masked pixel
            % value, and higher than 20% of the peak intenisty value
            % (along the axial direction)
            z_int_max = max(pxl_int_profile, [], 1);
            int_lm_Q = pxl_int_profile >= int_lm_Q & ...
                pxl_int_profile > max(options.threshold(1), (z_int_max * options.local_max_min_frac_max_int));
            [~, int_lm_idx] = max(int_lm_Q, [], 1);
            local_max_int = min(pxl_int_profile(...
                sub2ind(size(pxl_int_profile), int_lm_idx, 1:num_pxl)), options.max_peak_int);
            % Estimate background level. Use median or mean?
            est_bg_med = median(smip_shg(~smip_mask), 'omitnan');
            tmp_exp_th = (local_max_int - est_bg_med) * options.peak_int_fraction + est_bg_med;
            % Fractional max intensity point:
            tmp_exp_check_int_ab_th_Q = (pxl_int_profile > tmp_exp_th);
            hm_Q_diff = diff(cat(1, false(1, num_pxl), tmp_exp_check_int_ab_th_Q, ...
                false(1, num_pxl)));
            % Return the first greater than fractional maximum value position
            [max_Q_value, first_gt_fm_idx] = max(hm_Q_diff, [], 1);
            assert(all(max_Q_value==1));
            % Construct the output 
            result = struct;
            result.mask = smip_mask;
            result.pxl_ind = mask_pxl_ind;
            result.peak_int = z_int_max;
            result.local_max_int = local_max_int;            
            lm_int_sec_map = zeros(mask_yx_size);
            lm_int_sec_map(mask_pxl_ind) = int_lm_idx;
            result.lm_sec_map = lm_int_sec_map;
            % z-index of the surface in the mask 
            surf_sec_map = zeros(mask_yx_size);
            surf_sec_map(mask_pxl_ind) = first_gt_fm_idx;
            % Smooth the map
            surf_sec_map = medfilt2(surf_sec_map);
            if options.sf_map_imc_disk_r_pxl > 0
                surf_sec_map = imclose(surf_sec_map, strel('disk', options.sf_map_imc_disk_r_pxl));
            end            
            result.surface_sec_map = surf_sec_map;
            tmp_surf_map = surf_sec_map;
            tmp_surf_map(tmp_surf_map==0) = inf;
            surf_mask = repmat(permute(1:num_sections, [3,1,2]), mask_yx_size);
            % @ge for including the surface; @gt for excluding the surface
            surf_mask = bsxfun(@ge, surf_mask, tmp_surf_map);
            
            result.surface_mask = surf_mask;
            result.abl_mask = result.surface_mask & result.mask;         
            if options.abl_mask_min_cc_pxl_2d
                % Remove small component per layer
                result.abl_mask = bw_area_open_in_sections(result.abl_mask, options.abl_mask_min_cc_pxl_2d);
            end
            result.abl_mask_in_range_mip = any(result.abl_mask(:, :, sec_idx_list), 3);
            if options.visQ
                fig_hdl = figure;
                num_sp = 4;
                ax_0 = subplot(2,2,1);
                imagesc(ax_0, result.surface_sec_map);
                ax_0.Title.String = 'Surface layer';
                cmap = colormap(ax_0);
                cmap(1, :) = [0,0,0];
                colormap(ax_0, cmap);
                colorbar(ax_0);
                ax_1 = subplot(2,2,2);
                imagesc(ax_1, result.lm_sec_map);
                colorbar(ax_1);
                ax_1.Title.String = 'Peak intensity layer';
                colormap(ax_1, cmap);
                ax_2 = subplot(2,2,3);
                imagesc(ax_2, result.abl_mask_in_range_mip);
                colorbar(ax_2);
                ax_2.Title.String = 'Ablation mask in range';
                colormap(ax_2, 'gray');
                ax_3 = subplot(2,2,4);
                th_map = zeros(size(result.surface_sec_map));
                th_map(mask_pxl_ind) = tmp_exp_th;
                imagesc(ax_3, th_map);
                ax_3.Title.String = 'Threshold map';
                colorbar(ax_3);
                [ax_0.DataAspectRatio, ax_1.DataAspectRatio, ...
                    ax_2.DataAspectRatio, ax_3.DataAspectRatio] = deal([1,1,1]);
                if nargout > 1
                    varargout{1} = fig_hdl;
                end
            end            
        end

        function [shg_mask_ts, varargout] = mask_refinement_shg_thin_shell(s_smip_shg, shg_mask, options)
            arguments
                s_smip_shg 
                shg_mask logical
                options.min_dI (1,1) double = 250;
                options.max_thin_shell_ls_um (1,1) double = 8;
                options.axial_pixel_size_um (1,1) double = 5;
                options.max_refine_int (1,1) double = 8e3;
                % Only consider normalized gradient with original intensity greater than this value
                options.min_search_int (1,1) double = 1e3; 
                options.max_gradient_offset (1,1) double = 1;
                options.mask_sm_r_pxl (1,1) double = 2;
                options.visQ (1,1) logical = false;
            end
            % Derived parameter
            min_neg_dI_n = options.axial_pixel_size_um / options.max_thin_shell_ls_um;
            
            s_smip_size = size(s_smip_shg);
            num_pxl = prod(s_smip_size(1:2));
            num_z = s_smip_size(3);
            shg_mask_v = reshape(permute(shg_mask, [3,1,2]), num_z, []);         
            % Compute gradient and normalized gradient along the z
            % direction 
            [dI_n, dI] = WBIMAcqPostProcessing.compute_z_normalized_derivative(s_smip_shg);
            % Remove positive gradient and noise
            is_invalid_Q = dI_n > 0 | s_smip_shg < options.min_search_int | abs(dI) < options.min_dI;
            dI_n(is_invalid_Q) = nan;
            [max_neg_dI_n, max_idx] = max(-dI_n, [], 3);        
            % Refine masks: 
            %   Remove mask pixels below the maximum normalized gradient
            %   point, down to the next mask starting point (if the 1D mask
            %   has more than 1 ablation intervals)            
            shg_mask_2d = any(shg_mask, 3);
            s_smip_shg_v = reshape(permute(s_smip_shg, [3,1,2]), num_z, []);
            shg_mask_v_edge = diff(cat(1, false(1, num_pxl), shg_mask_v, false(1, num_pxl)));            
            thin_shell_mask_2d = shg_mask_2d & (max_neg_dI_n > min_neg_dI_n);
            thin_shell_ind = find(thin_shell_mask_2d);
            thin_shell_mask_3d = false(size(shg_mask_v));
            for i = 1 : numel(thin_shell_ind)
                tmp_ind = thin_shell_ind(i);
                tmp_max_g_depth_idx = max_idx(tmp_ind);
                for j = (tmp_max_g_depth_idx + options.max_gradient_offset) : num_z
                    if s_smip_shg_v(j, tmp_ind) <= options.max_refine_int
                        thin_shell_mask_3d(j, tmp_ind) = true;
                    end
                    if shg_mask_v_edge(j+1, tmp_ind) == 1
                        break
                    end
                end
            end
            thin_shell_mask_3d = permute(reshape(thin_shell_mask_3d, s_smip_size([3,1,2])), [2,3,1]);
            if options.mask_sm_r_pxl 
                thin_shell_mask_3d = imclose(thin_shell_mask_3d, ...
                    strel('disk', options.mask_sm_r_pxl));
            end
            shg_mask_ts = shg_mask & (~thin_shell_mask_3d);
            
            if options.visQ
                max_ls_um = options.axial_pixel_size_um ./max_neg_dI_n;
                max_idx(isnan(max_ls_um)) = nan;
                max_grad_depth_um = max_idx * options.axial_pixel_size_um;
                
                fig_hdl = figure;
                ax_1 = nexttile();
                imagesc(ax_1, max_ls_um);
                ax_1.DataAspectRatio = [1,1,1];
                c_hdl_1 = colorbar(ax_1);
                ax_1.Title.String = 'Minimum scattering legth (\mum)';
                ax_1.CLim = [0, 50];
                ax_2 = nexttile();
                imagesc(ax_2, max_grad_depth_um);
                ax_2.DataAspectRatio = [1,1,1];
                c_hdl_2 = colorbar(ax_2);
%                 ax_2.CLim = [0, 150];
                ax_2.Title.String = 'Maximum gradient depth (\mum)';
%                 [ax_1.XLabel.String, ax_2.XLabel.String] = deal(sprintf('X (%d \\mum)',  ...
%                     s_region_info.step_mip_pixel_yxz_um(2)));
%                 [ax_1.YLabel.String, ax_2.YLabel.String] = deal(sprintf('Y (%d \\mum)',  ...
%                     s_region_info.step_mip_pixel_yxz_um(1)));
                ax_3 = nexttile();
                imagesc(ax_3, any(thin_shell_mask_3d, 3));
                ax_3.DataAspectRatio = [1,1,1];
                ax_3.Title.String = sprintf('l_s < %.1f \\mum mask', options.max_thin_shell_ls_um);
                ax_4 = nexttile();
                h_bin = 0 : 1 : 50;
                h_hdl = histogram(ax_4, max_ls_um, 'BinEdges', h_bin);
                ax_4.XLabel.String = 'Minimum scattering length (\mum)';
                ax_4.YLabel.String = 'Counts';       
            else
                fig_hdl = [];
            end
            if nargout > 1
                varargout{1} = fig_hdl;
            end
        end
        

        %% Combined
        function result_str = compute_ablation_mask_yx_um_itp_vNs(mip_vsl, mip_shg, local_tile_info, options)
            arguments
                mip_vsl 
                mip_shg
                local_tile_info (1,1) struct
                options.est_vsl_th = [1e3, 3e3]
                options.est_shg_th = [7e2, 2e3]
                options.rm_shg_thin_shell_Q (1,1) logical = true;
                options.fill_hole_max_num_pxl (1,1) double = 50
                options.cmax_vslQ (1,1) logical = true;
                options.crmax_shgQ (1,1) logical = false;
                
                options.mask_spec_unmix_vsl_Q = true;
                options.mask_spec_unmix_shg_Q = true;
                % False positive detection will only waste time; False
                % negative exclusion will accumulate error out of control.
                % - Do not exclude
                options.exclude_shg_mask_from_vsl_mask_Q (1,1) logical = false;            
                options.min_cc_size_pxl (1,:) double = 27;
                options.visQ (1,1) logical = true;
            end
            [vsl_mask, shg_mask] = WBIMAcqPostProcessing.compute_ablation_mask_yx_um_vNs(...
                mip_vsl, mip_shg, 'est_vsl_th', options.est_vsl_th, ...
                'est_shg_th', options.est_shg_th, 'rm_shg_thin_shell_Q', options.rm_shg_thin_shell_Q, ...
                'fill_hole_max_num_pxl', options.fill_hole_max_num_pxl, ...
                'cmax_vslQ', options.cmax_vslQ, 'crmax_shgQ', options.crmax_shgQ, ...
                'mask_spec_unmix_shg_Q', options.mask_spec_unmix_shg_Q, ...
                'mask_spec_unmix_vsl_Q', options.mask_spec_unmix_vsl_Q, ...
                'exclude_shg_mask_from_vsl_mask_Q', options.exclude_shg_mask_from_vsl_mask_Q, ...
                'min_cc_size_pxl', options.min_cc_size_pxl, ...
                'visQ', options.visQ);
            % Iterpolation
            result_str = struct;
            result_str.region_bbox_ctr_yx_um = local_tile_info.region_bbox_ctr_yx_um;
            % Add in reversed order            
            result_str.itp(WBIMChannelName.SHG) = WBIMAcqPostProcessing.construct_ablation_mask_yx_um_itp(...
                shg_mask, local_tile_info.step_mip_pixel_yxz_um);
            
            result_str.itp(WBIMChannelName.Vessel) = WBIMAcqPostProcessing.construct_ablation_mask_yx_um_itp(...
                vsl_mask, local_tile_info.step_mip_pixel_yxz_um);
            
            result_str.abs_fp_z_um_mm = local_tile_info.abs_fp_z_um_done;
            result_str.step_mip_abs_fp_z_um = local_tile_info.step_mip_abs_fp_z_um;
            % For debug purpose
            result_str.mask = {vsl_mask, shg_mask};
        end        
        
        function [vsl_mask, shg_mask] = compute_ablation_mask_yx_um_vNs(mip_vsl, mip_shg, options)
            arguments
                mip_vsl 
                mip_shg
                options.est_vsl_th = [1e3, 3e3]
                options.est_shg_th = [7e2, 2e3]
                % Skull thin shell detection 
                options.rm_shg_thin_shell_Q (1,1) logical = true;                
                
                options.fill_hole_max_num_pxl (1,1) double = 50;
                options.cmax_vslQ (1,1) logical = true;
                options.crmax_shgQ (1,1) logical = false;
                
                options.mask_spec_unmix_vsl_Q = true;
                options.mask_spec_unmix_shg_Q = true;
                % False positive detection will only waste time; False
                % negative exclusion will accumulate error out of control.
                % - Do not exclude
                options.exclude_shg_mask_from_vsl_mask_Q (1,1) logical = false;               
                options.min_cc_size_pxl (1,:) double = 27;
                options.visQ (1,1) logical = false;
            end
            % Compute mip mask array x            
            vsl_mask = WBIMAcqPostProcessing.compute_vsl_mask_by_thresholding(...
                mip_vsl, mip_shg, 'threshold', options.est_vsl_th, ...
                'dilation_radius_pxl', 1);
            
            shg_mask = WBIMAcqPostProcessing.compute_shg_mask_by_thresholding(...
                mip_shg, mip_vsl, 'threshold', options.est_shg_th);
            
            if options.rm_shg_thin_shell_Q
               shg_mask = WBIMAcqPostProcessing.mask_refinement_shg_thin_shell(...
                   mip_shg, shg_mask, 'min_search_int', options.est_shg_th(1), ...
                   'visQ', false);
            end            
            % Fill small holes 
            if isfinite(options.fill_hole_max_num_pxl) && options.fill_hole_max_num_pxl > 0
                shg_mask = WBIMAcqPostProcessing.filter_sections(shg_mask,...
                    @bw_fill_small_holes, options.fill_hole_max_num_pxl);                
            end
            % Fill all the holes in the tissue mask 
            vsl_mask = WBIMAcqPostProcessing.filter_sections(vsl_mask, ...
                @bw_fill_small_holes_in_sections, inf);
            % Combine masks - remove vessels from the SHG channel 
            if options.mask_spec_unmix_vsl_Q || options.mask_spec_unmix_shg_Q                
                [shg_in_vsl_mask, vsl_in_shg_mask] = WBIMAcqPostProcessing.compute_step_mip_spectrum_unmixing_mask(...
                    mip_vsl, mip_shg);
                if options.mask_spec_unmix_vsl_Q
                    vsl_mask = vsl_mask & ~shg_in_vsl_mask;
                end                
                if options.mask_spec_unmix_shg_Q
                    shg_mask = shg_mask & ~vsl_in_shg_mask;
                end
            end
            % Remove small connected components. Operate in 3D mask
            % directly.
            if ~isempty(options.min_cc_size_pxl) && isfinite(options.min_cc_size_pxl)
                vsl_mask = bwareaopen(vsl_mask, options.min_cc_size_pxl);
                shg_mask = bwareaopen(shg_mask, options.min_cc_size_pxl);                
            end
            % Cumulated max projection for the vascular channel
            if options.cmax_vslQ
                % Potential problem - vessels / sinus in the skull ?
                % Assume there's always tissues underneath the tissues.
                % Might be less efficient when reaching the bottom of the
                % brain. 
                vsl_mask = cummax(vsl_mask, 3);
            end
            % Reverse cumulated max projection for the SHG channel 
            if options.crmax_shgQ
               shg_mask = cummax(shg_mask, 3, 'reverse');
            end
            
            if options.exclude_shg_mask_from_vsl_mask_Q
                vsl_mask = vsl_mask & ~shg_mask;
            end            
            if options.visQ
                vsl_depth_map = uint8(sum(vsl_mask, 3));
                shg_depth_map = uint8(sum(shg_mask, 3));
                f1 = WBIMAcqPostProcessing.vis_abl_thickness_map(vsl_depth_map);
                f2 = WBIMAcqPostProcessing.vis_abl_thickness_map(shg_depth_map);                
            end            
        end
        
        function result_str = compute_ablation_mask_yx_um_itp_vNs_FWHM(mip_vsl, mip_shg, local_tile_info, options)
            % Deprecated
            arguments
                mip_vsl 
                mip_shg
                local_tile_info (1,1) struct
                options.est_vsl_th = 4000
                options.est_shg_th = 4000
                
                options.fill_hole_max_num_pxl (1,1) double = 1256
                options.cmax_vslQ (1,1) logical = true;
                options.crmax_shgQ (1,1) logical = false;
                
                options.mask_spec_unmix_vsl_Q = true;
                options.mask_spec_unmix_shg_Q = true;
                % False positive detection will only waste time; False
                % negative exclusion will accumulate error out of control.
                % - Do not exclude
                options.exclude_shg_mask_from_vsl_mask_Q (1,1) logical = false; 
                options.min_cc_size_pxl (1,:) double = 100;
                options.visQ (1,1) logical = true;
            end
            % Compute mip mask array 
            abl_str_vsl = WBIMAcqPostProcessing.compute_step_mip_fractional_max_mask(...
                mip_vsl, 0.5, 'est_bg_th', options.est_vsl_th, 'visQ', false);
            vsl_mask = abl_str_vsl.mask;
            
            abl_str_shg = WBIMAcqPostProcessing.compute_step_mip_fractional_max_mask(...
                mip_shg, 0.25, 'est_bg_th', options.est_shg_th, 'visQ', false, ...
                'remove_distant_small_cc_Q', true);
            shg_mask = abl_str_shg.mask;
            % Fill small holes 
            if isfinite(options.fill_hole_max_num_pxl) && ...
                    options.fill_hole_max_num_pxl > 0
                shg_mask = bw_fill_small_holes_in_sections(shg_mask, options.fill_hole_max_num_pxl);
                % Fill all the holes in the tissue mask 
                vsl_mask = bw_fill_small_holes_in_sections(vsl_mask, inf);
            end
            % Combine masks - remove vessels from the SHG channel 
            
            if options.mask_spec_unmix_vsl_Q || options.mask_spec_unmix_shg_Q                
                [shg_in_vsl_mask, vsl_in_shg_mask] = WBIMAcqPostProcessing.compute_step_mip_spectrum_unmixing_mask(...
                    mip_vsl, mip_shg);
                if options.mask_spec_unmix_vsl_Q
                    vsl_mask = vsl_mask & ~shg_in_vsl_mask;
                end                
                if options.mask_spec_unmix_shg_Q
                    shg_mask = shg_mask & ~vsl_in_shg_mask;
                end
            end
            if ~isempty(options.min_cc_size_pxl) && isfinite(options.min_cc_size_pxl)
                vsl_mask = bwareaopen(vsl_mask, options.min_cc_size_pxl);
                shg_mask = bwareaopen(shg_mask, options.min_cc_size_pxl);                
            end
            % Remove small connected components? 
            if options.cmax_vslQ
                % Cumulated max projection for the vascular channel
                % Potential problem - vessels / sinus in the skull ?
                % Assume there's always tissues underneath the tissues.
                % Might be less efficient when reaching the bottom of the
                % brain. 
                vsl_mask = cummax(vsl_mask, 3);
            end
            if options.crmax_shgQ
               % Reverse cumulated max projection for the SHG channel 
               shg_mask = cummax(shg_mask, 3, 'reverse');
            end
            
            if options.exclude_shg_mask_from_vsl_mask_Q
                vsl_mask = vsl_mask & ~shg_mask;
            end
            % Add the second layer of SHG mask to the first layer
            shg_mask(:, :, 1) = shg_mask(:, :, 1) | shg_mask(:, :, 2);
            
            if options.visQ
                vsl_depth_map = uint8(sum(vsl_mask, 3));
                shg_depth_map = uint8(sum(shg_mask, 3));
                f1 = WBIMAcqPostProcessing.vis_abl_thickness_map(vsl_depth_map);
                f2 = WBIMAcqPostProcessing.vis_abl_thickness_map(shg_depth_map);                
            end            
            % Iterpolation
            result_str = struct;
            % Add in reversed order            
            result_str.itp(WBIMChannelName.SHG) = WBIMAcqPostProcessing.construct_ablation_mask_yx_um_itp(...
                shg_mask, local_tile_info.step_mip_pixel_yxz_um);
            
            result_str.itp(WBIMChannelName.Vessel) = WBIMAcqPostProcessing.construct_ablation_mask_yx_um_itp(...
                vsl_mask, local_tile_info.step_mip_pixel_yxz_um);
            
            result_str.abs_fp_z_um_mm = local_tile_info.abs_fp_z_um_done;
            result_str.step_mip_abs_fp_z_um = local_tile_info.step_mip_abs_fp_z_um;
            % For debug purpose
            result_str.mask = {vsl_mask, shg_mask};
        end        
        
        function result_str = analyze_explore_tile_smip(smip_cell,...
                local_tile_info, options)
            arguments
                smip_cell 
                local_tile_info (1,1) struct
                options.est_vsl_th = [1e3, 3e3]
                options.est_shg_th = [7e2, 2e3]
                
                options.surf2peak_int_ratio (1,1) double = 0.3;
                options.tissue_peak_int = 7e3;
                options.sec_range (1, :) double = [1, size(smip_cell{1}, 3)];
                % Reserve...
                options.mask_spec_unmix_vsl_Q = true;
                options.mask_spec_unmix_shg_Q = true;
                options.vsl_min_cc_pxl_2d (1,1) double = 16;
                options.visQ (1,1) logical = false;
            end
            num_ch = numel(smip_cell);
            abl_vol_cell = cell(num_ch, 1);
            for i_ch = 1 : num_ch
                if ~isempty(smip_cell{i_ch})
                    switch i_ch
                        case WBIMChannelName.Vessel
                            % Compute mip mask array
                            abl_vol_cell{i_ch} = WBIMAcqPostProcessing.estimate_surface_by_frac_z_peak_int(...
                                smip_cell{i_ch}, 'sec_range', options.sec_range, ...
                                'threshold', options.est_vsl_th, 'peak_int_fraction', options.surf2peak_int_ratio, ...
                                'abl_mask_min_cc_pxl_2d', options.vsl_min_cc_pxl_2d,...
                                'visQ', options.visQ , 'max_peak_int', options.tissue_peak_int);
                        case WBIMChannelName.SHG
                            abl_vol_cell{i_ch} = WBIMAcqPostProcessing.estimate_surface_by_frac_z_peak_int(...
                                smip_cell{i_ch}, 'sec_range', options.sec_range, 'threshold', options.est_shg_th, ...
                                'peak_int_fraction', options.surf2peak_int_ratio,  ...
                                'visQ', options.visQ);
                    end
                end
            end
            
            % Combined reuslt 
            result_str = struct;
            
            result_str.abs_fp_z_um_mm = local_tile_info.abs_fp_z_um_done;
            result_str.step_mip_abs_fp_z_um = local_tile_info.step_mip_abs_fp_z_um;
            result_str.region_bbox_ctr_yx_um = local_tile_info.region_bbox_ctr_yx_um;
            result_str.step_mip_pixel_yxz_um = local_tile_info.step_mip_pixel_yxz_um;
            
            result_str.local_tile_info = local_tile_info;
            % Combine / process the SHG and VSL mask here.         
            
            % Last step before interpolation generation 
            if isfield(local_tile_info, 'scan_roi_smip_bbox_pxl')
               % Remove mask outside the scan ROI
               scan_roi_mask = false(local_tile_info.step_mip_size([1,2]));
               scan_roi_mmxx = local_tile_info.scan_roi_smip_bbox_pxl;
               % Shrink the bounding box? 
               scan_roi_mmxx(1:2) = scan_roi_mmxx(1:2) + 1;
               scan_roi_mmxx(3:4) = scan_roi_mmxx(3:4) - 1;
               
               scan_roi_mask(scan_roi_mmxx(1) : scan_roi_mmxx(3), ...
                   scan_roi_mmxx(2) : scan_roi_mmxx(4)) = true;
               
               result_str.scan_roi_mask = scan_roi_mask;
               for i_c = 1 : num_ch
                   tmp_result = abl_vol_cell{i_c};
                   if ~isempty(tmp_result)
                       tmp_result.abl_mask = logical(tmp_result.abl_mask .* scan_roi_mask);
                       tmp_result.abl_mask_in_range_mip = tmp_result.abl_mask_in_range_mip .* ...
                           scan_roi_mask;
                       abl_vol_cell{i_c} = tmp_result;
                   end
               end               
            end   
            result_str.need_ablation_Q = any(cellfun(@(x) any(x.abl_mask_in_range_mip, 'all'), ...
                abl_vol_cell));
            result_str.surf_str = abl_vol_cell;
            % Generate ablation ROI interpolation object
            if result_str.need_ablation_Q                
                % Add in reversed order
                for i_c = num_ch : -1 : 1
                    result_str.itp(i_c) = WBIMAcqPostProcessing.construct_ablation_mask_yx_um_itp(...
                        result_str.surf_str{i_c}.abl_mask, local_tile_info.step_mip_pixel_yxz_um);
                end
            end
            
            if options.visQ
%                 fig_hdl = figure;
%                 ax_0 = subplot(1,3,1);
%                 imagesc(ax_0, any(shg_result.abl_mask, 3));
%                 ax_0.Title.String = 'Ablation mask MIP';
%                 cmap = colormap(ax_0);
%                 cmap(1, :) = [0,0,0];
%                 colormap(ax_0, cmap);
%                 colorbar(ax_0);
%                 ax_1 = subplot(1,3,2);
%                 imagesc(ax_1, result.lm_sec_map);
%                 colorbar(ax_1);
%                 ax_1.Title.String = 'Peak intensity layer';
%                 colormap(ax_1, cmap);
%                 ax_2 = subplot(1,3,3);
%                 imagesc(ax_2, result.surface_in_range_map);
%                 colorbar(ax_2);
%                 ax_2.Title.String = 'Surface in range';
%                 colormap(ax_2, 'gray');
%                 [ax_0.DataAspectRatio, ax_1.DataAspectRatio, ax_2.DataAspectRatio] = deal([1,1,1]);
%                 if nargout > 1
%                     varargout{1} = fig_hdl;
%                 end
            end
        end
        
        function inQ = check_tiles_in_mask_2D(mask_yx_um, local_bbox_mmxx_um)
            arguments
               mask_yx_um (:, :) logical
               local_bbox_mmxx_um (:, 4) double
            end
            mask_size = size(mask_yx_um);
            local_bbox_mmxx_um = round(min(max([1,1,1,1], local_bbox_mmxx_um), ...
                [mask_size, mask_size]));
            num_bbox = size(local_bbox_mmxx_um, 1);
            inQ = false(num_bbox, 1);
            for i = 1 : num_bbox
                tmp_bbox_mmxx = local_bbox_mmxx_um(i, :);
                tmp_mask = mask_yx_um(tmp_bbox_mmxx(1) : tmp_bbox_mmxx(3), ...
                    tmp_bbox_mmxx(2) : tmp_bbox_mmxx(4));
                inQ(i) = any(tmp_mask, 'all');                
            end            
        end
        %% Visualization
        function fig_hdl = vis_abl_thickness_map(depth_map)
            if ~isempty(depth_map)
                fig_hdl = figure;
                ax_hdl = axes(fig_hdl);
                imagesc(ax_hdl, depth_map);
                ax_hdl.DataAspectRatio = [1,1,1];
                cbar = colorbar(ax_hdl);
                cbar.Label.String = 'Number of ablation layers';
            end
        end
        
        function fig_hdl = vis_step_mip_with_mask(stitched_mip, mask, section_list)
            if nargin < 3
                section_list = 1 : size(mask, 3);
            end
            fig_hdl = figure;
            tiledlayout('flow');
            for i = 1 : numel(section_list)
                tmp_idx = section_list(i);
                tmp_im = fun_stretch_contrast(stitched_mip(:, :, tmp_idx));
                tmp_mask = mask(:, :, tmp_idx);
                ax_hdl = nexttile;
%                 im_1 = imagesc(ax_hdl, tmp_im);
%                 ax_hdl.DataAspectRatio = [1,1,1];
%                 hold(ax_hdl, 'on');
%                 im_2 = imagesc(tmp_mask, 
                imshowpair(tmp_im, tmp_mask);
            end            
        end
    end
    %% Utilities
    methods(Static)
        function exit_code = split_SI_raw_tiff_channels(im_fp, active_channel)
            persistent DM
            if isempty(DM)
                DM = WBIMFileManager;
            end
            num_channels = numel(active_channel);
            if num_channels == 1
                exit_code = 0;
                return;
            else
                assert(num_channels > 0);
                fprintf('Start processing %s\n', im_fp);
                raw_im_hdl = ScanImageTiffReader.ScanImageTiffReader(im_fp);
                [folder_path, file_name, file_ext] = fileparts(im_fp);
                raw_im_stack = raw_im_hdl.data();
                raw_im_stack = permute(raw_im_stack, [2,1,3]);
                for i = 1 : num_channels
                    tmp_data = raw_im_stack(:, :, i:num_channels:end);
                    tmp_fp = fullfile(folder_path, sprintf('%s_ch%d.%s', ...
                        file_name, active_channel(i), file_ext));
                    DM.write_tiff_stack(tmp_data, tmp_fp);
                    fprintf('Finish writing %s\n', tmp_fp);
                end
            end
        end
        
        function im3d = filter_sections(im3d, fun_hdl, varargin)
            num_sec = size(im3d, 3);
            for i = 1 : num_sec
                im3d(:, :, i) = fun_hdl(im3d(:, :, i), varargin{:});
            end            
        end
    end
end