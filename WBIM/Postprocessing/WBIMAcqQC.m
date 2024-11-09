classdef WBIMAcqQC < handle
    % Acquisition image quanlity control
    properties
        est_row_shift_per_tile (1, 1) double = 2;
        est_shift (1, 1) double = 0
        low_int_count (1, 1) double = 0
        
    end
    
    methods
        function obj = WBIMAcqQC()
            
            
        end
        
        function add_new_observation(obj, stack_cells)
            arguments
                obj WBIMAcqQC
                stack_cells (1, :) cell
            end
            
            
            
        end
    end
        
        %% Lineshift detection
        methods(Static)
            function rs_info = detect_row_shift_in_scan_stack(stack_cells, opts)
                arguments
                    stack_cells (1, :) cell
                    
                    % Number of columns to remove on each side. To elimnate
                    % the lineshift
                    opts.num_rm_col (1,1) double = 10; 
                    opts.max_possible_rs (1, 1) double = 30;
                    opts.visQ (1, 1) logical = false;
                end
                num_ch = numel(stack_cells);
                rs_info = struct('est_shift', nan(1, num_ch), ...
                    'est_shift_in_range', nan(1, num_ch));
                for i_ch = 1 : num_ch
                    tmp_stack = stack_cells{i_ch};
                    if ~isempty(tmp_stack)
                        tmp_mip = max(tmp_stack, [], 3);
                        tmp_mip = tmp_mip(:, (opts.num_rm_col + 1) : (end - opts.num_rm_col));
                        tmp_rs_stat = WBIMAcqQC.detect_row_shift_in_single_image(...
                            tmp_mip, 'visQ', opts.visQ);
                        if ~isempty(tmp_rs_stat.num_row_shift)
                            rs_info.est_shift(i_ch) = tmp_rs_stat.num_row_shift;
                            if abs(tmp_rs_stat.num_row_shift) <= opts.max_possible_rs
                                rs_info.est_shift_in_range(i_ch) = tmp_rs_stat.num_row_shift;
                            end
                        end
                    end
                end                
            end
            
            
            
            function row_shift_stat = detect_row_shift_in_single_image(im, opts)
                arguments
                    im 
                    opts.min_fg_int (1,1) double = 5e3;
                    opts.ipr (1,1) double = 1.5;
                    opts.min_corr_th (1,1) double = 0.05;
                    opts.min_int_grad_th (1,1) double = 1500;
                    opts.min_int_grad2m_th (1,1) double = 3e2;
                    opts.cc_ms_wd (1,1) double = 9;
                    % number of lines correspond to 3 ms, based on resonant scanner norminal
                    % frequency
                    opts.max_rs_duration_s (1,1) double = 3e-3;
                    opts.visQ (1,1) logical = false;
                end
                opts.max_oi_length = round(opts.max_rs_duration_s * 7920 * 2);
                assert(ismatrix(im), 'The input image must be 2D');
                bg_Q = im < opts.min_fg_int;
                row_shift_stat = struct('num_row_shift', [], 'low_bg_Q', false);
                
                if ~all(bg_Q, 'all')
                    im(bg_Q) = 0;
                    row_stat = WBIMAcqQC.compute_row_stat_of_image(im);
                    % Detect outliers
                    outlier_sm_se = strel('rectangle', [opts.cc_ms_wd, 1]);
                    [corr_th_range] = WBIMAcqQC.compute_percentile_inlier_range(row_stat.adj_corr_da, opts.ipr);
                    corr_th = max(corr_th_range(2), opts.min_corr_th);
                    corr_outlier_Q = imclose(row_stat.adj_corr_da > corr_th, outlier_sm_se);
                    
                    [int_grad_range] = WBIMAcqQC.compute_percentile_inlier_range(row_stat.row_abs_diff, opts.ipr);
                    int_grad_th = max(opts.min_int_grad_th, int_grad_range(2));
                    int_grad_outlier_Q = imclose(row_stat.row_abs_diff > int_grad_th, outlier_sm_se);
                    
                    row_shift_stat.low_bg_Q = all(row_stat.row_abs_diff <= opts.min_int_grad_th, 'all');
                    
                    is_outlier_Q = imclose(corr_outlier_Q & int_grad_outlier_Q,...
                        outlier_sm_se);
                    
                    % Select outliers
                    if any(is_outlier_Q)
                        % If detected - double check
                        % Seems better than isoutlier with movmean median.
                        % Need a minimum threshold? 
                        abs_diff_sm = movmean(row_stat.row_abs_diff, opts.max_oi_length * 2, ...
                            'omitnan', 'Endpoints', 'shrink');
                        abs_diff_dm = abs(row_stat.row_abs_diff - abs_diff_sm);
                        [int_graddm_range] = WBIMAcqQC.compute_percentile_inlier_range(abs_diff_dm, opts.ipr);
                        int_graddm_th = max(int_graddm_range(2), opts.min_int_grad2m_th);
                        is_mw_outlier_Q = imclose(abs_diff_dm > int_graddm_th, outlier_sm_se);
                        is_outlier_Q = is_outlier_Q & is_mw_outlier_Q;                        
                    end                    

                    if any(is_outlier_Q)
                        % Find outlier endpoint indices
                        o_ep_ind = reshape(find(diff([false; is_outlier_Q; false])), 2, []).';
                        o_ep_ind(:, 2) = o_ep_ind(:, 2) - 1;
                        
                        o_length = o_ep_ind(:, 2) - o_ep_ind(:, 1) + 1;
                        o_ep_ind = o_ep_ind((o_length < opts.max_oi_length), :);
                        
                        num_cc = size(o_ep_ind, 1);
                        if num_cc > 0
                            if num_cc > 1
                                % Further select cc
                                o_adj_cad_sum = arrayfun(@(x) sum(row_stat.adj_corr_da(o_ep_ind(x, 1):o_ep_ind(x, 2))), 1:num_cc);
                                [~, max_idx] = max(o_adj_cad_sum);
                                o_length = o_length(max_idx);
                                o_ep_ind = o_ep_ind(max_idx, :);
                            end
                            num_rows = size(im, 1);
                            
                            row_shift_stat.ep_ind = o_ep_ind;
                            row_shift_stat.avg_row_int = mean(row_stat.row_mean(o_ep_ind(1) : o_ep_ind(2)));
                            
                            
                            
                            % Determine the row shift
                            o_med_ind = (o_ep_ind(1) + o_ep_ind(2)) / 2;
                            if o_med_ind < (num_rows - o_med_ind)
                                row_shift_stat.num_row_shift = o_ep_ind(2);
                            else
                                row_shift_stat.num_row_shift = o_ep_ind(1) - num_rows;
                            end
                        end
                    else
                        o_ep_ind = [];
                    end
                else
                    row_shift_stat.low_bg_Q = true;
                end
                
                if opts.visQ
                    fig_hdl = figure;
                    fig_hdl.Position(3) = fig_hdl.Position(3) * 2;
                    ax1 = subplot(1,2,1);
                    % Statistics
                    yyaxis(ax1, 'left');
                    plt1 = plot(ax1, row_stat.adj_corr_da);
                    hold(ax1, 'on');
                    yline(ax1, corr_th, 'Color', plt1.Color);
                    ax1.YLabel.String = 'AdjCorrDA';
                    yyaxis(ax1, 'right');
                    plt2 = plot(ax1, row_stat.row_abs_diff);
                    yline(ax1, int_grad_th, 'Color', plt2.Color);
%                     plt3 = plot(ax1, row_stat.row_abs_diff2, '*');
%                     yline(ax1, int_grad2_th, 'Color', plt3.Color);
                    ax1.YLabel.String = 'AvgGradInt';
                    if ~isempty(o_ep_ind)
                        xline(ax1, o_ep_ind(1), 'r');
                        xline(ax1, o_ep_ind(2), 'r');         
                        ax1.Title.String = sprintf('Row shift artefact: [%d, %d]', o_ep_ind);
                    end

                    ax2 = subplot(1,2,2);
                    imagesc(ax2, im);
                    colorbar(ax2);
                    ax2.DataAspectRatio = [1,1,1];
                    if ~isempty(o_ep_ind)
                        hold(ax2, 'on');
                        yline(ax2, o_ep_ind(1), 'r');
                        yline(ax2, o_ep_ind(2), 'r');
                    end
                end
            end
            
            function im = zero_voxel_below_threashold(im, th)
                arguments
                    im
                    th (1,1) double
                end
                im(im < th) = 0;
            end
            
            function row_stat = compute_row_stat_of_image(im)
                assert(ismatrix(im), 'The input should be a 2D image (matrix)');

                if ~isfloat(im)
                    im = single(im);
                end
                row_stat = struct;
                row_stat.row_mean = mean(im, 2);
                row2_mean = mean(im.^2, 2);
                row_std = sqrt(row2_mean - row_stat.row_mean .^2);
                row_std(row_std == 0) = nan;
                row_std = fillmissing(row_std, 'linear', 'EndValues', 'nearest');
                row_stat.std = row_std;
                % Adjacent row correlation coefficient 
                im = im - row_stat.row_mean;
                adj_corr = mean(im(1:end-1, :) .* im(2:end, :), 2) ./ ...
                    (row_std(1:end-1) .* row_std(2:end));
                adj_corr = cat(1, adj_corr(1), max(adj_corr(1:end-1), adj_corr(2:end)), ...
                    adj_corr(end));
                row_stat.adj_corr = adj_corr;
                row_stat.adj_corr_d = gradient(row_stat.adj_corr);
                row_stat.adj_corr_da = abs(row_stat.adj_corr_d);
                
                % Average adjacent row intensity gradient absolute difference
                im_grad = abs(im(1:end-1, :) - im(2:end, :));
                im_grad = cat(1, im_grad(1, :), (im_grad(1:end-1, :) + im_grad(2:end, :)) ./ 2, ...
                    im_grad(end, :));
                % Second order derivative along the column: 
%                 im_grad2_c = im(3:end, :) + im(1:end-2, :) - 2 .* im(2:end-1, :);
%                 im_grad2_1 = 2 .* im(1, :) - 5 .* im(2, :) + 4 .* im(3, :) - im(4, :);
%                 im_grad2_end = 2 .* im(end, :) - 5 .* im(end-1, :) + 4 .* im(end-2, :) - im(end-3, :);
%                 im_grad2 = abs(cat(1, im_grad2_1, im_grad2_c, im_grad2_end));                
                
                
                row_stat.row_abs_diff = mean(im_grad, 2);
%                 row_stat.row_abs_diff2 = mean(im_grad2, 2);
            end
            
            function [th_range, varargout] = compute_percentile_inlier_range(data, ipr)
                if nargin < 2
                    ipr = 1.5;
                end
                p = prctile(data(:), [25, 50, 75]);
                half_width = (p(3) - p(1)) * ipr;
                th_low = p(1) - half_width;
                th_high = p(3) + half_width;
                th_range = [th_low, th_high];
                if nargout > 1
                    varargout{1} = half_width;
                end
            end
            
            function cc = bwconncomp1d(mask_1d)
                arguments
                    mask_1d (:, 1) logical
                end
                o_ep_ind = reshape(find(diff([false; mask_1d; false])), 2, []).';
                o_ep_ind(:, 2) = o_ep_ind(:, 2) - 1;
                
                cc = struct('Connectivity', 2, 'ImageSize', size(mask_1d), ...
                    'NumObjects', size(o_ep_ind, 1), ...
                    'Size', o_ep_ind(:, 2) - o_ep_ind(:, 1) + 1);
                cc.PixelIdxList = arrayfun(@(x) (o_ep_ind(x, 1) : o_ep_ind(x, 2)).', ...
                    1:cc.NumObjects, 'UniformOutput', false);
                cc.PixelEndpointList = o_ep_ind;
            end
            
            function mask_1d = remove_large_cc_1d(mask_1d, max_size)
               arguments
                   mask_1d (:, 1) logical
                   max_size (1, 1) double {mustBeNonnegative}
               end
               cc_1d = WBIMAcqQC.bwconncomp1d(mask_1d);
               rm_Q = cc_1d.Size > max_size;
               if any(rm_Q)
                   rm_ind = cat(1, cc_1d.PixelIdxList{rm_Q});
                   mask_1d(rm_ind) = false;
               end                
            end
                        
        end        
        
    end