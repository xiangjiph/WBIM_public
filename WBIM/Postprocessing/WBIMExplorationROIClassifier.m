classdef WBIMExplorationROIClassifier < handle
    
    
    methods(Static)
        function mip = parse_buffered_image_stack(buffered_data)
            num_ch = numel(buffered_data);
            mip = cell(1, num_ch);
            for i = 1 : num_ch
                tmp_data = buffered_data{i};
                num_sec = size(tmp_data, 3);
                for j = 1 : num_sec
                    % Reduce shot noise
                    tmp_data(:, :, j) = medfilt2(tmp_data(:, :, j));
                end
                tmp_data = max(tmp_data, [], 3);
                mip{i} = tmp_data;
            end            
        end
        
        function stat_str = compute_single_tile_mip_stat(im)
            assert(ismatrix(im));
            im_max_val = single(intmax(class(im)));
            im = max(0, single(im));
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
            stat_str.bin_edge = 0 : 0.1 : 1;
            stat_str.nmpdf = histcounts(rescale(im_f), stat_str.bin_edge, 'Normalization', 'probability');
            stat_str.mpdf = histcounts(im_f ./ im_max_val, stat_str.bin_edge, 'Normalization', 'probability');
            stat_str.mcdf = cumsum(stat_str.mpdf);
        end
        
        function is_bg_Q = classify_tiles_by_rules(mip, options)
            arguments
                mip (1, :) cell
                options.bg_10_ptrl_frac (1,1) double = 0.975
                options.bg_max_mean_int (1,1) double = 4e3;
            end
            num_ch = numel(mip);
%             stat_info = cell(num_ch, 1);
            is_bg_Q = true;
            for i = 1 : num_ch
                % Compute mip statistics
                tmp_stat = WBIMExplorationROIClassifier.compute_single_tile_mip_stat(...
                    mip{i});
                is_bg_Q = is_bg_Q && (tmp_stat.mcdf(1) > options.bg_10_ptrl_frac)...
                    && (tmp_stat.mean_raw < options.bg_max_mean_int);
            end
        end        
    end
    
end