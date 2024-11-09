classdef CoveragePathPlanning < handle
    
    % Input properties
    properties
        fast_acc_len_um (1,1) double
        slow_step_um (1,1) double
        mask_um_fs logical
    end
    % Derived properties
    properties
        % Path information
        num_stripe (1,1) double = 0
        slow_ep_um double
        
        fast_ep_um (:, 1) cell
        
        slow_ep_um_list (:, 1) double
        fast_ep_um_list (2, :) double
        
        num_line_per_stripe (:, 1) double
        stripe_line_idx_cell (:, 1) cell
        map_line_idx_2_stripe_idx (:, 1) double
        map_line_idx_2_int_idx (:, 1) double
        map_line_idx_2_cell_idx (:, 1) double
        num_line (1,1) double = 0
        % Merge nearby lines into a interval - based on spacing vs
        % acceleration time
        int_ep_line_idx (:, 1) cell
        int_ep_line_idx_list (2, :) double
        int_ep_line_local_idx (:, 1) cell
        int_ep_um (:, 1) cell
        num_int (1,1) double = 0
        num_int_per_stripe (:, 1) double
        stripe_int_idx_cell (:, 1) cell
        map_int_idx_2_stripe_idx (:, 1) double
        map_int_idx_2_cell_idx (:, 1) double
        % Merge nearby intervals into a cell - based on 1D connectivity
        num_cell (1,1) double = 0
        int_cell_label_list (:, 1) double
        num_int_per_cell
        cell_int_idx_cell (:, 1) cell
        cell_line_idx_cell (:, 1) cell
        cell_sfq(:, 1) cell
        % Path planning to tranverse cells
        
    end
    %%
    methods
        function obj = CoveragePathPlanning(fast_acc_len_um, ...
                slow_step_um, mask_um_fs)
            arguments
                fast_acc_len_um (1,1) double {mustBePositive}
                slow_step_um (1,1) double {mustBePositive}
                mask_um_fs {mustBeNumericOrLogical} % the fast axis is along the y direction (image)
            end
            obj.fast_acc_len_um = fast_acc_len_um;
            obj.slow_step_um = slow_step_um;
            obj.mask_um_fs = mask_um_fs;
            
        end
        %% Computation
        function obj = mask_to_stripe(obj, slow_offset_frac)
            arguments
                obj (1,1) CoveragePathPlanning
                slow_offset_frac (1,1) double = 0;
            end
            % Align each column with the fast axis.
            roi_size_um = size(obj.mask_um_fs);
            slow_valid_mask = any(obj.mask_um_fs, 1);
            slow_start_ind = find(slow_valid_mask, 1, 'first');
            slow_end_ind = find(slow_valid_mask, 1, 'last');
            
            slow_start_ind = slow_start_ind + slow_offset_frac * obj.slow_step_um;
            if (slow_end_ind - slow_start_ind) > obj.slow_step_um
                % More than 1 stripe 
                slow_ep_1 = round(slow_start_ind : obj.slow_step_um : ...
                    max(slow_start_ind, (slow_end_ind - obj.slow_step_um / 2)));
                slow_ep_2 = slow_ep_1 + round(obj.slow_step_um) - 1;
                slow_ctr = (slow_ep_1 + slow_ep_2) ./ 2;
            else
                slow_ctr = (slow_start_ind + slow_end_ind) / 2;
                slow_ep_1 = round(slow_ctr - obj.slow_step_um / 2);
                slow_ep_2 = round(slow_ctr + obj.slow_step_um / 2);
            end
            
            n_stripe = numel(slow_ep_1);
            num_valid_stripe = 0;
            obj.fast_ep_um = cell(n_stripe, 1);
            is_valid_stripe_Q = false(n_stripe, 1);
            obj.num_line_per_stripe = zeros(n_stripe, 1);
            for i_stripe = 1 : n_stripe
                col_idx = (max(1, slow_ep_1(i_stripe))) : ...
                    min((slow_ep_2(i_stripe)), roi_size_um(2));
                mask_stripe = obj.mask_um_fs(:, col_idx);
                mask_stripe = any(mask_stripe, 2);
                % Find intervals
                if any(mask_stripe)
                    is_valid_stripe_Q(i_stripe) = true;
                    num_valid_stripe = num_valid_stripe + 1;
                    % Determine the position to start and end ablation
                    mask_stripe_p = cat(1, 0, mask_stripe, 0);
                    inflect_pos = find(diff(mask_stripe_p)).';
                    assert(issorted(inflect_pos, 'ascend'), 'The indices are not in ascending order');
                    % Interval endpoint idx, each column is an interval
                    line_ep_sub_um = reshape(inflect_pos, 2, []);
                    line_ep_sub_um(2,:) = line_ep_sub_um(2,:) - 1;
                    assert(all(mask_stripe(line_ep_sub_um) == 1, 'all'), 'All mask value at the endpoints should be positive');
                    
                    obj.fast_ep_um{num_valid_stripe} = line_ep_sub_um;
                    obj.num_line_per_stripe(i_stripe) = size(line_ep_sub_um, 2);
                end
            end
            if num_valid_stripe > 0
                obj.num_stripe = num_valid_stripe;
                % Transform to the global coordinate
                %                 obj.slow_ep_um = slow_ctr(is_valid_stripe_Q) + obj.bbox_xy_mmxx_um(obj.slow_axis_id) - 1;
                obj.slow_ep_um = slow_ctr(is_valid_stripe_Q);
                
                obj.fast_ep_um = obj.fast_ep_um(1 : num_valid_stripe);
                obj.num_line_per_stripe = obj.num_line_per_stripe(is_valid_stripe_Q);
                obj.num_line = sum(obj.num_line_per_stripe);
                obj.stripe_line_idx_cell = mat2cell((1 : obj.num_line).', obj.num_line_per_stripe, 1);
                obj.map_line_idx_2_stripe_idx = repelem((1:num_valid_stripe).', ...
                    obj.num_line_per_stripe, 1);
                
                obj.slow_ep_um_list = repelem((obj.slow_ep_um).', ...
                    obj.num_line_per_stripe, 1);
                obj.fast_ep_um_list = cat(2, obj.fast_ep_um{:});
            end
        end
        
        function obj = merge_lines_by_distance(obj, max_merge_dist_um)
            arguments
                obj (1,1) CoveragePathPlanning
                max_merge_dist_um (1,1) double = 0;
            end
            if obj.num_line == 0
                return;
            end
            [obj.int_ep_line_local_idx, obj.int_ep_um, obj.int_ep_line_idx]...
                = deal(cell(obj.num_stripe, 1));
            obj.map_line_idx_2_int_idx = zeros(obj.num_line, 1);
            int_count = 0;
            for i_stripe = 1 : obj.num_stripe
                tmp_num_lines = obj.num_line_per_stripe(i_stripe);
                tmp_line_idx = obj.stripe_line_idx_cell{i_stripe};
                tmp_fast_ep_um = obj.fast_ep_um{i_stripe}.';
                % Determine if the adjacent segments should be merged
                tmp_int_dist_um = tmp_fast_ep_um(2:end, 1) - tmp_fast_ep_um(1:end-1, 2);
                tmp_split_Q = tmp_int_dist_um > max_merge_dist_um;
                
                tmp_num_int = nnz(tmp_split_Q) + 1;
                
                if tmp_num_int == 1
                    tmp_int_ep_local_idx = [1; tmp_num_lines];
                    tmp_int_ep_idx = obj.stripe_line_idx_cell{i_stripe}([1,end]);
                    if isrow(tmp_int_ep_idx)
                        tmp_int_ep_idx = tmp_int_ep_idx.';
                    end
                elseif tmp_num_int > 1
                    tmp_int_ep_local_idx = zeros(2, tmp_num_int);
                    tmp_split_ind = find(tmp_split_Q);
                    tmp_start_ind = 1;
                    for i_mint = 1 : (tmp_num_int - 1)
                        tmp_int_ep_local_idx(:, i_mint) = [tmp_start_ind; tmp_split_ind(i_mint)];
                        tmp_start_ind = tmp_split_ind(i_mint) + 1;
                    end
                    tmp_int_ep_local_idx(:, tmp_num_int) = [tmp_start_ind; tmp_num_lines];
                    tmp_int_ep_idx = obj.stripe_line_idx_cell{i_stripe}(tmp_int_ep_local_idx);
                end
                for j_int = 1 : tmp_num_int
                    obj.map_line_idx_2_int_idx(tmp_line_idx(...
                        tmp_int_ep_local_idx(1, j_int) : tmp_int_ep_local_idx(2, j_int))) = ...
                        j_int + int_count;
                end
                int_count = int_count + tmp_num_int;
                
                tmp_merged_int_ep_um = cat(1, tmp_fast_ep_um(tmp_int_ep_local_idx(1,:), 1).', ...
                    tmp_fast_ep_um(tmp_int_ep_local_idx(2,:), 2).'); % [[start, start...]; [end, end...];]
                
                obj.num_int_per_stripe(i_stripe) = tmp_num_int;
                obj.int_ep_line_local_idx{i_stripe} = tmp_int_ep_local_idx;
                obj.int_ep_line_idx{i_stripe} = tmp_int_ep_idx;
                obj.int_ep_um{i_stripe} = tmp_merged_int_ep_um;
            end
            obj.num_int = sum(obj.num_int_per_stripe);
            obj.stripe_int_idx_cell = mat2cell((1:obj.num_int).', obj.num_int_per_stripe, 1);
            obj.map_int_idx_2_stripe_idx = repelem((1:obj.num_stripe).', obj.num_int_per_stripe, 1);
            obj.int_ep_line_idx_list = cat(2, obj.int_ep_line_idx{:});
        end
        
        function obj = merge_interval_by_fracitonal_overlap(obj, min_merge_frac)
            arguments
                obj (1,1) CoveragePathPlanning
                min_merge_frac (1,1) double = 0.5
            end
            int_ep_um_list = cat(2, obj.int_ep_um{:}); % 2-by-N array
            obj.num_cell = 0;
            obj.int_cell_label_list = zeros(obj.num_int, 1);
            for i_int = 1 : obj.num_int
                tmp_stripe = obj.map_int_idx_2_stripe_idx(i_int);
                tmp_mint_ep_um = int_ep_um_list(:, i_int);
                if ~obj.int_cell_label_list(i_int)
                    obj.num_cell = obj.num_cell + 1;
                    obj.int_cell_label_list(i_int) = obj.num_cell;
                end
                % next neighboring stipes:
                if tmp_stripe < obj.num_stripe
                    tmp_nn_stripe = tmp_stripe + 1;
                    % Get the stripe endpoints
                    tmp_ns_mint_label = obj.stripe_int_idx_cell{tmp_nn_stripe};
                    for k_nint = 1 : numel(tmp_ns_mint_label)
                        tmp_nmint = tmp_ns_mint_label(k_nint);
                        tmp_nmint_ep_um = int_ep_um_list(:, tmp_nmint);
                        % Compute the overlapping ratio
                        tmp_overlap = min(tmp_mint_ep_um(2), tmp_nmint_ep_um(2)) - ...
                            max(tmp_mint_ep_um(1), tmp_nmint_ep_um(1));
                        tmp_overlap_1 = tmp_overlap / (tmp_mint_ep_um(2) - tmp_mint_ep_um(1));
                        tmp_overlap_2 = tmp_overlap / (tmp_nmint_ep_um(2) - tmp_nmint_ep_um(1));
                        if tmp_overlap_1 > min_merge_frac && tmp_overlap_2 > min_merge_frac
                            obj.int_cell_label_list(tmp_nmint) = obj.int_cell_label_list(i_int);
                        end
                    end
                end
            end
            obj.cell_int_idx_cell = fun_bin_data_to_idx_list(obj.int_cell_label_list);
            obj.num_int_per_cell = cellfun(@numel, obj.cell_int_idx_cell);
            obj.map_int_idx_2_cell_idx = zeros(obj.num_int, 1);
            obj.map_int_idx_2_cell_idx(cat(2, obj.cell_int_idx_cell{:})) = ...
                repelem(1:obj.num_cell, obj.num_int_per_cell);
            % Line label in each cell
            obj.map_line_idx_2_cell_idx = obj.map_int_idx_2_cell_idx(obj.map_line_idx_2_int_idx);
            obj.cell_line_idx_cell = fun_bin_data_to_idx_list(obj.map_line_idx_2_cell_idx);
        end
        
        function obj = compute_interval_point_sfq_for_each_cell(obj)
            % 2D coordinate in local coordiante; without stage delay
            % correction
            obj.cell_sfq = cell(obj.num_cell, 1);
            for i_c = 1 : obj.num_cell
                tmp_int_idx_list = obj.cell_int_idx_cell{i_c};
                assert(issorted(tmp_int_idx_list, 'ascend'), ...
                    'The interval index should be in ascending order');
                tmp_num_int = numel(tmp_int_idx_list);
                tmp_int_sfq_c = cell(tmp_num_int, 1);
                for j_i = 1 : tmp_num_int
                    tmp_int_idx = tmp_int_idx_list(j_i);
                    tmp_int_ep_line_idx = obj.int_ep_line_idx_list(:, tmp_int_idx);
                    tmp_int_line_idx = tmp_int_ep_line_idx(1) : tmp_int_ep_line_idx(2);
                    tmp_num_lines = numel(tmp_int_line_idx);
                    tmp_stripe_idx = obj.map_int_idx_2_stripe_idx(tmp_int_idx);
                    tmp_slow_ep_um = obj.slow_ep_um(tmp_stripe_idx);
                    tmp_fast_ep_um = obj.fast_ep_um_list(:, tmp_int_line_idx);
                    assert(issorted(tmp_fast_ep_um, 'ascend'), ...
                        'Fast axis endpoints should be in ascending order');
                    tmp_fast_ep_um_p = cat(1, tmp_fast_ep_um(1) - obj.fast_acc_len_um, ...
                        tmp_fast_ep_um(:), tmp_fast_ep_um(end) + obj.fast_acc_len_um);
                    %                     if mod(j_i, 2) == 0
                    %                         tmp_fast_ep_um_p = tmp_fast_ep_um_p(end:-1:1);
                    %                     end
                    tmp_trigger_Q = cat(1, 0, repmat([1;-1], tmp_num_lines, 1), 0);
                    tmp_sfq = zeros((tmp_num_lines * 2 + 2), 3);
                    tmp_sfq(:, 1) = tmp_slow_ep_um;
                    tmp_sfq(:, 2) = tmp_fast_ep_um_p;
                    tmp_sfq(:, end) = tmp_trigger_Q;
                    
                    tmp_int_sfq_c{j_i} = tmp_sfq;
                end
                obj.cell_sfq{i_c} = tmp_int_sfq_c;
            end
        end
        
        function trajectory_sfq = compute_optimal_trajectory(obj, options)
            arguments
                obj (1,1) CoveragePathPlanning
                options.start_sf_um (1, :) double = []
                options.fix_start_Q (1,1) logical = false
                options.visQ (1,1) logical = false
                options.visGIFFp char = []
            end
            cell_corner = cellfun(@obj.get_cell_corner_points, ...
                obj.cell_sfq, 'UniformOutput', false);
            if ~isempty(options.start_sf_um)
                assert(numel(options.start_sf_um) == 2, 'The starting position should be 2D');
                cell_corner = cat(1, {repmat(options.start_sf_um, 4, 1)}, cell_corner);
                options.fix_start_Q = true; % overwrite 
            end
            cell_dist = obj.compute_cell_dist_from_corners(cell_corner);  
            cell_label = 1 : size(cell_dist, 1);
            % Sort by ditance to the first cell
            if ~isempty(options.start_sf_um)
                [~, d21_idx] = sort(cell_dist(:, 1), 'ascend');
                assert(d21_idx(1) == 1);
                cell_label = cell_label(d21_idx);
                cell_corner = cell_corner(d21_idx);
                cell_dist = obj.compute_cell_dist_from_corners(cell_corner);
            end
            sorted_cell_idx = TSPSolver.get_ordered_position(cell_dist, ...
                'FixStartQ', options.fix_start_Q, 'NoReturnQ', true);
            sorted_cell_idx = cell_label(sorted_cell_idx);
            if ~isempty(options.start_sf_um)
                assert(sorted_cell_idx(1) == 1);
                sorted_cell_idx = sorted_cell_idx(2:end) - 1;
            end
            trajectory_sfq = obj.construct_trajectory(obj.cell_sfq, sorted_cell_idx);
            
            if options.visQ
                if ~isempty(options.visGIFFp)
                    [folder, fn, ~] = fileparts(options.visGIFFp);
                    if ~isempty(folder) && ~isfolder(folder)
                        mkdir(folder)
                    end
                end
                fig_hdl = figure;
                fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* 2;
                ax_hdl = axes(fig_hdl);
                imagesc(ax_hdl, obj.mask_um_fs);
                colormap('gray');
                hold(ax_hdl, 'on');
                ax_hdl.DataAspectRatio = [1,1,1];
                ax_hdl.XLabel.String = 'Slow axis (\mum)';
                ax_hdl.YLabel.String = 'Fast axis (\mum)';
                for i_c = 1 : obj.num_cell
                    tmp_t = trajectory_sfq{i_c};
                    plt_hdl = plot(ax_hdl, tmp_t(:, 1), tmp_t(:, 2), 'LineWidth', 1);
                    tmp_c = tmp_t([1, end], 1:2);
                    scatter(ax_hdl, tmp_c(:, 1), tmp_c(:, 2), 'filled', 'MarkerFaceColor', plt_hdl.Color);
%                     drawnow();
                    frame =  frame2im(getframe(fig_hdl));
                    [imind, cm] = rgb2ind(frame, 256);
                    if i_c == 1
                        imwrite(imind, cm, options.visGIFFp, 'gif', 'LoopCount', inf);
                    else
                        imwrite(imind, cm, options.visGIFFp, 'gif', 'WriteMode', 'append');
                    end                    
                    pause(0.5);
                end
            end
            trajectory_sfq = cat(1, trajectory_sfq{:});
        end
    end
        
    methods(Static)
        function corner_xy = get_cell_corner_points(int_xyq_cell)
            p1 = int_xyq_cell{1}(1, 1:2);
            p2 = int_xyq_cell{1}(end, 1:2);
            p3 = int_xyq_cell{end}(1, 1:2);
            p4 = int_xyq_cell{end}(end, 1:2);
            corner_xy = cat(1, p1, p2, p3, p4);
        end
        
        function sf_um_tQ = construct_path_from_intervals(int_sfq_cell, options)
            arguments
                int_sfq_cell
                options.start_xy_um (1, 2) double = [nan, nan];
                options.start_corner_idx (1,1) double = nan;
            end
            if isfinite(options.start_corner_idx)
                start_idx = options.start_corner_idx;
            end
            if all(isfinite(options.start_xy_um))
                corner_xy = CoveragePathPlanning.get_cell_corner_points(int_sfq_cell);
                start_idx = find(all(abs(corner_xy - start_xy_um) < eps, 2));
            end
            
            assert(isscalar(start_idx));
            if start_idx > 2
                % Reverse the interval order
                int_sfq_cell = int_sfq_cell(end:-1:1);
            end
            if mod(start_idx, 2) == 0
                % start from 2 or 4 - reverse the order
                int_sfq_cell(1:2:end) = cellfun(@(x) WBIMAblationPath2Dv2.reverse_trajectory(x), int_sfq_cell(1:2:end), 'UniformOutput', false);
            else
                int_sfq_cell(2:2:end) = cellfun(@(x) WBIMAblationPath2Dv2.reverse_trajectory(x), int_sfq_cell(2:2:end), 'UniformOutput', false);
            end
            sf_um_tQ = cat(1, int_sfq_cell{:});
        end
        
        function cell_dist = compute_cell_dist_from_corners(cell_corner)
            % Input:
            %   cell_corner: cell array, each cell contains the 4 x 2
            %   matrix where each row is the position of a corner (left
            %   top; left button; right top; right button, in the MATLAB
            %   image coordinate)
            num_cell = numel(cell_corner);
            cell_corner_list = cat(1, cell_corner{:});
            cell_dist = zeros(num_cell);
            for i_c = 1 : num_cell
                tmp_pdist = pdist2(cell_corner{i_c}, cell_corner_list);
                % Find the shortest distance for each column -> tmp_min_dist_idx_1(i)
                % is the index of i_c's corner that is closest to a given point
                tmp_pdist_m_1 = min(tmp_pdist, [], 1);
                % Reshape, each column is a cell
                tmp_pdist_m_1 = reshape(tmp_pdist_m_1, 4, []);
                % For each cell, find the index of its corner that is closest to cell
                % i_c
                tmp_pdist_m_2 = min(tmp_pdist_m_1, [], 1);
                cell_dist(:, i_c) = tmp_pdist_m_2;
            end
        end
        
        function start_corner_idx = get_paired_corner_idx(...
                num_int, end_corner_idx)
            parity = mod(num_int, 2);
            switch end_corner_idx
                case 1
                    if parity
                        start_corner_idx = 4;
                    else
                        start_corner_idx = 3;
                    end
                case 2
                    if parity
                        start_corner_idx = 3;
                    else
                        start_corner_idx = 4;
                    end
                case 3
                    if parity
                        start_corner_idx = 2;
                    else
                        start_corner_idx = 1;
                    end
                case 4
                    if parity
                        start_corner_idx = 1;
                    else
                        start_corner_idx = 2;
                    end
            end
        end
        
        function trajectory = construct_trajectory(cell_xyq, sorted_cell_idx)
            cell_corner = cellfun(@CoveragePathPlanning.get_cell_corner_points, ...
                cell_xyq, 'UniformOutput', false);
            num_cell = numel(cell_xyq);
            trajectory = cell(num_cell, 1);
            for i_c = 1 : num_cell
                tmp_curr_idx = sorted_cell_idx(i_c);
                tmp_curr_int_cell = cell_xyq{tmp_curr_idx};
                tmp_curr_num_int = numel(tmp_curr_int_cell);
                if i_c == 1
                    if num_cell > 1
                        tmp_next_idx = sorted_cell_idx(i_c + 1);
                        % Cloest endpoint pair:
                        tmp_corner_dist = pdist2(cell_corner{tmp_curr_idx}, cell_corner{tmp_next_idx});
                        [~, min_ind] = min(tmp_corner_dist, [], 'all', 'linear');
                        [tmp_curr_end_ep_idx, ~] = ind2sub(size(tmp_corner_dist), min_ind);
                        tmp_curr_start_idx = CoveragePathPlanning.get_paired_corner_idx(...
                            tmp_curr_num_int, tmp_curr_end_ep_idx);
                    else
                        tmp_curr_start_idx = 1;
                    end
                else
                    tmp_prev_idx = sorted_cell_idx(i_c - 1);
                    tmp_prev_end_ep_xy = cell_corner{tmp_prev_idx}(tmp_curr_end_ep_idx, :);
                    tmp_corner_dist = pdist2(tmp_prev_end_ep_xy, cell_corner{tmp_curr_idx});
                    [~, tmp_curr_start_idx] = min(tmp_corner_dist);
                    tmp_curr_end_ep_idx = CoveragePathPlanning.get_paired_corner_idx(...
                        tmp_curr_num_int, tmp_curr_start_idx);
                end
                tmp_curr_xyt = CoveragePathPlanning.construct_path_from_intervals(...
                    tmp_curr_int_cell, 'start_corner_idx', tmp_curr_start_idx);
                % Get the starting points
                trajectory{i_c} = tmp_curr_xyt;
                % Check: to be implemented
                CoveragePathPlanning.check_trajectory(tmp_curr_xyt);
            end
        end
        
        function exit_code = check_trajectory(trajectory_c)
            num_points = size(trajectory_c, 1);
            assert(all(trajectory_c([1, end], end) == 0), 'The first and last trigger should be 0');
            direction = 1;
            for i = 2 : (num_points-1)
                switch trajectory_c(i-1, end)
                    case 1
                        assert(trajectory_c(i, end) == -1, 'The trigger state following 1 should be -1');
                        
                    case -1
                        assert(any(trajectory_c(i, end) == [1, 0]), 'The trigger state following -1 should be either 0 or 1');
                    case 0
                        assert(any(trajectory_c(i, end) == [1, 0]), 'The trigger state following 0 should be either 0 or 1');
                end
            end           
            exit_code = 0;
        end
    end
    
end