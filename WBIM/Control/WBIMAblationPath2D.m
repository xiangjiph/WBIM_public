classdef WBIMAblationPath2D < matlab.mixin.Copyable & dynamicprops
    % Decompose arbitrary 2D regions into multiple parallel stripes. Each
    % stripe correspond to a laser ablation path
    % X, Y here are in stage coordinate
    properties(Constant)
        % axis indices in xyz order: 1 = x, 2 = y, 3 = z
        fast_axis_id (1,1) double = WBIMConfig.STAGE_AXIS_FAST_ID
        slow_axis_id (1,1) double = WBIMConfig.STAGE_AXIS_SLOW_ID
        z_axis_id (1,1) double = WBIMConfig.STAGE_AXIS_Z_ID   
        stream_axis_id = 1 : 3;
    end
    
    properties
        fast_v_um_s (1,1) double
        fast_step_um (1,1) double
        fast_acc_len_um (1,1) double
        fast_acc_m_s2 (1,1) double
        slow_step_um (1,1) double        
        abl_fp_abs_z_um (1, :) double
        abl_peak_fluence_J_cm2 (1, :) double
        
        % Ablation region
        bbox_xy_mmxx_um (1, 4) double
        bbox_xy_ll_um (1, 2) double
        mask_um_xy logical
        
        % Path information
        slow_ep_um double
        fast_seg_ep_um % could be double or cell? 
        fast_seg_ep_idx
        
        pretrigger_offset_distance_um (1,1) double
        num_line (1,1) double
    end
    
    properties(Constant, Hidden)
        fast_axis_id_0 (1,1) double = WBIMConfig.STAGE_AXIS_FAST_ID - 1
        slow_axis_id_0 (1,1) double = WBIMConfig.STAGE_AXIS_SLOW_ID - 1
        z_axis_id_0 (1,1) double = WBIMConfig.STAGE_AXIS_Z_ID - 1
    end
    %%
    methods
        function obj = WBIMAblationPath2D(fast_v_um_s, fast_acc_len_um, slow_step_um, ...
                abl_fp_abs_z_um, abl_peak_fluence_J_cm2, stage_bbox_xy_um, abl_mask_um_xy)
            % WBIMAblationPath2D converts ablation parameters to low-level
            % motion instructions for the stage controller. 
            % Input: 
            %   fast_v_um_s: speed of the translation stage in the fast
            %   axis (um/s)
            %   fast_acc_len_um: acceleration length for the fast stage to
            %   reach the specificed maximum speed (um)
            %   slow_step_um: spacing between the stripes, also the step
            %   size for the slow axis (um)
            %   abl_fp_abs_z_um: absolute z position of the ablation plane,
            %   defined as the z_stage + z_piezo (um)
            %   abl_bbox_xy_mmxx_um: bounding box of the ablation region in
            %   stage xy coordinate (um)           
            %   abl_mask_um_xy: 2D logical array of the ablation region in
            %   the boundign box, in stage xy coordinate (um)
            if nargin < 7
                abl_mask_um_xy = [];
            end
            validateattributes(fast_v_um_s, {'numeric'}, {'positive', 'scalar'});
            validateattributes(fast_acc_len_um, {'numeric'}, {'positive', 'scalar'});
            validateattributes(slow_step_um, {'numeric'}, {'positive', 'scalar'});
            validateattributes(stage_bbox_xy_um, {'numeric'}, {'nonnegative',...
                'numel', 4});
            validateattributes(abl_fp_abs_z_um, {'numeric'}, {'nonnegative'});
            validateattributes(abl_peak_fluence_J_cm2, {'numeric'}, {'nonnegative'});
            
            assert(all(stage_bbox_xy_um(1:2) <= stage_bbox_xy_um(3:4)), 'Empty stage bounding box');
            if ~isempty(abl_mask_um_xy)
                validateattributes(abl_mask_um_xy, {'logical'}, {'2d'});
                abl_bbox_ll_um = stage_bbox_xy_um(3:4) - ...
                    stage_bbox_xy_um(1:2) + 1;
                assert(all(size(abl_mask_um_xy) == abl_bbox_ll_um), ['...' ...
                    'Ablation bounding box size should be the same as the ablation mask size']);
            end
            if ~isscalar(abl_fp_abs_z_um)
                if isscalar(abl_peak_fluence_J_cm2)
                    abl_peak_fluence_J_cm2 = repelem(abl_peak_fluence_J_cm2, ...
                        1, numel(abl_fp_abs_z_um));
                else
                    assert(numel(abl_peak_fluence_J_cm2) == numel(abl_fp_abs_z_um), ...
                        'Number of ablation z positions should be the same as the number of ablation fluences');
                end
            end
            
            obj.fast_v_um_s = fast_v_um_s;
            obj.fast_step_um = fast_v_um_s / WBIMConfig.ABLATION_REPETITION_RATE_Hz;
            obj.fast_acc_len_um = fast_acc_len_um;
            obj.fast_acc_m_s2 = (obj.fast_v_um_s/1e6)^2 / 2 / (obj.fast_acc_len_um/1e6);            
            obj.slow_step_um = round(slow_step_um);
            obj.abl_fp_abs_z_um = abl_fp_abs_z_um;
            obj.abl_peak_fluence_J_cm2 = abl_peak_fluence_J_cm2;
            obj.bbox_xy_mmxx_um = stage_bbox_xy_um;
            obj.bbox_xy_ll_um = obj.bbox_xy_mmxx_um(3:4) - obj.bbox_xy_mmxx_um(1:2) + 1;
            obj.mask_um_xy = abl_mask_um_xy;
            obj.pretrigger_offset_distance_um = fast_v_um_s * WBIMConfig.ABLATION_ZABER_PRETRIGGER_TIME_s;
            
            if ~isempty(obj.mask_um_xy)
                obj.construct_from_mask();
            else
                obj.construct_from_bbox();
            end
        end
    end    
    %% Construct
    methods
        function obj = construct_from_bbox(obj, bbox_xy_mmxx_um)
            if nargin < 2
                bbox_xy_mmxx_um = obj.bbox_xy_mmxx_um;
            else
                obj.bbox_xy_mmxx_um = bbox_xy_mmxx_um;
            end
            bbox_size_xy_um = bbox_xy_mmxx_um(3:4) - bbox_xy_mmxx_um(1:2) + 1;
            bbox_mask_um_xy = true(bbox_size_xy_um);
            obj.construct_from_mask(bbox_mask_um_xy);
        end
        
        function obj = construct_from_mask(obj, mask_um_xy)
            % This XY is in the stage coordinate, not sample coordinate
           if nargin < 2
               mask_um_xy = obj.mask_um_xy;
           else
               obj.mask_um_xy = mask_um_xy;
           end
           max_movmax_wd_um = 10;
           mask_um_fs = permute(mask_um_xy, [obj.fast_axis_id, obj.slow_axis_id]);
           roi_size_um = size(mask_um_fs);
           slow_valid_mask = any(mask_um_fs, 1);
           slow_start_ind = find(slow_valid_mask, 1, 'first');
           slow_end_ind = find(slow_valid_mask, 1, 'last');
                      
           slow_ep_1 = slow_start_ind : obj.slow_step_um : (slow_end_ind - obj.slow_step_um / 2);
           slow_ep_2 = slow_ep_1 + obj.slow_step_um - 1;
           slow_ctr = (slow_ep_1 + slow_ep_2) ./ 2;
           
           num_stripe = numel(slow_ep_1);
           num_valid_stripe = 0;
           fast_ep = cell(num_stripe, 1);
           int_ep_idx = cell(num_stripe, 1);
           is_valid_stripe_Q = false(num_stripe, 1);
           for i_stripe = 1 : num_stripe
               col_idx = (slow_ep_1(i_stripe) + 1) : ...
                   min((slow_ep_2(i_stripe) + 1), roi_size_um(2));
               mask_stripe = mask_um_fs(:, col_idx);
               mask_stripe = any(mask_stripe, 2);
               mask_stripe = movmax(mask_stripe, min(max_movmax_wd_um, round(obj.fast_step_um) * 2));
               % Be conservative...
%                mask_stripe = movmean(mask_stripe, round(obj.slow_step_um)) > 0.75;
               if any(mask_stripe)
                   is_valid_stripe_Q(i_stripe) = true;
                   num_valid_stripe = num_valid_stripe + 1;
                   % Determine the position to start and end ablation
                   mask_stripe_p = cat(1, 0, mask_stripe, 0);
                   inflect_pos = find(diff(mask_stripe_p)).';
                   assert(issorted(inflect_pos, 'ascend'), 'The indices are not in ascending order');
                   num_inflect_pts = numel(inflect_pos);
                   assert(mod(num_inflect_pts, 2) == 0, 'Missing endpoint');          
                   % Position along the fast axis, in um, in ascending
                   % order
                   if isrow(inflect_pos)
                       inflect_pos = inflect_pos.';
                   end
                   tmp_fast_pos_um_r = cat(1, inflect_pos(1) - obj.fast_acc_len_um, ...
                       inflect_pos, inflect_pos(end) + obj.fast_acc_len_um);
                   if mod(num_valid_stripe, 2) == 0
                      tmp_fast_pos_um_r = flip(tmp_fast_pos_um_r);                       
                   end
                 
                   fast_ep{num_valid_stripe} = tmp_fast_pos_um_r + ...
                       obj.bbox_xy_mmxx_um(obj.fast_axis_id);
                   num_int = num_inflect_pts / 2;
                   
                   ep_idx = cat(1, 0, repmat([1;-1], num_int, 1), 0);
                   int_ep_idx{num_valid_stripe} = ep_idx;
               end
           end
           obj.num_line = num_valid_stripe;
           obj.slow_ep_um = slow_ctr(is_valid_stripe_Q) + obj.bbox_xy_mmxx_um(obj.slow_axis_id);
           obj.fast_seg_ep_um = fast_ep(1 : num_valid_stripe);
           obj.fast_seg_ep_idx = int_ep_idx(1 : num_valid_stripe);
        end
        
        function xyz_um_tQ = get_single_plane_trajectory(obj, z_um)
            if nargin < 2
                z_um = nan;
            end
            xyq = obj.get_single_plane_trajectory_xyq();
            num_pts = size(xyq, 1);
            xyz_um_tQ = cat(2, xyq(:, [1,2]), repelem(z_um, num_pts, 1), ...
                xyq(:, end));
        end
        
        function xy_um_tQ = get_single_plane_trajectory_xyq(obj)
            if iscell(obj.fast_seg_ep_um)
                num_pts_per_slow = cellfun(@numel, obj.fast_seg_ep_um);
            elseif ismatrix(obj.fast_seg_ep_um)
                num_pts_per_slow = size(obj.fast_seg_ep_um, 1);
            end
            slow_pos_um = repelem(obj.slow_ep_um(:), num_pts_per_slow, 1);
            fast_pos_um = cat(1, obj.fast_seg_ep_um{:});
            xy_um_tQ = zeros(numel(slow_pos_um), 3);
            xy_um_tQ(:, obj.fast_axis_id) = fast_pos_um;
            xy_um_tQ(:, obj.slow_axis_id) = slow_pos_um;
            xy_um_tQ(:, end) = cat(1, obj.fast_seg_ep_idx{:});
        end
        
        function xyz_um_tQ = get_multi_planes_trajectory(obj, z_um_list)
            xyz_um_tQ0 = obj.get_single_plane_trajectory(z_um_list(1));
            % Reverse order
            num_z = numel(z_um_list);
            num_pts = size(xyz_um_tQ0, 1);
            xyz_um_tQ = zeros(num_z * num_pts, 4);
            xyz_um_tQ(1 : num_pts, :) = xyz_um_tQ0;
            for i = 2 : num_z
                tmp_idx = (1 : num_pts) + num_pts * (i - 1);
                if mod(i, 2) == 0
                    xyz_um_tQ(tmp_idx, 1:2) = flip(xyz_um_tQ0(:, 1:2), 1);
                else
                    xyz_um_tQ(tmp_idx, 1:2) = xyz_um_tQ0(:, 1:2);
                end
                xyz_um_tQ(tmp_idx, 3) = z_um_list(i);
                % The order has been flipped above. No need to multiple by
                % -1 here
                xyz_um_tQ(tmp_idx, 4) = xyz_um_tQ0(:, 4);
            end
        end
        
        %%
        function xy_um_tQ_s = add_pretrigger_compensation(obj, xyz_um_tQ)
            % Fixed bug. Need to copy the array. 
            xy_um_tQ_s = xy_um_tQ;
            num_pts = size(xyz_um_tQ, 1);
            assert(xyz_um_tQ(1, end) == 0, 'The first point should not change the controller trigger output');
            for i = 2 : num_pts
                if xyz_um_tQ(i, end) ~= 0
                    assert(xyz_um_tQ(i-1, obj.slow_axis_id) == xyz_um_tQ(i, obj.slow_axis_id) && ...
                        xyz_um_tQ(i+1, obj.slow_axis_id) == xyz_um_tQ(i, obj.slow_axis_id), 'These three points should be colinear.')
                    % Change DIO
                    if xyz_um_tQ(i, obj.fast_axis_id) > xyz_um_tQ(i-1, obj.fast_axis_id)
                        % If moving in the positive direction 
                        xy_um_tQ_s(i, obj.fast_axis_id) = xyz_um_tQ(i, obj.fast_axis_id) + ...
                            obj.pretrigger_offset_distance_um;
                    elseif xyz_um_tQ(i, obj.fast_axis_id) < xyz_um_tQ(i-1, obj.fast_axis_id)
                        % Move in the negative direction 
                        xy_um_tQ_s(i, obj.fast_axis_id) = xyz_um_tQ(i, obj.fast_axis_id) - ...
                            obj.pretrigger_offset_distance_um;
                    end
                end
            end
        end 
    end
    %%
    methods(Static)
        function xyz_um_tQ_r = reverse_trajectory(xyz_um_tQ)
           xyz_um_tQ_r = xyz_um_tQ(end:-1:1, :);
           xyz_um_tQ_r(:,end) = - xyz_um_tQ_r(:,end);
        end
               
        %% Visualization
        function [fig_hdl, ax_hdl] = visualize_trajectory(xyz_um_tQ)
            num_dim = size(xyz_um_tQ, 2) - 1;
            
            fig_hdl = figure;
            ax_hdl = axes(fig_hdl);
            % imagesc(ax_hdl, roi_data.im);
            switch num_dim
                case 3
                    plot3(ax_hdl, xyz_um_tQ(:, 1), xyz_um_tQ(:, 2), xyz_um_tQ(:, 3),...
                        'LineWidth', 1);
                    hold(ax_hdl, 'on');
                    s_hdl_1 = scatter3(ax_hdl, xyz_um_tQ(xyz_um_tQ(:, end) == 1, 1), ...
                        xyz_um_tQ(xyz_um_tQ(:, end) == 1, 2), ...
                        xyz_um_tQ(xyz_um_tQ(:, end) == 1, 3), 100, 'r.');
                    s_hdl_2 = scatter3(ax_hdl, xyz_um_tQ(xyz_um_tQ(:, end) == -1, 1), ...
                        xyz_um_tQ(xyz_um_tQ(:, end) == -1, 2), ...
                        xyz_um_tQ(xyz_um_tQ(:, end) == -1, 3), 100, 'b.');
                    s_hdl_3 = scatter(ax_hdl, xyz_um_tQ(1,1), xyz_um_tQ(1,2), ...
                        xyz_um_tQ(1,3), 100, 'g.');
                    ax_hdl.ZLabel.String = 'Z (\mum)';
                case 2
                    plot(ax_hdl, xyz_um_tQ(:, 1), xyz_um_tQ(:, 2), ...
                        'LineWidth', 1);
                    hold(ax_hdl, 'on');
                    s_hdl_1 = scatter(ax_hdl, xyz_um_tQ(xyz_um_tQ(:, end) == 1, 1), ...
                        xyz_um_tQ(xyz_um_tQ(:, end) == 1, 2), 100, 'r.');
                    s_hdl_2 = scatter(ax_hdl, xyz_um_tQ(xyz_um_tQ(:, end) == -1, 1), ...
                        xyz_um_tQ(xyz_um_tQ(:, end) == -1, 2), 100, 'b.');
                    s_hdl_3 = scatter(ax_hdl, xyz_um_tQ(1,1), xyz_um_tQ(1,2), 100, 'g.');
            end
            ax_hdl.XLabel.String = 'X (\mum)';
            ax_hdl.YLabel.String = 'Y (\mum)';            
            legend(ax_hdl, [s_hdl_1, s_hdl_2, s_hdl_3], 'Trigger on', 'Trigger off',...
                'Start', 'Location', 'best');
            
%             fig_fp = fullfile(DataManager.DATA_ROOT_PATH, 'Visualization', 'Pipeline', ...
%                 'Ablation_path_exampel.png');
%             fun_print_image_in_several_formats(fig_hdl, fig_fp);
        end
    end
    %% Visualization 
    methods
        function fig_hdl = visualize_mask_with_trajectory(obj, im_stage_yx_um, fullViewQ)
            if nargin < 3
                fullViewQ = false;
            end
            abl_trajectory = obj.get_single_plane_trajectory_xyq();
            if isempty(abl_trajectory)
                fig_hdl = [];
                fprintf('The trajectory is empty\n');
                return;
            else
                if nargin < 2 || isempty(im_stage_yx_um)
                    im_stage_yx_um = obj.mask_um_xy.';
                end
                abl_trajectory_r = abl_trajectory - [obj.bbox_xy_mmxx_um(1:2), 0];
                
                fig_hdl = figure('Units', 'normalized', 'Position', [0.1, 0.1, 0.8, 0.8], 'Visible', 'off');
                ax_hdl = axes(fig_hdl);
                imagesc(ax_hdl, im_stage_yx_um);
                colormap(ax_hdl, 'gray');
                ax_hdl.Color = [0.5, 0.5, 0.5];
                ax_hdl.DataAspectRatio = [1,1,1];
                ax_hdl.YLabel.String = 'Stage 2 (\mum)';
                ax_hdl.XLabel.String = 'Stage 1 (\mum)';
                hold(ax_hdl, 'on');
                plt_hdl = plot(ax_hdl, abl_trajectory_r(:, 1), abl_trajectory_r(:, 2), ...
                    'LineWidth', 1, 'Color', '#D95319');
                sc_hdl_1 = scatter(ax_hdl, abl_trajectory_r(abl_trajectory(:, 3) == 1, 1),...
                    abl_trajectory_r(abl_trajectory(:, 3) == 1, 2), 200, 'r.');
                sc_hdl_2 = scatter(ax_hdl, abl_trajectory_r(abl_trajectory(:, 3) == -1, 1),...
                    abl_trajectory_r(abl_trajectory(:, 3) == -1, 2), 200, 'b.');
                sc_hdl_3 = line(ax_hdl, abl_trajectory_r(1, 1), abl_trajectory_r(1, 2), ...
                    'Color', 'g', 'Marker', '*', 'LineWidth', 1, 'LineStyle', 'none');
                if fullViewQ
                    ax_hdl.YLim = [min(abl_trajectory_r(:, 2)), max(abl_trajectory_r(:, 2))]...
                        + [-50, 50];
                    ax_hdl.XLim = ax_hdl.XLim + [-50, 50];
                end
                legend(ax_hdl, [plt_hdl, sc_hdl_3, sc_hdl_1, sc_hdl_2], ...
                    'Ablation path', 'Starting point', 'Ablation on', 'Ablation off');
                ax_hdl.XTickLabel = arrayfun(@(x) num2str(round(x), '%d'), ...
                    ax_hdl.XTick + obj.bbox_xy_mmxx_um(1), 'UniformOutput', false);
                ax_hdl.YTickLabel = arrayfun(@(x) num2str(round(x), '%d'), ...
                    ax_hdl.YTick + obj.bbox_xy_mmxx_um(1), 'UniformOutput', false);
                
                ax_hdl.Title.String = 'Planned ablation path in stage coordinate';
                fig_hdl.Visible = 'on';
            end
        end
    end
end