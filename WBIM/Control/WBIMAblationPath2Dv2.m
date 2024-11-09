classdef WBIMAblationPath2Dv2 < matlab.mixin.Copyable & dynamicprops
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
        fast_uncertainty_length_um (1,1) double
        slow_step_um (1,1) double        
        abl_fp_abs_z_um (1, 1) double
        power_W (1, 1) double
        diffuserQ (1,1) logical
        
        % Ablation region
        bbox_xy_mmxx_um (1, 4) double
        bbox_xy_ll_um (1, 2) double
        mask_um_xy logical
        
        pretrigger_offset_distance_um (1,1) double
        
        cpp (1, :) CoveragePathPlanning
        start_xy_um (1, :) double = []
    end
    
    properties
        trajectory (:, 3) {mustBeNumeric}
        is_not_empty_Q (1,1) logical
    end
    
    properties(Hidden)
        trajectory_i (:, 3) {mustBeNumeric} % idealized trajectory without stage delay correction
    end
    
    properties(Constant, Hidden)
        fast_axis_id_0 (1,1) double = WBIMConfig.STAGE_AXIS_FAST_ID - 1
        slow_axis_id_0 (1,1) double = WBIMConfig.STAGE_AXIS_SLOW_ID - 1
        z_axis_id_0 (1,1) double = WBIMConfig.STAGE_AXIS_Z_ID - 1
    end
    %%
    methods
        function obj = WBIMAblationPath2Dv2(fast_v_um_s, fast_acc_len_um, slow_step_um, ...
                abl_fp_abs_z_um, power_W, diffuserQ, stage_bbox_xy_um, abl_mask_um_xy)
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
            arguments
                fast_v_um_s (1,1) {mustBePositive, mustBeNumeric}
                fast_acc_len_um (1,1) {mustBePositive, mustBeNumeric}
                slow_step_um (1,1) {mustBePositive, mustBeNumeric}
                abl_fp_abs_z_um (1,1) {mustBeNonnegative, mustBeNumeric}
                power_W (1,1) {mustBeNumeric, mustBeNonnegative}
                diffuserQ (1,1) {logical}
                stage_bbox_xy_um (1, 4) {mustBeNumeric}
                abl_mask_um_xy = []
            end
            assert(all(stage_bbox_xy_um(1:2) <= stage_bbox_xy_um(3:4)), 'Empty stage bounding box');
            if ~isempty(abl_mask_um_xy)
                validateattributes(abl_mask_um_xy, {'logical'}, {'2d'});
                abl_bbox_ll_um = stage_bbox_xy_um(3:4) - ...
                    stage_bbox_xy_um(1:2) + 1;
                assert(all(size(abl_mask_um_xy) == abl_bbox_ll_um), ['...' ...
                    'Ablation bounding box size should be the same as the ablation mask size']);
            end
            
            obj.fast_v_um_s = fast_v_um_s;
            obj.fast_step_um = fast_v_um_s / WBIMConfig.ABLATION_REPETITION_RATE_Hz;
            obj.fast_acc_len_um = fast_acc_len_um;
            obj.fast_acc_m_s2 = (obj.fast_v_um_s/1e6)^2 / 2 / (obj.fast_acc_len_um/1e6);            
            obj.slow_step_um = round(slow_step_um);
            obj.abl_fp_abs_z_um = abl_fp_abs_z_um;
            obj.power_W = power_W;
            obj.diffuserQ = diffuserQ;
            obj.bbox_xy_mmxx_um = stage_bbox_xy_um;
            obj.bbox_xy_ll_um = obj.bbox_xy_mmxx_um(3:4) - obj.bbox_xy_mmxx_um(1:2) + 1;
            obj.mask_um_xy = abl_mask_um_xy;
            obj.pretrigger_offset_distance_um = fast_v_um_s * WBIMConfig.ABLATION_ZABER_PRETRIGGER_TIME_s;
            obj.fast_uncertainty_length_um = round(obj.fast_v_um_s *...
                WBIMConfig.ABLATION_ZABER_PRETRIGGER_TIME_UNCERTAINTY_ms/1e3);
            
            if isempty(obj.mask_um_xy)
                obj.mask_um_xy = true(obj.bbox_xy_ll_um);
            end
        end
    end 
    %% Get
    methods
        function val = get.is_not_empty_Q(obj)
            val = ~isempty(obj.trajectory);
        end
    end
    
    
    %% Construct
    methods                       
        %%      
        function obj = coverage_path_planning(obj, options)
            arguments
                obj (1,1) WBIMAblationPath2Dv2
                options.stripe_offset (1,1) double = 0
                options.max_merge_dist_um (1,1) double = max(1e3, 4 * obj.fast_acc_len_um)
            end
            mask_um_fs = permute(obj.mask_um_xy, [obj.fast_axis_id, obj.slow_axis_id]);
            % Dilate the mask along the fast direction based on the
            % uncertainty of stage acceleration time 
            mask_um_fs = imdilate(mask_um_fs, strel('rectangle', [obj.fast_uncertainty_length_um, 1]));
            
            obj.cpp = CoveragePathPlanning(obj.fast_acc_len_um, obj.slow_step_um, mask_um_fs);
            obj.cpp.mask_to_stripe(options.stripe_offset);
            % Merge lines in each stripe
            obj.cpp.merge_lines_by_distance(options.max_merge_dist_um);
            obj.cpp.merge_interval_by_fracitonal_overlap();
            obj.cpp.compute_interval_point_sfq_for_each_cell();
        end        

        function t_xyq = construct_optimal_trajectory(obj, options)
            arguments
                obj (1,1) WBIMAblationPath2Dv2
                options.start_xy_um (1, :) double = obj.start_xy_um; % in global stage coordinate
            end
            if ~isempty(options.start_xy_um)
                assert(numel(options.start_xy_um) == 2, 'The starting posiiton should be 2D')
                start_local_xy_um = options.start_xy_um - ...
                    obj.bbox_xy_mmxx_um(1:2) + 1;
                start_sf_um = start_local_xy_um([obj.slow_axis_id, obj.fast_axis_id]);
            else
                start_sf_um = [];
            end
            trajectory_sfq = obj.cpp.compute_optimal_trajectory('start_sf_um', start_sf_um);
            % Transform to global stage coordinate
            obj.trajectory = zeros(size(trajectory_sfq));
            obj.trajectory(:, obj.fast_axis_id) = trajectory_sfq(:, 2) + ...
                obj.bbox_xy_mmxx_um(obj.fast_axis_id) - 1;
            obj.trajectory(:, obj.slow_axis_id) = trajectory_sfq(:, 1) + ...
                obj.bbox_xy_mmxx_um(obj.slow_axis_id) - 1;
            obj.trajectory(:, end) = trajectory_sfq(:, end);
            obj.start_xy_um = options.start_xy_um;
            
            t_xyq = obj.trajectory;
        end
        
        %%
        function xy_um_tQ_s = add_pretrigger_compensation(obj, xy_um_tQ)
            % Zaber controller triggers digital output change based on the
            % idealized trajectory. 
            xy_um_tQ_s = xy_um_tQ;
            num_pts = size(xy_um_tQ, 1);
            assert(xy_um_tQ(1, end) == 0, 'The first point should not change the controller trigger output');
            for i = 2 : num_pts
                if xy_um_tQ(i, end) ~= 0 % || (xyz_um_tQ(i-1, end) == -1) - shift the last point in the stripe? 
                    assert(xy_um_tQ(i-1, obj.slow_axis_id) == xy_um_tQ(i, obj.slow_axis_id) && ...
                        xy_um_tQ(i+1, obj.slow_axis_id) == xy_um_tQ(i, obj.slow_axis_id), 'These three points should be colinear.')
                    % Change DIO
                    if xy_um_tQ(i, obj.fast_axis_id) > xy_um_tQ(i-1, obj.fast_axis_id)
                        % If moving in the positive direction
                        xy_um_tQ_s(i, obj.fast_axis_id) = xy_um_tQ(i, obj.fast_axis_id) + ...
                            obj.pretrigger_offset_distance_um;
                    elseif xy_um_tQ(i, obj.fast_axis_id) < xy_um_tQ(i-1, obj.fast_axis_id)
                        % Move in the negative direction
                        xy_um_tQ_s(i, obj.fast_axis_id) = xy_um_tQ(i, obj.fast_axis_id) - ...
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
            abl_trajectory = obj.trajectory_i;
            if isempty(abl_trajectory)
                fig_hdl = [];
                fprintf('The trajectory is empty\n');
                return;
            else
                if nargin < 2 || isempty(im_stage_yx_um)
                    im_stage_yx_um = obj.mask_um_xy.';
                end
                
                abl_trajectory_r = abl_trajectory - [obj.bbox_xy_mmxx_um(1:2), 0];
                
%                 fig_hdl = figure('Units', 'normalized', 'Position', [0.1, 0.1, 0.8, 0.8]);
                fig_hdl = figure; 
                fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* 2;
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
                sc_hdl_4 = line(ax_hdl, abl_trajectory_r(end, 1), abl_trajectory_r(end, 2), ...
                    'Color', 'y', 'Marker', '*', 'LineWidth', 1, 'LineStyle', 'none');
                if fullViewQ
                    ax_hdl.YLim = [min(abl_trajectory_r(:, 2)), max(abl_trajectory_r(:, 2))]...
                        + [-50, 50];
                    ax_hdl.XLim = ax_hdl.XLim + [-50, 50];
                end
                legend(ax_hdl, [plt_hdl, sc_hdl_3, sc_hdl_4, sc_hdl_1, sc_hdl_2], ...
                    'Ablation path', 'Starting point', 'End point', 'Ablation on', 'Ablation off');
                ax_hdl.XTickLabel = arrayfun(@(x) num2str(round(x), '%d'), ...
                    ax_hdl.XTick + obj.bbox_xy_mmxx_um(1), 'UniformOutput', false);
                ax_hdl.YTickLabel = arrayfun(@(x) num2str(round(x), '%d'), ...
                    ax_hdl.YTick + obj.bbox_xy_mmxx_um(2), 'UniformOutput', false);
                
                ax_hdl.Title.String = 'Planned ablation path in stage coordinate';
            end
        end
    end
end