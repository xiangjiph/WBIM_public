classdef WBIMAblationPath3D < matlab.mixin.Copyable
    
    properties(Constant)
        % axis indices in xyz order: 1 = x, 2 = y, 3 = z
        fast_axis_id (1,1) double = WBIMConfig.STAGE_AXIS_FAST_ID
        slow_axis_id (1,1) double = WBIMConfig.STAGE_AXIS_SLOW_ID
        z_axis_id (1,1) double = WBIMConfig.STAGE_AXIS_Z_ID
        stream_axis_id = 1 : 3;
    end
    
    properties
        fast_v_um_s (1, :) double
        fast_acc_len_um (1, :) double
        slow_step_um (1, :) double        
        abl_fp_abs_z_um (1, :) double
        num_plane (1,1) double = 0
        power_W (1, :) double
        diffuserQ (1, :) logical
        
        % Ablation region
        bbox_xy_mmxx_um (1, 4) double
        bbox_xy_ll_um (1, 2) double
        mask_xy_um_stack logical
        
        path (:, 1) WBIMAblationPath2Dv2
    end    
    
    methods
        function obj = WBIMAblationPath3D(fast_v_um_s, fast_acc_len_um, slow_step_um, ...
                abl_fp_abs_z_um, power_W, diffuserQ, stage_bbox_xy_um, mask_xy_um_stack)
            arguments
                fast_v_um_s (1, :) {mustBeNumeric, mustBePositive}
                fast_acc_len_um (1, :) {mustBeNumeric, mustBePositive}
                slow_step_um (1, :) {mustBeNumeric, mustBePositive}
                abl_fp_abs_z_um (1, :) {mustBeNumeric, mustBePositive}
                power_W (1, :) {mustBeNumeric, mustBePositive}
                diffuserQ (1, :) {logical}
                stage_bbox_xy_um (1, 4) {mustBeNumeric, mustBePositive}
                mask_xy_um_stack {mustBeNumericOrLogical} = [] % 3D mask or empty
            end
            
            obj.abl_fp_abs_z_um = abl_fp_abs_z_um;
            obj.num_plane = numel(obj.abl_fp_abs_z_um);
            if isscalar(fast_v_um_s)
                fast_v_um_s = repelem(fast_v_um_s, 1, obj.num_plane);
            else
                assert(numel(fast_v_um_s) == obj.num_plane, 'Mismatch array size');
            end
            if isscalar(fast_acc_len_um)
                fast_acc_len_um = repelem(fast_acc_len_um, 1, obj.num_plane);
            else
                assert(numel(fast_acc_len_um) == obj.num_plane, 'Mismatch array size');
            end
            if isscalar(slow_step_um)
                slow_step_um = repelem(slow_step_um, 1, obj.num_plane);
            else
                assert(numel(slow_step_um) == obj.num_plane, 'Mismatch array size');
            end
            if isscalar(power_W)
                power_W = repelem(power_W, 1, obj.num_plane);
            else
                assert(numel(power_W) == obj.num_plane, 'Mismatch array size');
            end            
            if isscalar(diffuserQ)
                diffuserQ = repelem(diffuserQ, 1, obj.num_plane);
            else
                assert(numel(diffuserQ) == obj.num_plane, 'Mismatch array size');
            end                        
            
            obj.fast_v_um_s = fast_v_um_s;
            obj.fast_acc_len_um = fast_acc_len_um;
            obj.slow_step_um = slow_step_um;
            
            obj.power_W = power_W;
            obj.diffuserQ = diffuserQ;
            
            obj.bbox_xy_mmxx_um = stage_bbox_xy_um;
            obj.bbox_xy_ll_um = obj.bbox_xy_mmxx_um(3:4) - obj.bbox_xy_mmxx_um(1:2) + 1;
            
            if isempty(mask_xy_um_stack)
                % Generate mask for bounding box ablation 
                mask_xy_um_stack = true([obj.bbox_xy_ll_um, obj.num_plane]);
            end
            if obj.num_plane > 1 && size(mask_xy_um_stack, 3) == 1
                mask_xy_um_stack = repmat(mask_xy_um_stack, 1, 1, obj.num_plane);
            end
            obj.mask_xy_um_stack = mask_xy_um_stack;            
            
            obj.validate_inputs();
        end
        %%
        function validate_inputs(obj)
            if ~isempty(obj.mask_xy_um_stack)
                mask_size_2D = size(obj.mask_xy_um_stack, [1, 2]);
                assert(all(obj.bbox_xy_ll_um == mask_size_2D), 'Mask size is inconsisdent with the bounding box');
                assert(numel(obj.abl_fp_abs_z_um) == size(obj.mask_xy_um_stack, 3));
            end            
        end
        
        %%
        function path_stack = construct_path_from_mask(obj, options)
            arguments
                obj (1,1) WBIMAblationPath3D
                options.fix_start_Q (1,1) logical = true;
            end
            path_stack = cell(obj.num_plane, 1);
            start_stage_xy_um = [];
            for i = 1 : obj.num_plane
                if any(obj.mask_xy_um_stack(:, :, i), 'all')
                    path_obj = WBIMAblationPath2Dv2(obj.fast_v_um_s(i), ...
                        obj.fast_acc_len_um(i), obj.slow_step_um(i), obj.abl_fp_abs_z_um(i), ...
                        obj.power_W(i), obj.diffuserQ(i), obj.bbox_xy_mmxx_um, ...
                        obj.mask_xy_um_stack(:, :, i));
                    stripe_offset = mod(i, 2) * -0.5;
                    path_obj.coverage_path_planning('stripe_offset', stripe_offset);
                    path_obj.trajectory_i = path_obj.construct_optimal_trajectory('start_xy_um', start_stage_xy_um);
                    path_obj.trajectory = path_obj.add_pretrigger_compensation(path_obj.trajectory_i);
                    if options.fix_start_Q
                        start_stage_xy_um = path_obj.trajectory(end, 1:2);
                    end
                    path_stack{i} = path_obj;
                else
                    fprintf('Section %d in the mask is emtpy. Skip...\n', i);
                end
            end
            path_stack = cat(1, path_stack{:});
            obj.path = path_stack;
        end
    end
        
    
end