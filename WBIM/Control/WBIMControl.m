classdef WBIMControl < WBIMControlBase
    %% WBIMControl implements the high-level operations of the microscope
    % 1. Switches between modes (explore, scan, ablate)
    % 2. Keep track of all the valid imaged tiles 
    % 3. 
        
    %%
    properties(Hidden, Constant)
       COMPONENT_NAME = 'WBIMControl'; 
    end
    properties        
        imaging WBIMImaging
        ablation WBIMAblation
    end
    %% State properties
    properties
        roi_detection_channel (1, :) double = [1];
        with_skull_Q (1,1) logical = true;
        parameters struct = struct % To be updated. For multiple ablation and imaging parameters        
        
        operation_mode (1,1) WBIMOperationMode = WBIMOperationMode.Manual
        % For handling e.g. large cc 
        warning_handle_mode (1,1) WBIMOperationMode = WBIMOperationMode.Manual
        
        continuous_operation_Q (1,1) logical = false;
        mode_sequence (1, :) WBIMMicroscopeMode = [WBIMMicroscopeMode.Explore, WBIMMicroscopeMode.Scan, WBIMMicroscopeMode.Ablation];
        ablation_detection_channel (1, :) double = [1];
        
        scan_roi (1,1) struct
        async_tile_Q (1,1) logical = true;
        % Add ablation step after exploration if the albation after scan was not clean. 
        exploration_ablation_refinement_Q (1,1) logical = true; 
        enable_reverse_refinement_Q (1, 1) logical = false;
    end
    
    properties(SetAccess=private)
        tile_manager WBIMTileManager      
        current_mode (1,1) WBIMMicroscopeMode = WBIMMicroscopeMode.Unknown
        

        running_reverse_refinement_Q (1,1) logical = false; % True if triggered reverse exploration ablation refinement
        reverse_refinement_times (1, 1) double = 0;
        reverse_refinement_done_Q (1, 1) logical = false;
        reverse_refinement_triggered_layer (1, :) = []
        rescan_during_reverse_refinement_Q (1,1) logical = false;
    end
    
    properties(Dependent)
        current_state (1,1) WBIMMachineState
    end
    
    properties(Hidden)
        tile_manager_list (:, 1) WBIMTileManager
        % History
        l_layer_idx (1, :) double
        l_layer_abs_z_um (1, :) double 
        previous_microscope_mode (1,1) WBIMMicroscopeMode = WBIMMicroscopeMode.Unknown
        exploration_ablated_tile_grid_ind (1, :) double = [];
        exploration_ablation_times (1,1) double = 0;
        l_num_reverse_refinement (1, :) double
        
        latest_ablation_str = []
        
        % Image processing properties
        detect_vsl_in_skull_Q (1, 1) logical = true;
        channel_unmixing_Q (1, 1) logical = true;
    end
    
    properties(Access=private, Transient)
        state_running_computation_Q (1,1) logical = false;
        roi_info_top
        roi_info_button
    end
    %% Listeners & Timers
    properties(Hidden, SetAccess=private)
        hPropListeners = []; 
        h_step_listeners = [];
        last_stage_xyz_um (1, 3) double
        h_state_checking_timer = []
    end
    %% Gets
    properties(Hidden, Dependent)
       fp_RI_log 
    end
    properties(Access=protected, Dependent)
       running_with_ablation_Q 
    end
    %% Major methods
    methods
        function obj = WBIMControl(exp_group, exp_name, im_para, ablation_para, hSI)
            obj @ WBIMControlBase(exp_group, exp_name, hSI);
            obj.initialize_imaging(im_para);
            try
                obj.initialize_ablation(ablation_para);
            catch ME
                obj.imaging.delete();
                obj.delete();
                rethrow(ME);
            end
            obj.parameters.imaging = im_para;
            obj.parameters.ablation = ablation_para;
            % TODO: decide when to save the parameters
%             obj.save_parameters();
            try
                obj.init_si_related_handles(); % Super class method
                obj.init_microscope_state_checker();
%             obj.set_microscope_mode(WBIMMicroscopeMode.Explore);
                if obj.run_with_SI_Q
                    obj.ablation.close_ablation_shutter();
                    obj.actuator_switch_to_dichroic();
                end
            catch ME
                obj.imaging.delete();
                obj.ablation.delete();
                obj.delete();
                rethrow(MW);
            end
        end
        
        function delete(obj)
%             if obj.run_with_SI_Q
%                 obj.actuator_switch_to_mirror();
%             end
            delete(obj.imaging);
            delete(obj.ablation);
            if ~isempty(obj.h_state_checking_timer)
                obj.h_state_checking_timer.stop();
            end
            delete(obj.h_state_checking_timer);
        end
                
        function save(obj)
            % To be implemented
            % Save all the controlling parameters 
        end
        
        function save_parameters(obj, overwriteQ)
            if nargin < 2
                overwriteQ = false;
            end
            info = struct;
            save_fields = {'experiment_group', 'experiment', 'parameters', ...
                'scan_roi', 'with_skull_Q', 'current_layer', 'l_layer_abs_z_um'};
            for i = 1 : numel(save_fields)
                tmp_fn = save_fields{i};
                info.(tmp_fn) = obj.(tmp_fn);                
            end
            info.imaging.grid = obj.imaging.grid;
            info.ablation.laser_output_power_W = obj.ablation.laser_output_power_W;
            info.date = datestr(now, WBIMConfig.LOGGER_TIME_FORMAT);
            file_count = 0;
            while true
                file_path = obj.file_manager.fp_acq_control_parameter(obj.experiment_group,...
                    obj.experiment, file_count);
                if isfile(file_path) && ~overwriteQ
                    % Do not overwrite 
                    file_count = file_count + 1;
                else
                    break;
                end
            end
            folder = fileparts(file_path);
            if ~isfolder(folder)
                mkdir(folder);
            end
            save(file_path, '-struct', 'info');        
            obj.h_logger.write("MESSAGE", sprintf("Saved current control parameters to %s", ...
                file_path));
        end
    end
    
    methods(Static)
        function obj = reload_obj(exp_group, exp_name, hSI)
            if nargin < 3
                hSI = [];
            end
            DataManager = WBIMFileManager;
            ctrl_para = load(DataManager.fp_acq_control_parameter(exp_group, exp_name));
            imp = ctrl_para.parameters.imaging;
            abl_p = ctrl_para.parameters.ablation;
            with_skull_Q = ctrl_para.with_skull_Q;
            obj = WBIMControl(exp_group, exp_name, imp, abl_p, hSI);
            obj.with_skull_Q = with_skull_Q;
            obj.scan_roi = ctrl_para.scan_roi;
            obj.imaging.exploration_sample_xyz_mmxx_um = obj.exp_xyz_mmxx_um;
            for i = 1 : numel(obj.tile_manager_list)
                obj.tile_manager_list(i).load_tile_info();
            end
            obj.l_layer_idx = 1 : obj.tile_manager_list(1).num_layer;
            obj.l_layer_abs_z_um = nan(size(obj.l_layer_idx));
            for i = 1 : obj.tile_manager_list(1).num_layer
                tmp_tiles = obj.tile_manager_list(1).lc_tile{i};
                if ~isempty(tmp_tiles)
                    obj.l_layer_abs_z_um(i) = tmp_tiles(1).sample_xyz_um_done(3);
                end
            end
        end
    end
%     %% Initialization
    methods(Hidden)        
        function initialize_imaging(obj, im_para)
            if isempty(obj.imaging) || ~isvalid(obj.imaging)
                obj.imaging = WBIMImaging(obj.experiment_group, obj.experiment, ...
                    im_para, obj.hSI);
                obj.imaging.hWBIM = obj;
            else
                obj.imaging.update_grid_parameters(im_para);
            end
            % Set BG detection channel here at the moment. Could be
            % decoupled: 
            obj.roi_detection_channel = im_para.active_channel;
            obj.init_tile_manager();
        end
        
        function initialize_ablation(obj, abl_para)
            % To be changed
            obj.ablation_detection_channel = obj.roi_detection_channel;
%             obj.ablation_detection_channel = uint8([abl_para.ch_id]);
            if isempty(obj.ablation) || ~isvalid(obj.ablation)
                obj.ablation = WBIMAblation(obj.experiment_group, obj.experiment, ...
                    obj.hSI);
                obj.ablation.hWBIM = obj;
            end
        end
        
        function init_tile_manager(obj)
            if isvalid(obj.tile_manager)
                obj.tile_manager.delete();
            end
            for i = 1 : numel(obj.imaging.grid)
                obj.tile_manager_list(i) = WBIMTileManager(obj.experiment_group, ...
                    obj.experiment, obj.imaging.grid(i));
            end
        end
        
        function reconnect_to_stage_handle_in_SI(obj)
            obj.imaging.h_zaber_controller.delete();
            obj.ablation.h_zaber_controller.delete();
            obj.h_zaber_controller.delete();
            obj.hSI.hMotors.hMotors{1}.reinit;
            obj.init_stages();
            obj.imaging.init_stages();
            obj.ablation.init_stages();      
            obj.ablation.h_zaber_controller.init_timer();
        end
    end
    %% Sets and Gets
    methods
        function set.current_mode(obj, val)
            validateattributes(val, {'WBIMMicroscopeMode'}, {});
            obj.current_mode = val;
            obj.update_current_mode();
        end
        
        function set.continuous_operation_Q(obj, val)
            arguments
                obj (1,1) WBIMControl
                val (1,1) logical
            end
            if obj.continuous_operation_Q ~= val
                obj.continuous_operation_Q = val;
                if val
                    % Start idle checker timer
%                     obj.h_step_listeners.start(); %#ok<MCSUP> 
                    obj.h_logger.write("MESSAGE", "Start continuous operation");
                else
                    % Stop idle checker timer 
%                     obj.h_step_listeners.stop(); %#ok<MCSUP>
                    obj.h_logger.write("MESSAGE", "Stop continuous operation");
                end                
            end
        end
        
        function set.operation_mode(obj, val)
            arguments 
                obj (1,1) WBIMControl
                val (1,1) WBIMOperationMode
            end
            obj.operation_mode = val;
            obj.imaging.operation_mode = val; %#ok<MCSUP>
            obj.ablation.operation_mode = val; %#ok<MCSUP>
%             if obj.operation_mode == WBIMOperationMode.Automatic && ...
%                     obj.current_mode == WBIMMicroscopeMode.Scan %#ok<MCSUP>
%                 obj.imaging.si_turn_off_display(); %#ok<MCSUP>
%             else
%                 obj.imaging.si_turn_on_display(); %#ok<MCSUP>
%             end                
        end
        
        function val = get.current_state(obj)
            if obj.imaging.current_state == WBIMMachineState.Busy || ...
                    obj.ablation.current_state == WBIMMachineState.Busy
                val = WBIMMachineState.Busy;
            elseif obj.imaging.current_state == WBIMMachineState.Idle && ...
                    obj.ablation.current_state == WBIMMachineState.Idle
                val = WBIMMachineState.Idle;
            else
                val = WBIMMachineState.Unknown;
                warning('Current machine state is unknwn!');
            end            
        end
        
        function val = get.fp_RI_log(obj)
            val = obj.file_manager.fp_RI_log_file(obj.experiment_group, obj.experiment);
        end
        
        function val = get.running_with_ablation_Q(obj)
            val = any(obj.mode_sequence == WBIMMicroscopeMode.Ablation);
        end
        
    end
    methods(Access=private)
        function update_current_mode(obj)
            if ismember(obj.current_mode, [WBIMMicroscopeMode.Explore, WBIMMicroscopeMode.Scan])
                obj.tile_manager = obj.tile_manager_list(obj.current_mode);
            end
        end
        
    end
    %% Integration
    methods
        

    end
    %% Main Elementary Functions
    methods
        function set_microscope_mode(obj, new_mode, varargin)
            if obj.current_state == WBIMMachineState.Idle
                obj.previous_microscope_mode = obj.current_mode;
                obj.current_mode = new_mode;
                if ismember(new_mode, [WBIMMicroscopeMode.Explore, WBIMMicroscopeMode.Scan])
                    obj.switch_to_imaging_mode();
                    if nargin > 2
                        im_z_r_um = varargin{1};
                    else
                        im_z_r_um = [];
                    end
                    if ~isnumeric(im_z_r_um)
                        switch im_z_r_um
                            case {'top', 'Top'} % Single frame, at the top of the stack
                                im_z_r_um = obj.imaging.grid(WBIMMicroscopeMode.Scan).piezo_z_list_um(1);
                            case {'button', 'Button'} % Single frame, at the button of the stack
                                im_z_r_um = obj.imaging.grid(WBIMMicroscopeMode.Scan).piezo_z_list_um(end);
                            otherwise
                                im_z_r_um = [];
                        end
                    end
                    obj.imaging.set_acq_mode(new_mode, im_z_r_um);
                    if ~isempty(obj.scan_roi) && all(isfield(obj.scan_roi, ...
                            {'exp_xyz_mmxx_um', 'scan_xyz_mmxx_um'}))
                        switch new_mode
                            case WBIMMicroscopeMode.Explore
                                obj.imaging.exploration_sample_xyz_mmxx_um = ...
                                    obj.scan_roi.exp_xyz_mmxx_um;
                            case WBIMMicroscopeMode.Scan
                                obj.imaging.exploration_sample_xyz_mmxx_um = ...
                                    obj.scan_roi.scan_xyz_mmxx_um;
                        end
                    else
                        obj.imaging.exploration_sample_xyz_mmxx_um = [0,0,0,...
                            WBIMConfig.STAGE_LIMIT_um];
                    end
                    
                    if new_mode == WBIMMicroscopeMode.Scan
                        obj.exploration_ablation_times = 0;
                    end                    
                elseif new_mode == WBIMMicroscopeMode.Ablation
                    obj.switch_to_ablation_mode();
                end
                if obj.continuous_operation_Q
                    obj.init_step_listener(new_mode);
                    if obj.operation_mode == WBIMOperationMode.Automatic
                        if isempty(obj.h_state_checking_timer)
                            obj.init_microscope_state_checker();
                        end
                        if ~strcmp(obj.h_state_checking_timer.Running, 'on')
                            obj.h_state_checking_timer.start();
                        else
                            obj.h_state_checking_timer.stop();
                            obj.h_state_checking_timer.start();
                        end
                    end
                end
            else
                error('The microscope is busy');
            end
        end
        
        function image_candidate_tiles(obj, confirmQ)
            if obj.current_state == WBIMMachineState.Idle
                if nargin < 2
                    confirmQ = (obj.operation_mode == WBIMOperationMode.Manual);
                end
                has_valid_Q = obj.imaging.add_tile_to_acq_queue(...
                    obj.tile_manager.candidate_map_to_grid_ind(true, true));
                if has_valid_Q
                    try
                        obj.imaging.init_scanning();
                        obj.imaging.start_scanning(confirmQ);
                    catch ME
                        obj.hSI.abort();
                        rethrow(ME);
                    end
                else
                    obj.h_logger.write("MESSAGE", "No valid candidate tiles. Terminated.");
                    
                end
            else
                error('The microscope is busy');
            end
        end
        
        function stop(obj)
            obj.delete_step_listener();
            if any(obj.current_mode == [WBIMMicroscopeMode.Explore, WBIMMicroscopeMode.Scan])
                obj.imaging.stop();
            elseif obj.current_mode == WBIMMicroscopeMode.Ablation
                obj.ablation.stop();
            end
        end
        
        function ablate_region_in_mask(obj, mask_sample_yx_um, region_bbox_ctr_yx_um, ...
                abl_fp_z_abs_um, confirmQ)
            % Deprecated
            if obj.current_state == WBIMMachineState.Idle
                if nargin < 5
                    confirmQ = (obj.operation_mode == WBIMOperationMode.Manual);
                end
                % mask_sample_yx_um might be 3D, with z position specified
                % by abl_fp_z_abs_um (1 x n vector)
                abl_path_obj = obj.ablation.convert_sample_roi_to_ablation_path_xy(...
                    abl_fp_z_abs_um, region_bbox_ctr_yx_um, ...
                    mask_sample_yx_um);     
                
                is_valid_planeQ = [abl_path_obj.num_line] > 0;
                abl_path_obj = abl_path_obj(is_valid_planeQ);
                if confirmQ
                    num_obj = numel(abl_path_obj);
                    for i = 1 : num_obj
                        fig_hdl = abl_path_obj(i).visualize_mask_with_trajectory([], true);
                    end
                end
                obj.ablation.execute_ablation_paths(abl_path_obj,...
                    confirmQ);
            else
                error('The microscope is busy');
            end
        end
        
        function explore_from_current_position(obj, z_r_um)
            if obj.current_state == WBIMMachineState.Idle
                if nargin < 2
                    z_r_um = [];
                end
                obj.set_microscope_mode(WBIMMicroscopeMode.Explore, z_r_um);
                obj.add_current_pos_to_candidate_map();
                obj.tile_manager.expand_candidate_map(1, 'all');
                obj.image_candidate_tiles(obj.operation_mode ~= WBIMOperationMode.Automatic);
            else
                error('The microscope is busy');
            end
        end
        
        function explore_layer_from_cache_roi(obj, z_r_um)
            if obj.current_state == WBIMMachineState.Idle
                if nargin < 2
                    z_r_um = [];
                end    
                obj.set_microscope_mode(WBIMMicroscopeMode.Explore, z_r_um);
                [explore_tile_ind, in_mask_fraction] = obj.tile_manager.get_tile_grid_ind_in_imaged_mask(...
                    obj.roi_info_top{1}, obj.roi_info_top{2});                
                obj.tile_manager.add_tile_by_grid_ind(explore_tile_ind);
                obj.image_candidate_tiles(obj.operation_mode ~= WBIMOperationMode.Automatic);
            else
                error('The microscope is busy');
            end                
        end
        
        function scan_explored_roi(obj, local_mask_yx_um, local_tile_info, multiple_acq_Q)
            arguments
                obj (1,1) WBIMControl
                local_mask_yx_um = []
                local_tile_info = []
                multiple_acq_Q (1,1) logical = false;
            end
            confirmQ = obj.operation_mode == WBIMOperationMode.Manual;
            if obj.current_state == WBIMMachineState.Idle
                if isempty(local_mask_yx_um) && isempty(local_tile_info)
                    [local_mask_yx_um, local_tile_info] = obj.get_imaged_roi_local_mask(...
                        'channel', obj.roi_detection_channel,...
                        'acq_mode', WBIMMicroscopeMode.Explore, ...
                        'visQ', confirmQ);
                end
                if ~any(local_mask_yx_um, 'all')
                        obj.h_logger.write("MESSAGE", "The ROI mask is empty. Terminate continuous operation");
                        obj.continuous_operation_Q = false;
                else
                    obj.set_microscope_mode(WBIMMicroscopeMode.Scan);
                    [scan_tile_ind, in_mask_fraction] = obj.tile_manager.get_tile_grid_ind_in_imaged_mask(...
                        local_mask_yx_um, local_tile_info);
                    % ad hoc threshold value here
                    min_in_mask_fraction = 5e-4;
                    scan_tile_ind = scan_tile_ind(in_mask_fraction > min_in_mask_fraction);
                    if ~isempty(obj.roi_info_top) && false
                        [scan_tile_ind_2, in_mask_fraction] = obj.tile_manager.get_tile_grid_ind_in_imaged_mask(...
                            obj.roi_info_top{1}, obj.roi_info_top{2});
                        scan_tile_ind_2 = scan_tile_ind_2(in_mask_fraction > min_in_mask_fraction);
                        scan_tile_ind = union(scan_tile_ind, scan_tile_ind_2);
                    end
                    if ~isempty(scan_tile_ind)
                        obj.h_logger.write("MESSAGE", ...
                            sprintf('Add %d tiles into the acquisition queue', numel(scan_tile_ind)));
                        obj.tile_manager.add_tile_by_grid_ind(scan_tile_ind, multiple_acq_Q);
                        obj.image_candidate_tiles(confirmQ);
                    else
                        obj.h_logger.write("MESSAGE", "No candidate tile is detected. Terminated");
                    end
                    % Clear cached roi data
                    %                 obj.roi_info_top = [];                    
                end
            else
                error('The microscope is busy');
            end
        end
        
        function park_for_loading_sample(obj)
            if obj.current_state == WBIMMachineState.Idle
                obj.stop();
                % TODO: If the laser is on, close the laser
                obj.actuator_switch_to_mirror();
                if obj.ablation.pump_state
                    obj.ablation.turn_off_pump();
                end
                obj.last_stage_xyz_um = obj.stage_xyz_um;
                obj.move_stage_along_axis_um(3, 0);
                obj.move_stage_along_axis_um(1, 0);      
            else
                error('The microscope is busy');
            end
        end
        
        function return_to_the_prior_position(obj)
            obj.move_stage_along_axis_um(1:3, obj.last_stage_xyz_um);
        end
        
    end
    
    methods(Hidden)
        function switch_to_imaging_mode(obj)
            if obj.run_with_SI_Q
                obj.ablation.close_ablation_shutter();
                if obj.h_zaber_controller.get_ablation_do_value()
                    obj.h_zaber_controller.set_ablation_do_value(...
                        obj.h_zaber_controller.val_do_ablation_end);
                end
                if obj.actuator_at_mirror_position_Q
                    obj.actuator_switch_to_dichroic;
                end
                if obj.ablation.pump_state() == WBIMDeviceStateOnOff.On
                    obj.ablation.turn_off_pump();
                end
            end
            
        end
        
        function switch_to_ablation_mode(obj)
            % For safity issue, disable the laser output before switching
            % the mirror position 
            if obj.run_with_SI_Q
                obj.ablation.disable_laser_output();
                if obj.actuator_at_dichroic_position_Q
                    obj.actuator_switch_to_mirror();
                end
                if ~obj.ablation.pump_state()
                    obj.ablation.turn_on_pump();
                end                
            end
        end
    end
    %% Helper functions
    methods
        function nearest_tile_info = move_to_nearest_tile(obj, acq_mode)
            arguments
                obj WBIMControl
                acq_mode (1,1) uint8 = obj.current_mode
            end
            grid_hdl = obj.imaging.grid(acq_mode);
            nearest_tile_info = grid_hdl.get_tile_info_at_xy_um(...
                obj.sample_xyz_um(1:2));
            center_xy_um = nearest_tile_info.center_xy_um;
            obj.move_sample_along_axis_um(1:2, center_xy_um);
        end
        
        function add_current_pos_to_candidate_map(obj)
           obj.tile_manager.add_tile_at_xy(obj.sample_xyz_um(1:2));  
        end
                
        function set_scan_roi(obj, scan_sample_xyz_mmxx_um)
            % Determine scan roi: 
            % The exploration bounding box should be cover the entire scan
            % bounding box. Set the exploration bounding to obj.imaging
            % for auto-exploration. Crop the step MIP of the exploration
            % mode by scan ROI. Only refine the ablation within the scan
            % ROI
            arguments
                obj (1,1) WBIMControl
                scan_sample_xyz_mmxx_um (1, 6) double
            end
            % All the coordinates by default are in the sample space 
            roi = struct;
            roi.scan_xyz_mmxx_um = scan_sample_xyz_mmxx_um;
            roi.scan_yx_mmxx_um = roi.scan_xyz_mmxx_um([2,1,5,4]);
            scan_grid_hdl = obj.imaging.grid(WBIMMicroscopeMode.Scan);
            roi.scan_grid_ind = scan_grid_hdl.get_tiles_completely_in_bbox_um(roi.scan_yx_mmxx_um);
            % Get all the exploration tiles that overlap with the scan ROI
            exp_grid_hdl = obj.imaging.grid(WBIMMicroscopeMode.Explore);
            exp_grid_ind = exp_grid_hdl.get_tiles_overlap_with_bbox_um(...
                roi.scan_yx_mmxx_um);
            exp_tile_info = exp_grid_hdl.get_stack_info(exp_grid_ind);
            roi.exp_yx_mmxx_um = [min(exp_tile_info.tile_mmxx_um(:, [1,2]), [], 1), ...
                max(exp_tile_info.tile_mmxx_um(:, [3,4]), [], 1)];
            % For cropping the exploration smip
            
            roi.scan_wrt_exp_yx_mmxx_um = roi.scan_yx_mmxx_um - ...
                roi.exp_yx_mmxx_um([1,2,1,2]) + 1;
            assert(all(roi.scan_wrt_exp_yx_mmxx_um > 0), 'Scan bbox should be completely in the exploration bounding box');
            
            roi.exp_xyz_mmxx_um = [roi.exp_yx_mmxx_um([2,1]), roi.scan_xyz_mmxx_um(3), ...
                roi.exp_yx_mmxx_um([4,3]), roi.scan_xyz_mmxx_um(6)];
            roi.exp_grid_ind = exp_tile_info.grid_ind;
            
            obj.scan_roi = roi;
            obj.imaging.exploration_sample_xyz_mmxx_um = roi.exp_xyz_mmxx_um;
%             obj.imaging.exploration_sample_xyz_mmxx_um
%             exp_bbox_yx_mmxx_um = 
            % Compute the corresponding bounding box for the
            % exploration tiles.             
        end
    end
    %% Pipeline control
    methods(Hidden)
        function next_mode = get_next_microscope_mode(obj, current_mode)
            if nargin < 2
                current_mode = obj.current_mode;
            end
            current_idx = find(obj.mode_sequence == current_mode);
            if ~isempty(current_idx)
                next_mode_idx = mod(current_idx, numel(obj.mode_sequence)) + 1;
                next_mode = obj.mode_sequence(next_mode_idx); 
            else
                obj.h_logger.write("WARNING", sprintf("%s mode is not in the sequence", ...
                    WBIMMicroscopeMode(current_mode)));
                next_mode = WBIMMicroscopeMode(0);
            end
        end
    end
    
    methods        
        function reach_the_end_of_exploration(obj, varargin)
            obj.delete_step_listener();
            confirmQ = (obj.operation_mode == WBIMOperationMode.Manual);
            if obj.continuous_operation_Q   
                try
                    obj.state_running_computation_Q = true;                    
                    if confirmQ
                        response = questdlg('Finish exploration. Start next step?', obj.COMPONENT_NAME, ...
                            'OK', 'Cancel', 'Cancel');
                        switch response
                            case 'OK'
                                move_forward_Q = true;
                            case 'Cancel'
                                move_forward_Q = false;
                        end
                    else
                        move_forward_Q = true;
                    end                    
                    
                    if move_forward_Q
                        next_mode = obj.get_next_microscope_mode(WBIMMicroscopeMode.Explore);
                        assert(next_mode ~= WBIMMicroscopeMode.Ablation, 'To be implemented')
                        if obj.running_with_ablation_Q && obj.exploration_ablation_refinement_Q
                            abl_vol_hdl = obj.analyze_explored_tiles('visQ', false);
                            assert(next_mode == WBIMMicroscopeMode.Scan);
                            next_mode = obj.state_transition_selection_after_exploration(abl_vol_hdl);
                        end
                        
                        switch next_mode
                            case WBIMMicroscopeMode.Explore
                                if obj.enable_reverse_refinement_Q && obj.running_reverse_refinement_Q
                                    if obj.reverse_refinement_done_Q
                                        % After ablating the previous layer
                                        obj.state_transition_exploration_to_next_layer_exploration();
                                    else
                                        % Start to explore the previous
                                        % layer
                                        obj.state_transition_exploration_to_previous_layer_exploration();
                                    end
                                else
                                    % Explore-only cycle
                                    assert(isscalar(obj.mode_sequence) && (obj.mode_sequence == WBIMMicroscopeMode.Explore));
                                    obj.state_transition_to_next_layer_imaging(WBIMMicroscopeMode.Explore);
                                end
                            case WBIMMicroscopeMode.Scan
                                % Explore-scan cycle or reach the
                                % refinement time limit or confirm clean. 
                                obj.scan_explored_roi();
                                obj.exploration_ablation_times = 0;
                            case WBIMMicroscopeMode.Ablation
                                obj.exploration_ablation_times = obj.exploration_ablation_times + 1;
                                obj.state_transition_imaging_to_ablation(abl_vol_hdl, WBIMMicroscopeMode.Explore);
                            case WBIMMicroscopeMode.Unknown
                                obj.h_logger.write("WARNING", "Unknown next mode.");
                            otherwise 
                                error('Unrecognized next mode');
                        end
                    end
                    obj.state_running_computation_Q = false;
                catch ME
                    obj.state_running_computation_Q = false;
                    obj.continuous_operation_Q = false;
                    err_report = getReport(ME, 'extended', 'hyperlinks', 'off');
                    obj.h_logger.write("ERROR", err_report);
                    obj.h_notify.send_email("WBIM ERROR", err_report);
                    rethrow(ME);
                end
            else
                obj.h_logger.write("MESSAGE", "Warning: continuous_operation_Q is false");
            end
        end

        function reach_the_end_of_scan(obj, varargin)
            obj.delete_step_listener();
            confirmQ = (obj.operation_mode == WBIMOperationMode.Manual);
            if obj.continuous_operation_Q
                try
                    obj.state_running_computation_Q = true;
                    % Need to check if any tiles have been acquired
                    obj.cache_imaged_local_mask();                    
                    if confirmQ
                        response = questdlg('Start analyzing scanned tile?', obj.COMPONENT_NAME, ...
                            'Yes', 'No', 'No');
                        switch response
                            case 'Yes'
                                move_forward_Q = true;
                            case 'No'
                                move_forward_Q = false;
                        end
                    else
                        move_forward_Q = true;
                    end
                    
                    if move_forward_Q
                        next_mode = obj.get_next_microscope_mode(WBIMMicroscopeMode.Scan);
                        switch next_mode
                            case WBIMMicroscopeMode.Ablation
                                abl_vol_hdl = obj.analyze_scaned_tiles('visQ', confirmQ && false);
                                obj.state_transition_imaging_to_ablation(abl_vol_hdl, WBIMMicroscopeMode.Scan);
                            case WBIMMicroscopeMode.Explore
                                % Start exploration without ablation 
                                obj.move_to_next_layer();
                                obj.explore_layer_from_cache_roi();
                            case WBIMMicroscopeMode.Scan
                                % Scan-only cycle
                                assert(isscalar(obj.mode_sequence) && obj.mode_sequence == WBIMMicroscopeMode.Scan);
                                obj.state_transition_to_next_layer_imaging(WBIMMicroscopeMode.Scan);
                        end
                    end
                    obj.state_running_computation_Q = false;
                catch ME
                    obj.state_running_computation_Q = false;
                    obj.continuous_operation_Q = false;
                    err_report = getReport(ME, 'extended', 'hyperlinks', 'off');
                    obj.h_logger.write("ERROR", err_report);
                    obj.h_notify.send_email("WBIM ERROR", err_report);
                    rethrow(ME);                    
                end
            else
                obj.h_logger.write("MESSAGE", "Warning: continuous_operation_Q is false");
            end
        end

        function reach_the_end_of_ablation(obj, varargin)
            obj.delete_step_listener();
            confirmQ = (obj.operation_mode == WBIMOperationMode.Manual); 
            if obj.continuous_operation_Q       
                try
                    obj.state_running_computation_Q = true;
                    if confirmQ
                        response = questdlg('Finish ablation. Start next operation mode?', obj.COMPONENT_NAME, ...
                            'OK', 'Cancel', 'Cancel');
                        switch response
                            case 'OK'
                                move_forward_Q = true;
                            case 'Cancel'
                                move_forward_Q = false;
                        end
                    else
                        move_forward_Q = true;
                    end
                    if ~isempty(obj.scan_roi) && isfield(obj.scan_roi, 'scan_xyz_mmxx_um')
                        if obj.sample_xyz_um(3) > (obj.scan_roi.scan_xyz_mmxx_um(6) + 10)
                            move_forward_Q = false;
                        end
                    end
                    
                    if move_forward_Q
                        next_mode = obj.get_next_microscope_mode(WBIMMicroscopeMode.Ablation);
                        switch next_mode
                            case WBIMMicroscopeMode.Explore
                                switch obj.previous_microscope_mode % mode before the ablation mode
                                    case WBIMMicroscopeMode.Scan
                                        % Switch to the ablation mode from
                                        % the scan mode. Move on to
                                        % exploration mode.
                                        obj.move_to_next_layer();
                                        obj.explore_layer_from_cache_roi();
                                    case WBIMMicroscopeMode.Explore
                                        % Switch to the ablation mode from
                                        % the explore mode. Go back to
                                        % check ablation result. 
                                        % Do not go into scan mode directly
                                        % from here as there is no
                                        % gaurantee for the ablation
                                        % result. 
                                        obj.state_transition_ablation_to_exploration_reimaging();
                                    otherwise
                                        error('Unrecoganize previous state: %s', ...
                                            char(obj.previous_microscope_mode));
                                end
                            case WBIMMicroscopeMode.Scan
                                % This only happen in the [Scan, Abaltion]
                                % sequence
                                assert(numel(obj.mode_sequence) == 2 && ...
                                    any(obj.mode_sequence == WBIMMicroscopeMode.Scan) && ...
                                    any(obj.mode_sequence == WBIMMicroscopeMode.Ablation));
                                obj.state_transition_to_next_layer_imaging(WBIMMicroscopeMode.Scan);
                        end
                    end
                    obj.state_running_computation_Q = false;
                catch ME
                    obj.state_running_computation_Q = false;
                    obj.continuous_operation_Q = false;
                    err_report = getReport(ME, 'extended', 'hyperlinks', 'off');
                    obj.h_logger.write("ERROR", err_report);
                    obj.h_notify.send_email("WBIM ERROR", err_report);
                    rethrow(ME);
                end
            end
        end
    end
    
    methods(Hidden) % Should be private later
        function init_microscope_state_checker(obj)
            if isempty(obj.h_state_checking_timer)
                t_obj = timer();
                t_obj.Name = "WBIM Machine State Monitor";
                t_obj.BusyMode = 'drop';
                t_obj.ExecutionMode = 'fixedSpacing';
                t_obj.Period = 30; 
                t_obj.TimerFcn = {@obj.timer_microscope_state_checker, obj};
                t_obj.StartFcn = {@obj.timer_start_microscope_state, obj};
                obj.h_state_checking_timer = t_obj;                
            end            
        end
                
        function init_step_listener(obj, microscope_mode)
            obj.delete_step_listener();
            if nargin < 2
                microscope_mode = obj.current_mode;
            end
            if microscope_mode == WBIMMicroscopeMode.Explore
                obj.h_step_listeners = addlistener(obj.imaging, 'EAcqExploreDone', @obj.reach_the_end_of_exploration);
            elseif microscope_mode == WBIMMicroscopeMode.Scan
                obj.h_step_listeners = addlistener(obj.imaging, 'EAcqScanDone', @obj.reach_the_end_of_scan);
            elseif microscope_mode == WBIMMicroscopeMode.Ablation
                obj.h_step_listeners = addlistener(obj.ablation, 'EAblationModeDone', @obj.reach_the_end_of_ablation);
            end
        end
        
        function delete_step_listener(obj)
            if ~isempty(obj.h_step_listeners)
                arrayfun(@(x) x.delete(), obj.h_step_listeners, 'UniformOutput', false);
                obj.h_step_listeners = [];
            end
        end        
        
        %% state transition 
        function state_transition_to_next_layer_imaging(obj, acq_mode)
            arguments
                obj (1, 1) WBIMControl
                acq_mode (1, 1) WBIMMicroscopeMode
            end
            assert(any(acq_mode == [WBIMMicroscopeMode.Scan, WBIMMicroscopeMode.Explore]), 'acq_mode should be either scan or explore');
%             assert(numel(obj.mode_sequence) <= 2, 'This method is for [scan, ablation] or [explore, ablation] operation');
            % In Exploration - ablation mode
            if obj.imaging.auto_exploration && (acq_mode == WBIMMicroscopeMode.Explore)
                if all(obj.mode_sequence ~= WBIMMicroscopeMode.Scan)
                    % If the sequence contains the scan mode, scan image
                    % has been cached for initializing exploration
                    obj.cache_imaged_local_mask();
                end
                obj.move_to_next_layer();
                obj.explore_layer_from_cache_roi();
            else
                tile_manager_h = obj.tile_manager_list(acq_mode);
                latest_imaged_tile = tile_manager_h.get_latest_tiles_info();
                roi_tile_ind = [latest_imaged_tile.grid_ind];
                obj.move_to_next_layer();
                obj.set_microscope_mode(acq_mode);
                obj.tile_manager.add_tile_by_grid_ind(roi_tile_ind);
                obj.image_candidate_tiles();
            end
        end
        
        function state_transition_imaging_to_ablation(obj, abl_vol_hdl, current_mode)
            arguments
                obj (1, 1) WBIMControl
                abl_vol_hdl (1, 1) WBIMAVD
                current_mode (1, 1) WBIMMicroscopeMode = obj.current_mode
            end
            
            if current_mode == WBIMMicroscopeMode.Explore
                obj.exploration_ablated_tile_grid_ind = abl_vol_hdl.get_exploration_refined_tile_ind;
            end
            obj.latest_ablation_str = abl_vol_hdl;
            if ~isempty(abl_vol_hdl.ablation_path)
                obj.set_microscope_mode(WBIMMicroscopeMode.Ablation);
                obj.ablation.execute_ablation_paths(abl_vol_hdl.ablation_path);
                switch current_mode
                    case WBIMMicroscopeMode.Explore
                        msg = sprintf("Start the %d round of exploration ablation refinement", ...
                            obj.exploration_ablation_times);
                    case WBIMMicroscopeMode.Scan
                        msg = sprintf("Start ablation after the scan mode");
                end
                obj.h_logger.write("MESSAGE", msg);                
            else
                msg_txt = sprintf("The ablation mask in %s mode is empty!", ...
                    char(obj.current_mode));
                obj.h_logger.write("MESSAGE", msg_txt);
                obj.h_notify.send_email("WBIM WARNING", msg_txt);
            end
            obj.start_sync_tiles_to_server();
        end
        
        function state_transition_ablation_to_exploration_reimaging(obj)
            obj.set_microscope_mode(WBIMMicroscopeMode.Explore);
            obj.tile_manager.add_tile_by_grid_ind(...
                obj.exploration_ablated_tile_grid_ind);
            obj.image_candidate_tiles(false);
            obj.h_logger.write("MESSAGE", "Switching from ablation mode to exploration mode for reimaging the ablated tiles");
        end
        
        function state_transition_exploration_to_previous_layer_exploration(obj)
            % After the last round of exploration inspection 
            assert(obj.current_mode == WBIMMicroscopeMode.Explore);
            assert(obj.enable_reverse_refinement_Q && obj.running_reverse_refinement_Q  ...
                && ~obj.reverse_refinement_done_Q);
            % Must move before setting the mode. Setting the mode
            % recalculated the z position of the stack, which is defined
            % using the absolute position. 
            obj.move_to_previous_layer();
            obj.set_microscope_mode(WBIMMicroscopeMode.Explore);
            obj.exploration_ablation_times = 0;
            obj.tile_manager.add_tile_by_grid_ind(obj.exploration_ablated_tile_grid_ind);
            obj.image_candidate_tiles(false);
            obj.h_logger.write("MESSAGE", "Reverse refinment. Return to the previous layer and explore the ablated tiles. ");
        end        
        
        function state_transition_exploration_to_next_layer_exploration(obj)
            % After the last round of exploration inspection 
            assert(obj.current_mode == WBIMMicroscopeMode.Explore);
            assert(obj.enable_reverse_refinement_Q && obj.running_reverse_refinement_Q  ...
                && obj.reverse_refinement_done_Q);
            % Must move before setting the mode. Setting the mode
            % recalculated the z position of the stack, which is defined
            % using the absolute position. 
            obj.move_to_next_layer();
            obj.set_microscope_mode(WBIMMicroscopeMode.Explore);
            obj.exploration_ablation_times = 0;
            obj.tile_manager.add_tile_by_grid_ind(obj.exploration_ablated_tile_grid_ind);
            obj.image_candidate_tiles(false);
            obj.h_logger.write("MESSAGE", "Reverse refinment done. Return to the previous layer and explore the ablated tiles. ");
        end                        
        
        function next_mode = state_transition_selection_after_exploration(obj, abl_vol_hdl)
            next_mode = WBIMMicroscopeMode.Scan;
            if abl_vol_hdl.need_ablation_Q
                if obj.exploration_ablation_times < ...
                        obj.ABLATION_MAX_NUM_EXPLORATION_REFINEMENT_TIMES
                    % Normal refinement
                    next_mode = WBIMMicroscopeMode.Ablation;
                else % Check nonedge cc statistics
                    obj.h_logger.write("MESSAGE", "Exploration ablation refinment time reaches limit.");
                    if abl_vol_hdl.exist_surface_nonedge_cc_Q()
                        obj.h_logger.write("WARNING", "Scan ROI contain surface nonedge cc.");
                        if obj.warning_handle_mode == WBIMOperationMode.Manual
                            obj.h_notify.send_email("WBIM ACTION REQUIRED", 'Detect surface nonedge cc. Need manual operation!');
                            keyboard
%                             next_mode = WBIMMicroscopeMode.Unknown;
                        elseif obj.warning_handle_mode == WBIMOperationMode.Automatic
                            if obj.enable_reverse_refinement_Q && ...
                                    obj.reverse_refinement_times < 1 && ...
                                    ~obj.running_reverse_refinement_Q
                                % Limit to one layer at most. 
                                % State control is complicated... 
                                % Go back to the previous layer, explore,
                                % and ablate.
                                obj.state_control_set_up_reverse_refinement();
                                next_mode = WBIMMicroscopeMode.Explore;
                            else
                                obj.h_notify.send_email("WBIM WARNING", 'Detect surface nonedge cc. Next mode cannot be determined. Pause operation');
                                obj.continuous_operation_Q = false;
                                next_mode = WBIMMicroscopeMode.Unknown;
                            end
                        end
                    else
                        obj.h_logger.write("MESSAGE", "Scan ROI does not contain surface nonedge cc. Move on to imaging.");
                        if obj.enable_reverse_refinement_Q && obj.running_reverse_refinement_Q
                            if obj.reverse_refinement_done_Q
                                % next_mode = WBIMMicroscopeMode.Scan;
                            end
                        end
                    end
                end
            else
                % Clean in the inspected volume
                obj.h_logger.write("MESSAGE", "Exploration confirm clean scan sample surface.");
            end
            
            if next_mode == WBIMMicroscopeMode.Scan 
                if obj.enable_reverse_refinement_Q && obj.running_reverse_refinement_Q
                    if ~obj.reverse_refinement_done_Q                   
                        % This should only be triggered when the refinement
                        % is done in the previous layer. 
                        assert(obj.reverse_refinement_triggered_layer == (obj.current_layer + 1), ...
                            'Not in the reversed layer');
                        assert(~obj.reverse_refinement_done_Q);
                        obj.reverse_refinement_done_Q = true;
                        next_mode = WBIMMicroscopeMode.Explore;
                        % Move to the next layer and explore again. Maybe
                        % can to next layer scan directly? More risky, and
                        % what if later we want to reverse to more than 1
                        % layers? 
                    else
                        % This should only be triggered when the refinement
                        % is done and the machine return to the layer where
                        % the reverse ablation refinement started.
                        assert(obj.current_layer == obj.reverse_refinement_triggered_layer)
                        obj.state_control_turn_off_reverse_refinement();
                        % Continue with normal scan of the current layer
                    end
                end                
            end
        end
        
        function state_control_turn_off_reverse_refinement(obj)
            obj.running_reverse_refinement_Q = false;
            obj.reverse_refinement_done_Q = false;
            obj.reverse_refinement_times = 0;
            obj.reverse_refinement_triggered_layer = [];
        end        
        
        function state_control_set_up_reverse_refinement(obj)
            assert(obj.enable_reverse_refinement_Q);
            assert(~obj.running_reverse_refinement_Q);
            obj.running_reverse_refinement_Q = true;
            obj.exploration_ablation_times = 0;
            obj.reverse_refinement_done_Q = false;
            obj.reverse_refinement_times = obj.reverse_refinement_times + 1;
            obj.reverse_refinement_triggered_layer = obj.current_layer;
        end
    end    
        %% Timers related
    methods(Static)
        function timer_microscope_state_checker(t_obj, evnt, obj)
            if obj.operation_mode == WBIMOperationMode.Automatic
                switch obj.current_state
                    case WBIMMachineState.Idle
                        if ~obj.state_running_computation_Q
                            if isfinite(t_obj.UserData.last_busy_time)
                                idle_time = toc(t_obj.UserData.last_busy_time);
                                if idle_time > t_obj.UserData.maximum_idle_wait_time_s
                                    t_obj.stop();
                                    % Stop continuous operation
                                    obj.continuous_operation_Q = false;
                                    warm_msg = sprintf("WBIM has been idle for %d seconds. Stop continuous operation", ...
                                        round(t_obj.UserData.maximum_idle_wait_time_s));
                                    obj.h_logger.write("WARNING", warm_msg);
                                    obj.h_notify.send_email("WBIM WARNING", warm_msg);
                                end
                            end
                        else
                            obj.h_logger.write("INFO", "Microscope state checker: processing images");
                        end
                    case WBIMMachineState.Busy
                        t_obj.UserData.last_busy_time = tic;
                    otherwise
                        % Do nothing at the moment
                end
            end
        end
        
        function timer_start_microscope_state(t_obj, evnt, obj)
            obj.h_logger.write("DEBUG", "Microscope state checker timer start");
            info = struct;
            info.maximum_idle_wait_time_s = WBIMConfig.MACHINE_MAX_IDLE_TIME_s;
            info.last_busy_time = tic;
            t_obj.UserData = info;
        end
    end
    %% Listeners
    methods(Access=private)
        function init_listeners(obj)
            
        end        
    end
    
    %% Albation related computation
    methods                
        function abl_vol_dect = analyze_scaned_tiles(obj, options)
            arguments
                obj (1,1) WBIMControl
                options.write_vis_Q (1,1) logical = true;
                options.visQ (1,1) logical = false;
            end
            tiles_in_layer = obj.tile_manager.get_latest_tiles_info();
            [smip_cell, local_tile_info] = WBIMTileManager.get_stitched_step_mip(...
                tiles_in_layer, obj.ablation_detection_channel);
            
            abl_vol_dect = WBIMAVD(smip_cell, local_tile_info, 'acquisition_mode', ...
                WBIMMicroscopeMode.Scan, 'has_skull_Q', obj.with_skull_Q);
            abl_vol_dect.compute_single_channel_masks();
            abl_vol_dect.construct_labeled_mask('detect_vsl_in_skull_Q', ...
                obj.detect_vsl_in_skull_Q);
            if abl_vol_dect.detect_ablation_object_Q
                abl_vol_dect.construct_mask_yx_um_itp();
                abl_vol_dect.construct_merged_ablation_mask(obj.parameters.ablation);
                abl_vol_dect.ablation_path = obj.ablation.construct_ablation_path_from_AVD(...
                    abl_vol_dect);
                
                if options.write_vis_Q
                    vis_folder = obj.file_manager.fp_ablation_volume_visualization(...
                        obj.experiment_group, obj.experiment, obj.current_layer, ...
                        WBIMMicroscopeMode.Scan);
                    for ch = obj.ablation_detection_channel
                        abl_vol_dect.vis_single_channel_image_and_mask(ch, 'displayQ', ...
                            options.visQ, 'save_dir', vis_folder);
                    end
                    obj.h_logger.write("DEBUG", "Finish writing detected ablation volume visualization");
                end
            else
                obj.h_logger.write("MESSAGE", "Empty ablation mask");
            end
        end
        
        function abl_str = construct_ablation_path_from_AVD(obj, abl_vol_dect)
            arguments
                obj (1,1) WBIMControl
                abl_vol_dect (:, 1) WBIMAVD
            end
            num_abl_paras = numel(abl_vol_dect.ablation_parameter);
            abl_str = cell(num_abl_paras, 1);
            for i = 1 : num_abl_paras
                tmp_abl_para = abl_vol_dect.ablation_parameter(i);
                tmp_abl_mask = abl_vol_dect.ablation_mask{i};
                tmp_abl_z_r_um = abl_vol_dect.ablation_z_r_um{i};
                tmp_abl_fp_z_um = abl_vol_dect.abs_fp_z_um_mm + tmp_abl_z_r_um;
                abl_str{i} = obj.ablation.construct_path_for_sample_space_mask(...
                    tmp_abl_mask, abl_vol_dect.region_bbox_ctr_yx_um, tmp_abl_fp_z_um, ...
                    tmp_abl_para.site_scale_factor_fast, tmp_abl_para.site_scale_factor_slow, ...
                    tmp_abl_para.peak_fluence_J_cm2, tmp_abl_para.with_diffuser_Q);
            end
            abl_str = cat(1, abl_str{:});
            if isempty(abl_str)
               abl_str = WBIMAblationPath2Dv2.empty();
            end
            % Sort ablation path in ascending order: 
            % Concern: mechanical error of the rotational stage 
%             [~, z_idx] = sort([abl_str.abl_fp_abs_z_um], 'ascend');
%             abl_str = abl_str(z_idx);
            % Not sure if sorting is a good idea. The brain also expands
            % during ablation. Might want to minimize the time difference
            % between adjacent ablation section of the same material -
            % tissue. 
        end           
    end
        %% Exploration result analysis 
    methods
        function abl_vol_dect = analyze_explored_tiles(obj, opt)
            arguments
                obj (1,1) WBIMControl
                
                opt.diff_amp_exp_lz_um (1,1) double = 1000;
                opt.line_amp_exp_lz_um (1,1) double = 200;
                opt.line_amp_baseline (1,1) double = 1.25;
                opt.sec_range (1, :) double = [];
                
                opt.write_vis_Q (1,1) logical = true;
                opt.visQ (1,1) logical = false;
            end
            tiles_in_layer = obj.tile_manager_list(WBIMMicroscopeMode.Explore).get_latest_tiles_info();
            exp_z0_offset_um = obj.imaging.grid(WBIMMicroscopeMode.Scan).piezo_z_list_um(1) - ...
                obj.imaging.grid(WBIMMicroscopeMode.Explore).piezo_z_list_um(1);
            
            % Get step mip
            [smip_cell, local_tile_info] = WBIMTileManager.get_stitched_step_mip(...
                tiles_in_layer, obj.ablation_detection_channel, 'nanTo0Q', true);
                                    
            % Crop regions within the scan ROI
            if ~isempty(obj.scan_roi) && isfield(obj.scan_roi, 'scan_yx_mmxx_um')% DEV
                local_tile_info.local_scan_yx_mmxx_um = obj.scan_roi.scan_yx_mmxx_um - ...
                    local_tile_info.region_bbox_mm_um([1,2,1,2]);
                local_tile_info.scan_roi_smip_bbox_pxl = round(local_tile_info.local_scan_yx_mmxx_um ./ ...
                    local_tile_info.step_mip_pixel_yxz_um([1,2,1,2]));
                % Prevent overflow 
                local_tile_info.scan_roi_smip_bbox_pxl = min(max(1, local_tile_info.scan_roi_smip_bbox_pxl), ...
                    local_tile_info.step_mip_size([1,2,1,2]));                
            else

            end
            
            num_layers_above_scan_plane = round(exp_z0_offset_um / local_tile_info.step_mip_pixel_yxz_um(3));
            if isempty(opt.sec_range)
                if obj.running_reverse_refinement_Q
                    opt.sec_range = [1, num_layers_above_scan_plane];
                    % The offset is the overlap
%                     opt.sec_range = [1, num_layers_above_scan_plane * 2];
%                     opt.sec_range(2) = min(opt.sec_range(2), local_tile_info.step_mip_size(3));
%                     assert(WBIMSPAblation.is_safe_for_referse_refinement(obj.parameters.ablation), ...
%                         'The ablation parameter sets are not compatable with current implementation of the reverse ablation refinement process');
                else
                    opt.sec_range = [1, num_layers_above_scan_plane];
                end
            else
                assert(numel(opt.sec_range) == 2 && all(opt.sec_range > 0), 'Invalid sec_range');
            end            
            abl_vol_dect = WBIMAVD(smip_cell, local_tile_info, 'acquisition_mode', ...
                WBIMMicroscopeMode.Explore, 'has_skull_Q', obj.with_skull_Q, ...
                'sec_range', opt.sec_range);
            
            % To be improved here. Too many leftover vessel segments above
            % the last ablation plane at the moment...
            abl_vol_dect.compute_single_channel_masks('channel_unmixing_Q', ...
                obj.channel_unmixing_Q);
            
            abl_vol_dect.construct_labeled_mask('detect_vsl_in_skull_Q', ...
                obj.detect_vsl_in_skull_Q); 
            if abl_vol_dect.detect_ablation_object_Q
                abl_vol_dect.construct_mask_yx_um_itp();
                % Enhance the ablation here? 
                % Check CC statistics
                cc_stat = abl_vol_dect.compute_section_cc_statistics();
                % Adjust exploration ablation parameters based on cc
                % statistics
                diff_amp_exp_lz_um = opt.diff_amp_exp_lz_um / (1 + obj.exploration_ablation_times);
                line_amp_exp_lz_um = opt.line_amp_exp_lz_um / (1 + obj.exploration_ablation_times);
                line_amp_baseline = opt.line_amp_baseline ^ (1 + obj.exploration_ablation_times);
                
                
                exp_abl_p = abl_vol_dect.adjust_ablation_parameters_based_on_section_cc_stat(...
                    obj.parameters.ablation.copy(), cc_stat, 'diff_amp_exp_lz_um', diff_amp_exp_lz_um, ...
                    'line_amp_exp_lz_um', line_amp_exp_lz_um, 'line_amp_baseline', line_amp_baseline);
                
                abl_vol_dect.construct_merged_ablation_mask(exp_abl_p);
                abl_vol_dect.ablation_path = obj.ablation.construct_ablation_path_from_AVD(...
                    abl_vol_dect);
                obj.h_logger.write("MESSAGE", "Exploration found objects above the first scanning plane");                
            else
                obj.h_logger.write("MESSAGE", "Exploration confirmed clean ablation result");
            end
            
            if opt.write_vis_Q
                vis_folder = obj.file_manager.fp_ablation_volume_visualization(...
                    obj.experiment_group, obj.experiment, obj.current_layer, ...
                    WBIMMicroscopeMode.Explore);
                for ch = obj.ablation_detection_channel
                    abl_vol_dect.vis_single_channel_image_and_mask(ch, 'displayQ', ...
                        opt.visQ, 'save_dir', vis_folder, ...
                        'section_list', opt.sec_range(1) : opt.sec_range(2));
                end
                obj.h_logger.write("DEBUG", "Finish writing detected ablation volume visualization");
            end
            % Cache for debuging purpose
            obj.latest_ablation_str = abl_vol_dect;
        end
        
    end 
    %% Imaging related utility
    methods
        function [local_mask, local_tile_info] = get_imaged_roi_local_mask(obj, options)
            arguments
               obj 
               options.channel (1, :) uint8 = obj.roi_detection_channel;
               options.acq_mode (1,1) WBIMMicroscopeMode = obj.current_mode;
               options.visQ (1,1) logical = false;
            end
            if obj.current_mode ~= options.acq_mode
                tm = obj.tile_manager_list(options.acq_mode);
            else
                tm = obj.tile_manager;
            end
            num_ch = numel(options.channel);
            for i = 1 : num_ch
                if i == 1
                    % To be replaced by channel-specific ROI detection 
                    [local_mask, local_tile_info] = tm.get_roi_mask_from_tiles(...
                        tm.tile_info, options.channel(i), options.visQ);
                else
                    [tmp_mask, ~] = tm.get_roi_mask_from_tiles(...
                        tm.tile_info, options.channel(i), options.visQ);
                    local_mask = local_mask | tmp_mask;
                end
            end            
        end
    end
    %% Layer related
    methods        
        function set_current_layer_abs_z_um(obj, abs_z)
            if nargin < 2
                abs_z = obj.sample_xyz_um(3);
            end
            obj.layer_abs_z_um = abs_z;
            obj.l_layer_abs_z_um(obj.current_layer) = abs_z;
            obj.imaging.layer_abs_z_um = abs_z;
            obj.ablation.layer_abs_z_um = abs_z;
        end
        
        function move_to_next_layer(obj, delta_z_um)
            if obj.current_state == WBIMMachineState.Idle
                if nargin < 2
                    delta_z_um = obj.imaging.grid(WBIMMicroscopeMode.Scan).stack_size_um(3) - ...
                        obj.imaging.grid(WBIMMicroscopeMode.Scan).stack_overlap_um(3);
                end
                next_z_abs_um = obj.layer_abs_z_um + delta_z_um;
                next_sample_xyz_um = obj.sample_xyz_um;
                next_sample_xyz_um(3) = next_z_abs_um;
                if obj.sample_xyz_um_is_in_c_space_Q(next_sample_xyz_um)
                    obj.move_sample_along_axis_um(1:3, next_sample_xyz_um);
                    obj.current_layer = obj.current_layer + 1;
                    obj.set_current_layer_abs_z_um();
                else
                    error("The target position is outside the configuration space");
                end
            else
               error('The microscope is busy'); 
            end
        end
        
        function move_to_previous_layer(obj)
            if obj.current_layer >= 2
                obj.move_to_layer(obj.current_layer - 1);
            else
                fprintf('Already in the first layer\n');
            end
        end
        
        function move_to_layer(obj, idx)
            if obj.current_state == WBIMMachineState.Idle
                next_sample_xyz_um = obj.sample_xyz_um;
                next_sample_xyz_um(3) = obj.l_layer_abs_z_um(idx);
                if obj.sample_xyz_um_is_in_c_space_Q(next_sample_xyz_um)
                    obj.move_sample_along_axis_um(1:3, next_sample_xyz_um);
                    obj.current_layer = idx;
                    obj.set_current_layer_abs_z_um();
                end
            else
                error('The microscope is busy');
            end
        end
        
        function cache_imaged_local_mask(obj, channel, acq_mode)
            arguments
                obj (1,1) WBIMControl
                channel (1, :) = obj.roi_detection_channel
                acq_mode WBIMMicroscopeMode = obj.current_mode
            end
            [local_mask_yx_um, local_tile_info] = obj.get_imaged_roi_local_mask(...
                'channel', channel, 'acq_mode', acq_mode);
            obj.roi_info_top = {local_mask_yx_um, local_tile_info};
        end
        
    end
    
    methods(Access=protected)
        function update_layer_related_info(obj)
            % Overwrite superclass method called upon set.current_layer
            obj.ablation.current_layer = obj.current_layer;
            obj.imaging.current_layer = obj.current_layer; 
            for i = 1 : numel(obj.tile_manager_list)
                obj.tile_manager_list(i).layer = obj.current_layer;
            end
        end
    end
    %% Utilities    
        %% Recovering
    methods
        function obj = reload_tiles_from_disk(obj)
            for i = 1 : numel(obj.tile_manager_list)
                obj.tile_manager_list(i).load_tile_info();
            end           
        end
        
        function obj = reset_stage_controller(obj)
            % Deinit the stage from SI handle
            if obj.run_with_SI_Q
                obj.prepare_xy_stage_for_system_reset();
                obj.hSI.hMotors.hMotors{1}.deinit();
                obj.h_logger.write("DEBUG", "SI releases Zaber handle");
            end            
            import zaber.motion.ascii.Connection;
            zaber_hdl = ZaberController();
            try
                zaber_hdl.h_all_axes.park();
                zaber_hdl.connections.genericCommand('system reset');
                zaber_hdl.delete();
                obj.h_logger.write("DEBUG", "Reset Zaber controller");
                pause(5);
                zaber_hdl = ZaberController();
                obj.h_logger.write("DEBUG", "Reconnect to the Zaber controller");
                % TODO: Confirm that the z stages and the actuator are
                % parked
                zaber_hdl.h_all_axes.unpark();      
                obj.h_logger.write("DEBUG", "Unpark stages");            
                % Home the xy stages here? 
%                 zaber_hdl.home(); % - working
                zaber_hdl.home([WBIMConfig.STAGE_AXIS_SLOW_ID, ...
                    WBIMConfig.STAGE_AXIS_FAST_ID])
                obj.h_logger.write("DEBUG", "Home xy stages");      
                zaber_hdl.delete();
                obj.h_logger.write("DEBUG", "Zaber controller handle is released successfully");
            catch exception
                zaber_hdl.delete();
                rethrow(exception);
            end
            % Reinit 
            if obj.run_with_SI_Q
                obj.hSI.hMotors.hMotors{1}.reinit();
                obj.reconnect_to_stage_handle_in_SI();
                obj.h_logger.write("DEBUG", "SI reconnect to Zaber handle");
                assert(all(obj.h_zaber_controller.isHomed), 'Some stages have not been homed yet');
                obj.return_xy_stage_after_sytem_reset();
                obj.h_logger.write("DEBUG", "Return to the position before stage controller reset.");
            end
        end        
        
        function obj = prepare_xy_stage_for_system_reset(obj, homeZQ)
            if nargin < 2
                homeZQ = true;
            end
            % Move the xy stages to a position greater than the home index
            % position if allowed
            % In stage coordiante
            obj.reset_stage_xyz_um = obj.stage_xyz_um;
            if homeZQ
                obj.move_stage_along_axis_um(3, 0);
            end
            reset_pos_um = [WBIMConfig.XY_STAGE_HOME_MIN_POSITION_um, obj.stage_xyz_um(3)];
            assert(obj.stage_xyz_um_is_in_c_space_Q(reset_pos_um), ...
                "Target reset position is not in the configuration space");
            obj.move_stage_along_axis_um(1:2, reset_pos_um(1:2));
            obj.h_logger.write("MESSAGE", "Move xy stages to reset position");
        end
        
        function obj = return_xy_stage_after_sytem_reset(obj)
            % This is assuming the xy stages do not move during system
            % reset.
            % TODO: Need to renew the zaber controller handle first
            obj.h_logger.write("MESSAGE", "Finish homing xy stages");
            assert(obj.stage_xyz_um_is_in_c_space_Q(obj.reset_stage_xyz_um), ...
                "The preset position is outside the configuration space");
            obj.move_stage_along_axis_um(1:3, obj.reset_stage_xyz_um);
            obj.h_logger.write("MESSAGE", "Move xyz stages back to the position before system reset");
        end
    end
    
    methods(Static)                
%         function save_tile_manager(obj, acq_mode)
%             validateattributes(acq_mode, 'WBIMMicroscopeMode', {'scalar'});
%             assert(any(acq_mode == [WBIMMicroscopeMode.Explore, WBIMMicroscopeMode.Scan]));
%             if nargin < 2
%                 acq_mode = obj.current_mode;
%             end
%             if isempty(obj.file_manager)
%                 obj.tile_manager = WBIMFileManager;
%             end            
%             fp = obj.file_manager.fp_tile_manager(obj.experiment_group, ...
%                 obj.experiment, char(acq_mode));
%             builtin('save', fp, obj.tile_manager_list(acq_mode));            
%         end        
    end
    %% Emergency shutdown
    % Park the stage: allAxes.park() Parks the device in anticipation of
    % turning the power off. It can later be powered on, unparked, and
    % moved without first having to home it.
    
    %% File IO 
    methods
        function start_sync_tiles_to_server(obj)
            if obj.async_tile_Q
                folder_pair = obj.tile_manager.get_tile_sync_folder_pairs();
                obj.imaging.h_server.add_sync_folders(folder_pair);
                obj.imaging.h_server.running_sync_Q = true;
            end
        end        
    end
end