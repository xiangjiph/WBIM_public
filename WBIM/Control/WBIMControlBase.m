classdef WBIMControlBase < WBIMConfig
    
    properties
       experiment_group string = "Test"
       experiment string = sprintf("TEST_%s", datestr(now, 'YYYYMMDD'));
       current_layer (1,1) double = 1
    end    
    
    properties(SetAccess=protected)
        layer_abs_z_um (1,1) double = nan;
    end        
    
    properties(Transient, SetAccess=protected, Hidden)
        h_zaber_controller ZaberController
        file_manager WBIMFileManager = WBIMFileManager();
    end
    properties(Transient, SetAccess=protected)
        h_logger mlog.Logger
        h_notify WBIMNotification
        run_with_SI_Q (1,1) logical = false;
        hSI scanimage.SI
    end
    
    properties(Dependent, Hidden)
        microscope_log_filepath string
        log_folder string
    end
    
    properties(Access=protected)
        % Parent handle
        hWBIM WBIMControl
        
        % Piezo position and park position
        piezo_z_r_um (1,1) double
        piezo_park_z_r_um (1,1) double
    end
    %% 3-axes stage parameters
    properties(Dependent)
        stage_xyz_um (1,3) double
        sample_xyz_um (1,3) double
        stage_fp_xyz_um (1, 3) double
    end
    
    properties(Hidden)
        reset_stage_xyz_um (1,3) double
    end
    %% Zaber
    properties(Access=protected)
        zaber_unit_um % um
        zaber_unit_um_per_s % um/s
        zaber_unit_m_per_s2 % m/s2
    end
    %% Constructor
    methods
        function obj = WBIMControlBase(exp_group, exp_name, hSI)
            obj@WBIMConfig();            
            obj.experiment_group = exp_group;
            obj.experiment = exp_name;
            if nargin == 3 && ~isempty(hSI) && isvalid(hSI)
                obj.hSI = hSI;
                obj.run_with_SI_Q = true;
            else
                obj.run_with_SI_Q = false;
            end
            obj.init_logger_and_notifier();
            obj.initialize_zaber_constants();
        end
        
        function delete(obj)
            delete(obj.h_zaber_controller);
            delete(obj.h_logger);
        end
        
        function save(obj)
            % To be implemented
        end
    end
    %% Initializations
    methods(Hidden)        
        function initialize_zaber_constants(obj)
            import zaber.motion.Units;
            import zaber.motion.Measurement;
            import zaber.motion.ascii.StreamBuffer;
            
            obj.zaber_unit_um = Units.LENGTH_MICROMETRES;
            obj.zaber_unit_um_per_s = Units.VELOCITY_MICROMETRES_PER_SECOND;
            obj.zaber_unit_m_per_s2 = Units.ACCELERATION_METRES_PER_SECOND_SQUARED;
        end
        
        function init_logger_and_notifier(obj)
            [folder, ~] = fileparts(obj.microscope_log_filepath);
            if ~isfolder(folder)
                mkdir(folder);
            end
            obj.h_logger = mlog.Logger(obj.LOGGER_NAME, obj.microscope_log_filepath);
            obj.h_logger.FileThreshold = obj.LOGGER_FILE_THRESHOLD;
            obj.h_logger.CommandWindowThreshold = obj.LOGGER_COMMAND_WINDOW_THRESHOLD;
            obj.h_logger.MessageReceivedEventThreshold = obj.LOGGER_EVENT_THRESHOLD;
            obj.h_logger.write("DEBUG", "Finish initializing logger");
            obj.h_notify = WBIMNotification();
        end
        
        function init_stages(obj)
            if ~isempty(obj.hSI)
                obj.h_zaber_controller = ZaberController(obj.hSI.hMotors.hMotors{1}.hAxes);
                obj.h_logger.write("DEBUG", "Finish initializing zaber controller handle");
                obj.h_zaber_controller.set_max_speed_um_s(obj.h_zaber_controller.axis_slow_id, ...
                    obj.SLOW_AXES_MAX_SPEED_um_s);
                obj.h_zaber_controller.set_max_speed_um_s(obj.h_zaber_controller.axis_fast_id, ...
                    obj.FAST_AXES_MAX_SPEED_um_s);
                obj.h_zaber_controller.set_max_acceleration_m_per_s2(obj.h_zaber_controller.axis_slow_id, ...
                    obj.SLOW_AXIS_MAX_ACCELERATION_m_s2);
                obj.h_zaber_controller.set_max_acceleration_m_per_s2(obj.h_zaber_controller.axis_fast_id, ...
                    obj.FAST_AXIS_MAX_ACCELERATION_m_s2);
                
                obj.h_logger.write("DEBUG", "Finish setting slow axis maximum speed");
            else
                obj.h_logger.write("ERROR", "hSI is empty. Unable to initialize ZaberController");
            end
        end
        
        function init_si_related_handles(obj, hSI)
            if nargin > 1 && isvalid(hSI)
                obj.h_logger.write("MESSAGE", "Overwrite ScanImage handle");
                obj.hSI = hSI;               
                obj.run_with_SI_Q = true;
            end
            if obj.run_with_SI_Q
                obj.init_stages();
            end
        end
    end
    %% Get values for the dependent properties
    methods
        function fp = get.microscope_log_filepath(obj)
           fp = obj.file_manager.fp_microscope_log_file(obj.experiment_group, ...
               obj.experiment, obj.current_layer);
        end         
        function fp = get.log_folder(obj)
           fp = obj.file_manager.fp_log_folder(obj.experiment_group, ...
               obj.experiment, obj.current_layer);
        end        
    end
    %% Coordinate
    methods
        function val = get.stage_xyz_um(obj)
           if ~isempty(obj.hSI)
               val = obj.h_zaber_controller.get_axes_pos_um(1:3);
           else
               val = nan(1,3);
           end
        end
        
        function val = get.sample_xyz_um(obj)
            if ~isempty(obj.hSI)
                val = round(obj.stage_to_sample_xyz_um(obj.stage_xyz_um.').');
            else
                val = nan(1,3);
            end
        end
        
        function val = get.stage_fp_xyz_um(obj)
           val = obj.stage_xyz_um;
           val(3) = val(3) + obj.piezo_z_r_um;
        end
    end
    %% Coordiante and Configuration space 
    %% Layer
    methods
        function set.current_layer(obj, layer_idx)
           obj.current_layer = layer_idx;
           obj.update_layer_related_info();
        end
    end
    
    methods(Access=protected)
        function update_layer_related_info(obj)
            new_log_folder = obj.log_folder;
            if ~isfolder(new_log_folder)
                mkdir(new_log_folder)
            end            
            % TODO: Rsync the log folder to the disk 
            % Sync WBIM log file
            if ~strcmpi(obj.h_logger.LogFile, obj.microscope_log_filepath)
                try
                    % Copy the log file to HDD
                    archieve_log_fp = strrep(obj.h_logger.LogFile, ...
                        obj.file_manager.fp_acquisition_scratch,...
                        obj.file_manager.fp_acquisition_disk);
                    obj.file_manager.copy_file_local_wsl(obj.h_logger.LogFile, ...
                        archieve_log_fp, false);
                catch ME
                    obj.h_logger.write("DEBUG", "Failed to move log file to local HDD store folder");
                end                
                % Update log filepath 
                obj.h_logger.LogFile = obj.microscope_log_filepath;
            end                         
        end
    end
    %% Actuator        
    methods        
        function exit_code = actuator_switch_to_dichroic(obj)
            if ~obj.actuator_at_dichroic_position_Q
                obj.h_zaber_controller.move_to_absolute_pos_um(...
                    obj.h_zaber_controller.axis_actuator_id, ...
                    obj.ACTUATOR_DICHROIC_POS_um, true);
                obj.h_logger.write("MESSAGE", "Finish switching to the imaging dichroic");
            end
            exit_code = 0;
        end
        
        function exit_code = actuator_switch_to_mirror(obj)
            if ~obj.actuator_at_mirror_position_Q
                assert(obj.ACTUATOR_POS_MAX_um >= obj.ACTUATOR_MIRROR_POS_um); 
                obj.h_zaber_controller.move_to_absolute_pos_um(...
                    obj.h_zaber_controller.axis_actuator_id, ...
                    obj.ACTUATOR_MIRROR_POS_um, true);
                obj.h_logger.write("MESSAGE", "Finish switching to the ablation mirror");
            end
            exit_code = 0;
        end
        
        function using_mirror_Q = actuator_at_mirror_position_Q(obj)
            if obj.run_with_SI_Q
                pos_um = obj.h_zaber_controller.get_axes_pos_um(...
                    obj.h_zaber_controller.axis_actuator_id);
                using_mirror_Q = abs(pos_um - obj.ACTUATOR_MIRROR_POS_um) < ...
                    obj.ACTUATOR_POS_TOLERENCE_um;
            else
                using_mirror_Q = false;
            end
        end
        
        function using_dichroic_Q = actuator_at_dichroic_position_Q(obj)
            if obj.run_with_SI_Q
                pos_um = obj.h_zaber_controller.get_axes_pos_um(...
                    obj.h_zaber_controller.axis_actuator_id);
                using_dichroic_Q = abs(pos_um - obj.ACTUATOR_DICHROIC_POS_um) < ...
                    obj.ACTUATOR_POS_TOLERENCE_um;
            else
                using_dichroic_Q = false;
            end
        end        
    end
    %% Stage movement        
    methods
        function move_sample_along_axis_um(obj, axis_id, target_pos_um, asyncQ)
            if nargin < 4
                asyncQ = false;
            end
            assert(numel(axis_id) == numel(target_pos_um));
            next_pos_um = obj.sample_xyz_um;
            next_pos_um(axis_id) = target_pos_um;
            next_pos_um_axis = obj.sample_to_stage_xyz_um(next_pos_um.').';
            % Need to move all the axis of the stage, as the axis_id only
            % specify the axis id of the smaple space, not the stage space
            current_stage_xyz_um = obj.stage_xyz_um;
            mv_ax_ind = find(next_pos_um_axis ~= current_stage_xyz_um);   
            obj.move_stage_along_axis_um(mv_ax_ind, next_pos_um_axis(mv_ax_ind), ...
                asyncQ);
        end
        
        function move_stage_along_axis_um(obj, axis_id, target_pos_um, asyncQ)
            if nargin < 4
                % i.e. wait for one axis to finish movement before the next
                % axis starts
                asyncQ = false; 
            end
            assert(numel(axis_id) == numel(target_pos_um));
            % TODO: check configuration space here
            next_pos_um = obj.stage_xyz_um;
            next_pos_um(axis_id) = target_pos_um;          
            next_in_c_spaceQ = obj.stage_xyz_um_is_in_c_space_Q(next_pos_um);
            if next_in_c_spaceQ
                obj.h_zaber_controller.move_to_absolute_pos_um(axis_id, target_pos_um, ...
                    ~asyncQ);
            else
                err_txt = sprintf("Target position (%d, %d, %d) um is outside the configuration space", ...
                    round(next_pos_um));
                obj.h_logger.write("ERROR", err_txt);
                error(err_txt); %#ok<SPERR>
            end
            % Update SI coordiante
            obj.hSI.hMotors.queryPosition();
        end
        
        function obj = stop(obj)
            % Interrupt stage motion
            obj.h_zaber_controller.stop();
            obj.h_logger.write("MESSAGE", "Zaber stage emergency stop");
        end
        
        function move_focal_plane_to_abs_z_um(obj, target_fp_z_abs_um)
            next_fp_pos_um = obj.stage_fp_xyz_um;
            next_fp_pos_um(3) = target_fp_z_abs_um;
            % This is just an approximation 
            next_in_c_spaceQ = obj.stage_xyz_um_is_in_c_space_Q(next_fp_pos_um);
            if next_in_c_spaceQ
                dz = next_fp_pos_um(3) - obj.stage_xyz_um(3);
                if dz > obj.PIEZO_RANGE_um(1) && dz < obj.PIEZO_RANGE_um(2)
                    obj.move_piezo_to_z_r_um(dz);                    
                else
                    obj.move_piezo_to_z_r_um(0); 
                    obj.move_stage_along_axis_um(obj.h_zaber_controller.axis_z_id, next_fp_pos_um(3));
                end                
            else
                err_txt = sprintf("Target position (%d, %d, %d) um is outside the configuration space", ...
                    round(next_fp_pos_um));
                obj.h_logger.write("ERROR", err_txt);
                error(err_txt); %#ok<SPERR>
            end
        end
    end
    %% Piezo
    methods
        function move_piezo_to_z_r_um(obj, z_r_um)
            validateattributes(z_r_um, {'numeric'}, {'scalar','finite','real', ...
                '>=', obj.PIEZO_RANGE_um(1), '<=', obj.PIEZO_RANGE_um(2)});
            if obj.piezo_z_r_um ~= z_r_um
                % Only work for single peizo 
                obj.hSI.hFastZ.move(obj.hSI.hFastZ.hFastZs{1}, z_r_um);
            end
        end        
        
        function set_peizo_park_z_r_um(obj, z_r_um)
            validateattributes(z_r_um, {'numeric'}, {'scalar','finite','real', ...
                '>=', obj.PIEZO_RANGE_um(1), '<=', obj.PIEZO_RANGE_um(2)});
            if obj.hSI.hFastZ.hFastZs{1}.parkPosition ~= z_r_um
                obj.hSI.hFastZ.hFastZs{1}.parkPosition = z_r_um;
                obj.park_piezo();
                obj.h_logger.write("DEBUG", sprintf("Set peizo park position to %d um", ...
                    round(z_r_um)));
            end
            % This is redundant... but the previous park does not work
            % sometimes...
%             if abs(obj.hSI.hFastZ.position - z_r_um) > 1
%                 obj.hSI.hFastZ.move(z_r_um);
%             end
        end
        
        function park_piezo(obj)
            obj.hSI.hFastZ.hFastZs{1}.park();
        end
        
        function val = get.piezo_z_r_um(obj)
            if obj.run_with_SI_Q
                val = obj.hSI.hFastZ.position;
            else
                val = nan;
            end
        end        
        
        function val = get.piezo_park_z_r_um(obj)
            if obj.run_with_SI_Q
                val = obj.hSI.hFastZ.hFastZs{1}.parkPosition;
            else
                val = nan;
            end
        end
    end
    %% Utilities
    methods
        function record(obj, level, message)
            obj.h_logger.write(level, message);
            
        end
    end
    methods(Static)

    end
    %% Handling stage buffer memory overflow error
    methods

    end
    
    methods

    end
end