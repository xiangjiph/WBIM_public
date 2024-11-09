classdef WBIMImaging < WBIMControlBase
    properties(Hidden, Constant)
        COMPONENT_NAME = 'WBIMImaging';
    end
    %% Acquisition parameters
    properties
        imaging_wavelength_nm (1,1) double = 950
        auto_exploration (1,1) logical = true;        
        auto_scan_qc_sync_Q (1,1) logical = true;
        auto_scan_qc_async_Q (1,1) logical = false;
        exploration_sample_xyz_mmxx_um (1, 6) double
        
        overwrite_file_Q (1,1) logical = false;
        active_channel (1, :) double % Set method
        bg_detection_channel (1, :) double
        im_qc_channel (1, :) double = [1,2];
        operation_mode (1,1) WBIMOperationMode = WBIMOperationMode.Manual
        image_z_r_um (:, 1) double % Set method
        target_imaging_power_mW (1,1) double = 40;
        imaging_power_profile_list (1, :) cell
    end
    
    properties(Dependent)
        imaging_power_mW
    end
    
    properties(Hidden)
        max_reimage_time (1, 1) double = 1;
        state_finish_initialization_Q (1,1) logical
        mamual_pause_Q (1,1) logical = false;
    end
    %% Directory
    properties(Dependent, Hidden)
        scratch_root_folder char
        save_root_folder char
        laser_log_filepath char

        current_scratch_folder char
        current_save_folder char        
    end

    properties(NonCopyable, Transient, Access=private)
        scratch_root_folder_
        save_root_folder_
    end
    %% Laser power control
    properties(Transient)
        h_laser LaserStateMonitor
        h_hwp WBIMPowerHWP
        h_server WBIMTCPServer = WBIMTCPServer.empty();
    end
    %% Handles
    properties(Access=private)
       fig_hdl
       ax_hdl
       im_hdl
    end
    %% Pipeline state properties
    properties(SetAccess=protected)
        current_mode (1,1) WBIMMicroscopeMode = WBIMMicroscopeMode.Unknown
        
        
        current_grid_ind (1,1) double = 0;
        current_grid_xy_um (1, 2) double
        scan_deque = java.util.ArrayDeque;
        acq_count (1,1) double = 0;
        reimage_current_tile_Q = false;
        reimage_count (1,1) double = 0;        
%         acq_piezo_z_r_um (1,1) double
        init_laser_output_power_mW (1, 1) double = nan;
    end
    
    properties(Hidden)
        current_state (1,1) WBIMMachineState
        grid (:, 1) WBIMImagingGrid2D
        current_grid WBIMImagingGrid2D
        current_tile_info WBIMTileMetadata
    end
    %% Pipeline control
    properties(Hidden) % These should all be private later
        % [Listener] Array of Listeners for controlling tile scanning - Specifically done and abort calls.
        hScanControlListeners = [];
        h_listener_laser_monitor = [];
        h_listener_server = [];
        
        h_state_timer = [];
        h_delay_timer = [];
        
    end
    
    properties(SetObservable, Transient, Hidden)
        normalAcqDone (1,1)logical = false;
        scanning_in_progress_Q (1,1)logical = false;
        continuous_acquisition_Q (1,1) logical = true;
        
        % For lineshift detection and correction
        need_restart_for_abnormal_hardwares_Q (1,1) logical = false;
    end
    %% Events
    events
        EAcqDone % Fires when the acquisition is completed.
        EAcqScanDone % Fires when the scanning acquisition is completed.
        EAcqExploreDone % Fires when the explore mode is completed
    end
    %% Life cycle
    methods
        function obj = WBIMImaging(exp_group, exp_name, im_para, hSI)
            obj @ WBIMControlBase(exp_group, exp_name, hSI);
            obj.update_grid_parameters(im_para);
            obj.initialization();
            obj.current_state = WBIMMachineState.Idle;
        end
        
        function delete(obj)
            if obj.run_with_SI_Q
                obj.h_logger.write("MESSAGE", "Resetting imaging HWP to 50%");
                obj.h_hwp.set_fraction_power(0.5);
            end
            obj.sync_server_log_file_to_disk();
            delete@WBIMControlBase(obj);
            most.idioms.safeDeleteObj(obj.hScanControlListeners);
            delete(obj.h_listener_laser_monitor);
            delete(obj.h_state_timer);
            delete(obj.h_hwp);
            delete(obj.h_laser);            
            delete(obj.h_server);
        end
    end
    %% Main functions
    methods
        function stop(obj)
           stop@WBIMControlBase(obj);
           obj.hSI.abort();  
           obj.normalAcqDone = false;
           obj.h_state_timer.stop();
           obj.h_logger.write("WARNING", "Acquisition interrupted");
        end
        
        function start_scanning(obj, ask_for_confirmation_Q)
            arguments
                obj
                ask_for_confirmation_Q (1, 1) logical = true;
            end
            
            if obj.state_finish_initialization_Q
                if ask_for_confirmation_Q
                    response = questdlg('Start scanning?', obj.COMPONENT_NAME, ...
                        'OK', 'Cancel', 'Cancel');
                    switch response
                        case 'OK'
                            obj.scan_next_tile();
                        case 'Cancel'
                            obj.hSI.abort();
                            return;
                    end
                else
                    obj.scan_next_tile();
                end
            else
                obj.h_logger.write("WARNING", "Scanning has not been properly initialized. Abort");
            end
        end
        
        function resume_scanning(obj, confirmQ)
            arguments
                obj
                confirmQ (1,1) logical = obj.operation_mode == WBIMOperationMode.Manual;
            end
            obj.init_scanning(true, false);
            obj.start_scanning(confirmQ);            
        end
        
    end    
    %% Initialize parameters for scanning mode and explore mode
    methods
        function obj = initialization(obj)
            obj.exploration_sample_xyz_mmxx_um = [0,0,0, obj.STAGE_LIMIT_um];
            if strcmpi(obj.file_manager.HOSTNAME, "PIA") && obj.run_with_SI_Q
                % Initialize hardware handles
                try
                    obj.init_imaging_laser();
                    obj.init_imaging_HWP();
                    if ~isempty(obj.hSI)
                        obj.init_si_related_handles();
                    else
                        obj.h_logger.write("WARNING", "Invalid hSI handle. Skip SI-related initialization");
                    end
                catch ME
                    obj.h_logger.write("ERROR", getReport(ME, 'extended', 'hyperlinks', 'off'));
                end
                
                try
                    obj.h_laser.set_GDD_fs2(WBIMConfig.IMAGING_LASER_GDD);
                    obj.h_laser.start_recording();
                catch ME
                    obj.h_logger.write("ERROR", getReport(ME, 'extended', 'hyperlinks', 'off'));
                end                
            end
            obj.init_server();
%             obj.set_acq_mode(WBIMMicroscopeMode.Explore);
        end
        
        function obj = init_si_related_handles(obj, hSI)
            if nargin < 2
                init_si_related_handles@WBIMControlBase(obj);
            else
                init_si_related_handles@WBIMControlBase(obj, hSI);
            end
            obj.init_timer();
            obj.si_set_channel_properties();
            obj.hSI.hScan2D.flybackTimePerFrame = obj.SI_FRAME_FLYBACK_TIME_s;
            try
                obj.hSI.hFastZ.hFastZs{1}.parkPosition = obj.PIEZO_REFERENCE_POS_um;
                obj.park_piezo();
            catch ME
                obj.h_logger.write("WARNING", 'Fail to set the piezo park position. Error message: %s', ...
                    getReport(ME, 'extended', 'hyperlinks', 'off'));
            end
%             obj.hSI.hFastZ.discardFlybackFrames = true;
%             obj.hSI.hFastZ.flybackTime = obj.SI_PIEZO_FLYBACK_TIME_s;
            obj.h_logger.write("MESSAGE", "Finish initializing scanImage");
        end
        
    end
    
    methods(Hidden)
        function update_grid_parameters(obj, im_para)
            % Constructing grid for image acquisition and fast stitching. 
            % These grids are in sample coordiante:
            %   X: positive direction of resonant scanner (fast axis),
            %   horizontal
            %   Y: positive direction of the galvo (slow axis), vertical
            obj.construct_2D_grid(WBIMMicroscopeMode.Scan, im_para.scan.stack_xyz_size, ...
                im_para.scan.pixel_xyz_size_um, im_para.scan.overlap_xyz_size, ...
                im_para.scan.z_offset_um);
            obj.construct_2D_grid(WBIMMicroscopeMode.Explore, im_para.explore.stack_xyz_size,...
                im_para.explore.pixel_xyz_size_um, im_para.explore.overlap_xyz_size, ...
                im_para.explore.z_offset_um);
            obj.update_imaging_power_parameters(im_para);
            if isfield(im_para, 'active_channel')
                obj.active_channel = im_para.active_channel;
            end
            if isfield(im_para, 'bg_detection_channel')
                obj.bg_detection_channel = im_para.bg_detection_channel;
            end
            if isfield(im_para, 'im_qc_channel')
                obj.im_qc_channel = im_para.im_qc_channel;
            end
        end
        
        function update_imaging_power_parameters(obj, im_para)
            if isfield(im_para.scan, 'imaging_power')
                obj.imaging_power_profile_list{WBIMMicroscopeMode.Scan} = im_para.scan.imaging_power;
                obj.imaging_power_profile_list{WBIMMicroscopeMode.Explore} = im_para.explore.imaging_power;
            end            
        end
        
        function init_imaging_laser(obj)
            try
                obj.h_laser = LaserStateMonitor(...
                    obj.IMAGING_LASER_PORT_ID, ...
                    obj.IMAGING_LASER_BAUD_RATE, obj.imaging_wavelength_nm, ...
                    obj.ILSM_PERIOD_S, obj.ILSM_BUFFER_SIZE, ...
                    obj.laser_log_filepath);
                obj.h_laser.normal_min_power_mW = WBIMConfig.IMAGING_LASER_MIN_NORMAL_POWER_mW;
                
                
                if ~isempty(obj.h_listener_laser_monitor)
                    delete(obj.h_listener_laser_monitor)
                end
                obj.h_logger.write("MESSAGE", "Finish initializing laser");
            catch ME
                obj.h_logger.write("ERROR", getReport(ME, 'extended', 'hyperlinks', 'off'));
                obj.h_logger.write("Info", "Fail to initialize the imaging laser");
            end
        end
             
        function init_imaging_HWP(obj)
           obj.h_hwp = WBIMPowerHWP(obj.IMAGING_HWP_STAGE_SERIAL_NUMBER, ...
               obj.IMAGING_HWP_DELTA_THETA_deg, ...
               obj.IMAGING_HWP_DEFAULT_ENERGY_FRACTION);
           obj.h_logger.write("MESSAGE", "Finish initializing imaging laser half-wave plate");
        end
        
        function obj = si_set_channel_properties(obj)
            obj.hSI.hChannels.channelSave = obj.active_channel;
            obj.hSI.hDisplay.volumeDisplayStyle = 'Current';
%             obj.hSI.hChannels.channelSubtractOffset = false(1, 4);            
            obj.hSI.hChannels.channelSubtractOffset = true(1, 4);            
            obj.hSI.hChannels.channelMergeColor{2} = 'gray';
            obj.hSI.hChannels.channelMergeColor{3} = 'blue';
            obj.hSI.hChannels.channelMergeColor{4} = 'green';
        end        
        
    end
    %% Server
    methods
        function obj = init_server(obj)
            log_filepath = obj.file_manager.fp_tcp_server_log_file(...
                obj.experiment_group, obj.experiment);
            obj.h_server = WBIMTCPServer(log_filepath);
            obj.h_server.launch_client();
        end
        
        function restart_server_and_client(obj)
            obj.h_server.shut_down_client();
            obj.h_server.delete();
            obj.init_server();           
        end
    end
    %% Set and Get
    methods
        function fp = get.scratch_root_folder(obj)
            if isempty(obj.scratch_root_folder_)
                obj.scratch_root_folder_ = fullfile(obj.file_manager.fp_acquisition_scratch, ...
                    obj.file_manager.fpr_experiment(obj.experiment_group, ...
                    obj.experiment));
            end
            fp = obj.scratch_root_folder_;
        end
        
        function fp = get.save_root_folder(obj)
            if isempty(obj.save_root_folder_)
                obj.save_root_folder_ = fullfile(obj.file_manager.fp_acquisition_disk, ...
                    obj.file_manager.fpr_experiment(obj.experiment_group, ...
                    obj.experiment));
            end
            fp = obj.save_root_folder_;
        end
        
        function fp = get.current_scratch_folder(obj)
            fp = fullfile(obj.file_manager.fp_acquisition_scratch, ...
                obj.file_manager.fpr_layer(obj.experiment_group, ...
                obj.experiment, obj.current_layer), char(obj.current_mode));
        end
        
        function fp = get.current_save_folder(obj)
            fp = fullfile(obj.file_manager.fp_acquisition_disk, ...
                obj.file_manager.fpr_layer(obj.experiment_group, ...
                obj.experiment, obj.current_layer), char(obj.current_mode));
        end
        
        function fp = get.laser_log_filepath(obj)
            fp = obj.file_manager.fp_imaging_laser_log_file(obj.experiment_group, ...
                obj.experiment, obj.current_layer);
        end
        
        function set_imaging_wavelength_nm(obj, val)
            obj.h_laser.set_wavelength_nm(val);
            obj.imaging_wavelength_nm = val;
        end
        
        function set.active_channel(obj, val)
           obj.active_channel = val;
           if obj.run_with_SI_Q
               obj.hSI.hChannels.channelSave = val;
               obj.hSI.hChannels.channelDisplay = val;
           end
        end
        
        function set.image_z_r_um(obj, val)
            assert(all(val >= obj.PIEZO_RANGE_um(1)) && ...
                all(val <= obj.PIEZO_RANGE_um(2)), 'The relative z-position to be imaged is out of range');
            obj.image_z_r_um = val;
            if obj.run_with_SI_Q
               obj.si_set_stack_manager(obj.image_z_r_um); 
            end
        end    
        
        function val = get.imaging_power_mW(obj)
            if ~isempty(obj.h_laser)
                val = obj.h_laser.power_mW * obj.hSI.hBeams.powerFractions * ...
                    obj.IMAGING_LASER_POWER_TRANSMISSION_RATE * ...
                    obj.h_hwp.fractional_power;
            else
                val = nan;
            end
        end
        
    end
    %% ScanImage
    methods
        function obj = si_set_to_focus_mode(obj)
            obj.hSI.hStackManager.enable = false;
            obj.hSI.extTrigEnable = false;
        end
        
        function si_set_beam_adjust(obj, method, varargin)
            % If method = 'exp', varargin = {fractional_intensity_z0, Lz_um};
            switch method
                case 'none'
                    if obj.hSI.hBeams.pzAdjust ~= scanimage.types.BeamAdjustTypes.None
                        obj.hSI.hBeams.pzAdjust = scanimage.types.BeamAdjustTypes.None;
                    end
                    if obj.hSI.hBeams.powerFractions ~= varargin{1}
                        obj.hSI.hBeams.powerFractions = varargin{1};
                    end
                case 'exp'
                    if obj.hSI.hBeams.pzAdjust ~= scanimage.types.BeamAdjustTypes.Exponential
                        obj.hSI.hBeams.pzAdjust = scanimage.types.BeamAdjustTypes.Exponential;
                    end
                    if obj.hSI.hBeams.powerFractions ~= varargin{1}
                        obj.hSI.hBeams.powerFractions = varargin{1};
                    end
                    if obj.hSI.hBeams.lengthConstants ~= varargin{2}
                        obj.hSI.hBeams.lengthConstants = varargin{2};
                    end
            end            
        end
        
        function si_turn_off_display(obj)
            obj.si_set_display_ch([]);
        end
        
        function si_turn_on_all_availabel_display(obj)
            obj.si_set_display_ch(obj.active_channel);
        end
        
        function si_set_display_ch(obj, disp_list)
            arguments
                obj WBIMImaging
                disp_list (1, :) double = obj.active_channel
            end
            if isempty(disp_list) && ~isempty(obj.hSI.hChannels.channelDisplay)
                % Don't use []. The get.channelsDisplay() in Acquisition.m
                % checks for assert(isequal(val,unique(val))), while the
                % unique function by default returns a row vector. Notice
                % that isequal([], unique([])) = false!
                obj.hSI.hChannels.channelDisplay = zeros([0, 1]);
            else
                assert(isempty(setdiff([1,2], 1:4)), 'Display channel list out of range');
                if numel(obj.hSI.hChannels.channelDisplay) ~= numel(disp_list) || ...
                        any(obj.hSI.hChannels.channelDisplay ~= disp_list)
                    obj.hSI.hChannels.channelDisplay = disp_list;
                end
            end            
        end
        
        
    end    
    methods(Hidden)
        %% SI initialization        
        function obj = si_set_save_dir(obj, folder_name, file_name_prefex, file_count)
            if nargin < 4
                file_count = 0;
            end
            % Somehow these two properties have to be char
            if isstring(folder_name)
                folder_name = char(folder_name);
            end
            if isstring(file_name_prefex)
                file_name_prefex = char(file_name_prefex);
            end
            
            if ~strcmpi(obj.hSI.hScan2D.logFilePath, folder_name)
                obj.hSI.hScan2D.logFilePath = folder_name;
            end
            if ~isfolder(folder_name)
                mkdir(folder_name);
            end
            if ~strcmpi(obj.hSI.hScan2D.logFileStem, file_name_prefex)
                obj.hSI.hScan2D.logFileStem = file_name_prefex;
            end
            if obj.hSI.hScan2D.logFileCounter ~= file_count
                obj.hSI.hScan2D.logFileCounter = file_count;
            end
        end
                
        function fp = si_get_current_stack_filepath(obj, file_counter)
            if nargin < 2
                file_counter = obj.hSI.hScan2D.logFileCounter;
            end
            if isfinite(obj.hSI.hScan2D.logFramesPerFile)
                fn = obj.file_manager.fn_si_tile_name_fix_num_frames(...
                    obj.hSI.hScan2D.logFileStem, file_counter);
            else
                fn = obj.file_manager.fn_si_tile_name(obj.hSI.hScan2D.logFileStem, ...
                    file_counter);
            end
            fp = fullfile(obj.hSI.hScan2D.logFilePath, fn);
        end
        %% Set nonmROI imaging
        function obj = si_set_non_mroi_imaging(obj, scan_angle_deg, im_size)
            % scan angle = [angle_fast, angle_slow], which by default
            % is the x and y axis of the sample coordinate in ScanImage
            
            % Set 2D acquisition parameters in scanimage. 
            % Input: 
            %   scan_angle_deg = [fast scanner, slow scanner]
            %   im_size = [num_row, num_column]
            if obj.hSI.hRoiManager.mroiEnable
                obj.hSI.hRoiManager.mroiEnable = false;
            end
            % Set up scan field
            if obj.hSI.hRoiManager.forceSquarePixelation
                obj.hSI.hRoiManager.forceSquarePixelation = false; % If true, num_row = num_col
            end
            if obj.hSI.hRoiManager.pixelsPerLine ~= im_size(2)
                obj.hSI.hRoiManager.pixelsPerLine = im_size(2);
            end
            if obj.hSI.hRoiManager.linesPerFrame ~= im_size(1)
                obj.hSI.hRoiManager.linesPerFrame = im_size(1);
            end
            % Compute zoom factors            
            zoom_factor_fs = round(obj.SCANNER_RANGE_DEG .* obj.hSI.hScan2D.fillFractionSpatial ./ ...
                scan_angle_deg, 2);
            assert(all(zoom_factor_fs >=1), 'Scan angle is too large'); % Well not really as here we are limited by the galvo. Can be % loosen later. TODO.
            [overall_zoom_factor, min_idx]= min(zoom_factor_fs);
            if obj.hSI.hRoiManager.scanZoomFactor ~= overall_zoom_factor
                obj.hSI.hRoiManager.scanZoomFactor = overall_zoom_factor;
            end
            min_scan_angle = obj.SCANNER_RANGE_DEG(min_idx) / obj.hSI.hRoiManager.scanZoomFactor;
            % Deal with internal rounding error - why does scanimage only
            % keep two digits?
            zoom_mag_fs = scan_angle_deg / obj.hSI.hScan2D.fillFractionSpatial / min_scan_angle;
%             zoom_mag_fs = obj.hSI.hRoiManager.scanZoomFactor ./ zoom_factor_fs;
            if obj.hSI.hRoiManager.forceSquarePixels
                obj.hSI.hRoiManager.forceSquarePixels = false;
            end
            if obj.hSI.hRoiManager.scanAngleMultiplierFast ~= zoom_mag_fs(1)
                obj.hSI.hRoiManager.scanAngleMultiplierFast = zoom_mag_fs(1);
            end
            if obj.hSI.hRoiManager.scanAngleMultiplierSlow ~= zoom_mag_fs(2) 
                obj.hSI.hRoiManager.scanAngleMultiplierSlow = zoom_mag_fs(2);
            end            
            obj.h_logger.write("DEBUG", "Finish setting RoiManager");
        end
        
        function obj = si_set_stack_manager(obj, z_rel_um)
            % Z-position in the piezo coordinate (w.r.t. 0, not the park
            % position)
            validateattributes(z_rel_um, {'numeric'}, {'vector'});
            % For single axes piezo, the relative z-position must be a
            % column vector.
            if isrow(z_rel_um)
                z_rel_um = z_rel_um.';
            end
            % General setting
            if ~strcmpi(obj.hSI.hStackManager.stackMode, 'fast')
                obj.hSI.hStackManager.stackMode = 'fast';
            end
            if ~strcmpi(obj.hSI.hStackManager.stackFastWaveformType, 'step')
                % For arbitrary fast stacks, the waveform type must be set
                % to 'step'. Need to change... ask SI
                obj.hSI.hStackManager.stackFastWaveformType = 'step';
            end
            if ~strcmpi(obj.hSI.hStackManager.stackDefinition, 'arbitrary')
                obj.hSI.hStackManager.stackDefinition = 'arbitrary';
            end
            obj.hSI.hStackManager.framesPerSlice = 1;            
            % This ROI handle is the same for all the tiles in the same
            % 2D grid on the same layer Does scanimage move the z-stage
            % based on the arbitraryZs here?
            z_rel_um = unique(z_rel_um, 'sorted');% ascending
            current_z_um = obj.sample_xyz_um(3);
            z_abs_um = current_z_um + z_rel_um;
            num_z = numel(z_abs_um);
            obj.hSI.hStackManager.numSlices = num_z;
            obj.move_piezo_to_z_r_um(z_rel_um(1));
            % arbitraryZs is the absolute z-position in microm.
            obj.hSI.hStackManager.arbitraryZs = z_abs_um;            
            switch obj.current_mode
                case WBIMMicroscopeMode.Scan
                    vol_fly_back_time_s = obj.IMAGING_SCAN_VOLUME_FLYBACK_TIME_s;
                case WBIMMicroscopeMode.Explore
                    vol_fly_back_time_s = obj.IMAGING_EXPLORE_VOLUME_FLYBACK_TIME_s;
            end
            if obj.hSI.hFastZ.flybackTime ~= vol_fly_back_time_s
                obj.hSI.hFastZ.flybackTime = vol_fly_back_time_s;               
%                 obj.hSI.hFastZ.discardFlybackFrames = true;
            end
            
            obj.h_logger.write("DEBUG", "Finish setting stateManager", 3);
        end
    end
    %% Single layer tile acquisition control
    methods
        function obj = set_acq_mode(obj, acq_mode, im_z_r_um)            
            validateattributes(acq_mode, {'WBIMMicroscopeMode'}, {});
            obj.current_mode = acq_mode;
            obj.current_grid = obj.grid(acq_mode);        
            
            if nargin < 3 || isempty(im_z_r_um)
                obj.image_z_r_um = obj.current_grid.piezo_z_list_um;
            else
                obj.image_z_r_um = im_z_r_um;
            end
            
            if ~obj.scan_deque.isEmpty
                obj.clear_scan_queue();
            end            
            % Settings depends on the acquision mode: scan or explore
            % Might not need to save the 
            % Scan image configuration
            if obj.run_with_SI_Q
                switch acq_mode
                    case WBIMMicroscopeMode.Scan
                        obj.si_set_display_ch(obj.im_qc_channel);                        
                    case WBIMMicroscopeMode.Explore
                        obj.si_set_display_ch(obj.bg_detection_channel);
                end
                if obj.hSI.hDisplay.displayRollingAverageFactor ~= 1
                    obj.hSI.hDisplay.displayRollingAverageFactor = 1;
                end
                if ~obj.hSI.hChannels.loggingEnable
                    obj.hSI.hChannels.loggingEnable = true;
                end
                
                % This setting changes the file name in SI. 
                if obj.hSI.hScan2D.logFramesPerFile ~= inf
                    obj.hSI.hScan2D.logFramesPerFile = inf;
                end
                obj.si_set_save_dir(obj.current_scratch_folder, obj.SI_FILESTEM, 0);
                % scan angle = [angle_fast, angle_slow], which by default
                % is the x and y axis of the sample coordinate in scan
                % image
                obj.si_set_non_mroi_imaging(obj.current_grid.scan_angle_xy, ...
                    obj.current_grid.stack_size(1:2));
            else
                obj.h_logger.write("DEBUG", "hSI is not valid. Skip setting SI for stack acquisition");
            end
            obj.h_logger.write("MESSAGE", sprintf("Finish setting %s mode", ...
                acq_mode));
        end
                
        function clear_scan_queue(obj, confirmQ)
            if nargin < 2
                confirmQ = (obj.operation_mode == WBIMOperationMode.Manual);
            end
            
            if ~confirmQ
                obj.scan_deque.clear();
            else
                response = questdlg('The scan queue is not empty. Clear scan queue?', ...
                    obj.COMPONENT_NAME, 'OK', 'Cancel', 'Cancel');
                switch response
                    case 'OK'
                        obj.h_logger.write("MESSAGE", "Clear the non-empty acquisition queue");
                        obj.scan_deque.clear();
                    case 'Cancel'
                        return;
                end
            end
        end
        
        function has_valid_Q = add_tile_to_acq_queue(obj, grid_ind, checkCSQ)
            if nargin < 3
                checkCSQ = true;
            end
            % Check configuration space and exploration space
            if checkCSQ
                valid_Q = obj.check_grid_center_is_in_c_space(grid_ind);
                grid_ind = grid_ind(valid_Q);
                if ~all(valid_Q)
                    obj.h_logger.write("MESSAGE", ...
                        "Some of the candidate tiles are are not in the configuration space");
                end                
            end
            valid_exp_Q = obj.check_grid_center_is_in_e_space(grid_ind);
            grid_ind = grid_ind(valid_exp_Q);
            has_valid_Q = any(valid_exp_Q);
            if ~has_valid_Q
                obj.h_logger.write("MESSAGE", "All tiles are outside the exploraltion space");
               return 
            elseif ~all(has_valid_Q)
                obj.h_logger.write("MESSAGE", "Tiles with centers outside the exploraltion space are removed");
            end
            grid_ind = uint32(obj.current_grid.sort_grid_ind(grid_ind));
            for i = 1 : numel(grid_ind)
               if ~obj.scan_deque.contains(grid_ind(i))
                   obj.scan_deque.add(grid_ind(i));
               else
                   obj.h_logger.write("DEBUG", ...
                       sprintf('Tile %d is already in the acquisition queue.', grid_ind(i)));
               end
            end            
        end        
        
        function has_feasible_tile_Q = load_next_feasible_tile(obj)
            has_feasible_tile_Q = ~obj.scan_deque.isEmpty();
            if has_feasible_tile_Q
                obj.current_grid_ind = double(obj.scan_deque.getFirst());
                obj.current_tile_info = obj.get_tile_metadata_str();
                obj.current_grid_xy_um = obj.current_tile_info.center_xy_um;                
            end
        end
        
        function init_scanning(obj, continuous_acq_Q, reset_file_count_Q)
            arguments
                obj WBIMImaging
                continuous_acq_Q (1,1) logical = true;
                reset_file_count_Q (1,1) logical = true;
            end

            if obj.scan_deque.isEmpty()
                obj.h_logger.write("MESSAGE", 'The acquisition queue is empty. Return.');
                obj.state_finish_initialization_Q = false;
            else
                % Stop syncing 
                if obj.h_server.running_sync_Q
                    obj.h_server.running_sync_Q = false;
                end
                
                if ~isfolder(obj.current_scratch_folder)
                    mkdir(obj.current_scratch_folder);
                end
                % Set peizo park position to the first acquisition plane
                % This setting should be here as the parking position will be
                % set to the default value at the end of each acqAbort
                %             obj.set_peizo_park_z_r_um(obj.image_z_r_um(1));
                
                % Check actuator position
                if ~obj.actuator_at_dichroic_position_Q
                    obj.actuator_switch_to_dichroic();
                end
                
                obj.init_listeners();
                % Use looped acquisition
                obj.continuous_acquisition_Q = continuous_acq_Q;
                if ~obj.hSI.hStackManager.enable
                    obj.hSI.hStackManager.enable = true;
                end
                if ~obj.hSI.extTrigEnable
                    obj.hSI.extTrigEnable = true;
                end
                if ~obj.hSI.hScan2D.keepResonantScannerOn
                    obj.hSI.hScan2D.keepResonantScannerOn = true;
                end
                if reset_file_count_Q 
                    if obj.hSI.hScan2D.logFileCounter ~= 0
                        % Reset the acquisition file count
                        obj.hSI.hScan2D.logFileCounter = 0;
                    end
                end
                if obj.hSI.acqsPerLoop ~= 1e4
                    % This is an overestimate - to accomodate for re-imaging
                    obj.hSI.acqsPerLoop = 1e4;
                end
                
                if ~isempty(obj.imaging_power_profile_list)
                    if any(obj.current_mode == [WBIMMicroscopeMode.Scan, WBIMMicroscopeMode.Explore])
                        if isnan(obj.init_laser_output_power_mW)
                            obj.init_laser_power_value();
                        end                        
                        % {profile_fun, pwr_percentage}
                        pwr_profile = obj.imaging_power_profile_list{obj.current_mode};
                        if ~isnan(obj.init_laser_output_power_mW)
                            pwr_amp_factor = obj.get_laser_power_amplification_factor();
                            pwr_profile{2} = pwr_profile{2} * pwr_amp_factor;
                        end                        
                        obj.si_set_beam_adjust(pwr_profile{:});
                    else
                        obj.h_logger.write("WARNING", ...
                            "The microscope is not in either imaging nor exploration mode when setting the imaging power profile");
                    end                        
                end                
                obj.need_restart_for_abnormal_hardwares_Q = false;
                obj.mamual_pause_Q = false;
                obj.state_finish_initialization_Q = true;
                obj.current_state = WBIMMachineState.Busy;
                % Start timer here
                if ~strcmp(obj.h_state_timer.Running, 'on')
                    obj.h_state_timer.start();
                else
                    obj.h_state_timer.stop();
                    obj.h_state_timer.start();
                end
                %             obj.scan_next_tile();
            end
        end
        
        function scan_next_tile(obj)
%             obj.current_state = WBIMMachineState.Idle;
            try
                assert(~obj.scanning_in_progress_Q, ...
                    'Scanning in progress, but trying to initialize the scanning of the next tile');
                obj.normalAcqDone = false;
                check_counter = 0;
                laser_normal_Q = false;
                while ~laser_normal_Q && (check_counter < WBIMConfig.IMAGING_LASER_WAIT_COUNT)
                    laser_normal_Q = obj.laser_operation_before_imaging();
                    if ~laser_normal_Q
                        check_counter = check_counter + 1;
                        pause(WBIMConfig.IMAGING_LASER_WAIT_DURATION_s);
                        obj.h_logger.write("WARNING", "Laser abnormal during scan initialization. Wait...");
                    end
                end
                
                % Load next tile                
                if ~obj.load_next_feasible_tile() || ~laser_normal_Q
                    obj.hSI.abort();
                    if ~laser_normal_Q
                        error("Laser abnormal during scan initialization. Wait time out.");
                    end
                else
                    obj.current_state = WBIMMachineState.Busy;
                    % This information cannot be updated in loop
                    % acquisition
%                     si_info_str = sprintf('SI.WBIMTileMetadata = "%s"', ...
%                         obj.current_tile_info.fp_info_json);
%                     obj.hSI.extCustomProps = {si_info_str};
                    
                    % TODO: Set the z-dependent intensity profile here
                    
                    % Move to next tile
                    obj.h_logger.write("DEBUG", sprintf('Moving stage to grid (%d, %d) at (%.2f, %.2f) um',...
                        obj.current_tile_info.grid_sub(1), ...
                        obj.current_tile_info.grid_sub(2), obj.current_grid_xy_um));
                    % The default moveSample moves stages synchronously. No
                    % need to wait.
                    obj.move_sample_along_axis_um([1,2], obj.current_grid_xy_um, false);
                    obj.h_logger.write("DEBUG", "Finish moving stage");
                    obj.set_peizo_park_z_r_um(obj.image_z_r_um(1));
                    if obj.current_mode == WBIMMicroscopeMode.Explore
                        pause(2);
                    else
                        if numel(obj.active_channel) > 2
                            pause_t_s = 1.0;
                        else
                            pause_t_s = 2.0;
                        end
                        obj.h_logger.write("DEBUG", sprintf("Pause %.2f seconds during scan mode", ...
                            pause_t_s));
                        pause(pause_t_s);
                    end
                    obj.current_tile_info.piezo_z0_r_um = obj.image_z_r_um(1);
                    % It seems that the file count has not been updated
                    % here. The second tile in the acquisition queue has
                    % the same filepath as the first tile... 
%                     obj.current_tile_info.SI_filepath = obj.si_get_current_stack_filepath();    
                    
                    obj.scanning_in_progress_Q = true;
                    obj.reimage_current_tile_Q = false;
                    
                    obj.acq_count = obj.acq_count + 1;
                    if ~strcmp(obj.hSI.acqState, 'loop')
                        obj.hSI.startLoop();                        
                    end
                    obj.hSI.hScan2D.trigIssueSoftwareAcq();                    
                    obj.current_tile_info.t_init = datestr(now, obj.LOGGER_TIME_FORMAT);
                    pause(0.05);
                end
            catch ME
                obj.hSI.abort();
                rethrow(ME);
            end
        end
        
        function local_acq_post_processing(obj)
            try                
                tiles_info = obj.hWBIM.tile_manager.tile_info;
                processedQ = [tiles_info.file_processed_Q];
                for i = 1 : numel(processedQ)
                    tmp_tile = tiles_info(i);
                    if ~tmp_tile.file_processed_Q && ...
                            (tmp_tile.processing_state ~= WBIMProcessingState.Running)
                        if ~tmp_tile.init_SI_file_exist_Q
                            obj.h_logger.write("WARNING", sprintf("Lost raw file for tile %s", ...
                                tmp_tile(i).fp_info));
                            continue;
                        else
                            obj.h_logger.write("DEBUG", "Start post processing the acquired tile");
                            obj.h_server.submit_tile_for_processing(tmp_tile);
                            obj.h_logger.write("DEBUG", ...
                                sprintf("Submit %s tile %d to python processor", char(obj.current_mode), i));
                            if ~strcmpi(obj.hSI.acqState, 'idle')
                                % Only process one tile during acquisition
                                break;
                            else
                                switch tiles_info(i).acq_mode
                                    case WBIMMicroscopeMode.Explore
                                        continue;
                                    case WBIMMicroscopeMode.Scan
                                        pause(6);
                                end
                            end
                        end                                                
                    end
                end                
            catch ME
                obj.h_logger.write("ERROR", sprintf("Error occured when post-processing the acquired image\nError message: %s", ...
                    getReport(ME, 'extended', 'hyperlinks', 'off')));
            end
        end
        
        function pause_acquisition(obj, auto_pause_Q)
            if nargin < 2
                auto_pause_Q = false;
            end
            obj.continuous_acquisition_Q = false;
            obj.mamual_pause_Q = ~auto_pause_Q;
            obj.h_logger.write("MESSAGE", "Acquisition will stop after the current tile scan.");
        end
        
    end
    % Helper functions for continuous acquisition
    methods(Hidden)
        function laser_normal_Q = laser_operation_before_imaging(obj)
            laser_normal_Q = false;
            if ~isempty(obj.h_laser) && isvalid(obj.h_laser) && obj.h_laser.laser_connected
                obj.h_laser.stop_recording();
                if obj.h_laser.soft_key == 1
                    if obj.h_laser.target_wavelength_nm ~= obj.h_laser.wavelength_nm
                        obj.h_laser.set_wavelength_nm(obj.h_laser.target_wavelength_nm);
                        obj.h_logger.write("MESSAGE", "Setting imaging laser wavelength... wait for 10 seconds");
                        pause(10);
                    end
                    if obj.h_laser.tunable_shutter_state ~= ChameleonTunableShutter.Open
                        obj.h_laser.open_tunable_output_shutter();
                        pause(1);
                        if obj.h_laser.tunable_shutter_state == ChameleonTunableShutter.Open
                            obj.h_logger.write("MESSAGE", "Open Discovery tunable shutter");
                        else
                            obj.h_logger.write("WARNING", 'Unable to open the tunable output shutter');
                            return
                        end                        
                    end            
                    
                    laser_normal_Q = obj.h_laser.is_power_normal_Q();
                    
                    if isequal(obj.h_laser.h_timer.Running, 'off')
                        obj.h_laser.start_recording();
                    end
                else
                    obj.h_logger.write("WARNING", "The imaging laser is connected, but the soft key is off");
                    return
                end
                obj.h_laser.start_recording();
                % TODO:
                % Adjust half wave plate if applicatble
                % Record the laser state - for later adjustment
            else
                obj.h_logger.write("WARNING", "The imaging laser is not connected. Please check the shutter state and wavelength");
            end
        end
        
        function laser_operation_after_imaging(obj)
            if ~isempty(obj.h_laser) && isvalid(obj.h_laser) && obj.h_laser.laser_connected
%                 if obj.h_laser.soft_key == 1
%                     if obj.h_laser.tunable_shutter_state
%                         % Do not close the laser shutter
% %                         obj.h_laser.close_tunable_output_shutter();
% %                         obj.h_logger.write("MESSAGE", "Close Discovery tunable shutter");
%                     end
%                     % Do not stop recording the laser state? 
% %                     if isequal(obj.h_laser.h_timer.Running, 'on')
% %                         obj.h_laser.stop_recording();
% %                         obj.h_logger.write("MESSAGE", "Stop recording laser state");
% %                     end
%                 end
            end            
        end
        
        function init_listeners(obj)
            obj.delete_acq_listeners();
            % SI's implementation
            obj.hScanControlListeners = [most.ErrorHandler.addCatchingListener(obj.hSI.hUserFunctions, 'acqModeDone', @obj.end_of_acq_mode), ...
                most.ErrorHandler.addCatchingListener(obj.hSI.hUserFunctions, 'acqAbort', @obj.cleanup_acq_abort), ...
                most.ErrorHandler.addCatchingListener(obj.hSI.hUserFunctions, 'acqDone', @obj.single_tile_acq_done)];
            obj.h_logger.write("DEBUG", "Added SI tile acquisition listeners");
            obj.h_listener_laser_monitor = addlistener(obj.h_laser, 'ELaserAbnormal', ...
                @obj.imaging_laser_anomaly_detected);
            obj.h_logger.write("DEBUG", "Added laser state listener");
%             obj.h_listener_server = addlistener(obj.h_server, 'ERowShiftDetected', @obj.rowshift_detected);

            obj.h_logger.write("DEBUG", "Add pipeline control listener");
        end
        
        function single_tile_acq_done(obj, varargin)
            assert(obj.scanning_in_progress_Q, 'Scanning not in progress but acqDone is triggered. Debug');
            % This function is called on acqDone
            obj.h_logger.write("DEBUG", "Single tile acquisition done");
            obj.normalAcqDone = true;
            obj.scanning_in_progress_Q = false;
            % The following line relaies on acqDone being triggered before
            % the SI file counter increments. 
            obj.current_tile_info.SI_filepath = obj.si_get_current_stack_filepath();                
            % Update tile metadata
            if ~obj.current_tile_info.init_SI_file_exist_Q
                err_str = sprintf("Initial SI file %s does not exsit", ...
                    obj.current_tile_info.SI_filepath);
                obj.h_logger.write("ERROR", err_str);
                if obj.operation_mode == WBIMOperationMode.Manual
                    keyboard();
                else
                    obj.stop();
                end
            end
            obj.current_tile_info.normal_acq_Q = ~obj.reimage_current_tile_Q;
            obj.current_tile_info.t_done = datestr(now, obj.LOGGER_TIME_FORMAT);
            obj.current_tile_info.stage_xyz_um_done = obj.stage_xyz_um;
            obj.current_tile_info.sample_xyz_um_done = obj.sample_xyz_um;
            obj.current_tile_info.save();
            
            obj.hWBIM.tile_manager.add_imaged_tile(obj.current_tile_info.copy());
            % Tile index is removed from the queue only if it does not
            % need to be re-imaged.
            if obj.current_mode == WBIMMicroscopeMode.Scan 
                % Check row shift artefact 
                if obj.scan_image_quality_check() &&...
                        obj.reimage_count < obj.max_reimage_time
                    obj.reimage_current_tile_Q = true;
                    obj.max_reimage_time = WBIMConfig.SI_MAX_REIMAGE_TIME_ROW_SHIFT;
                    obj.continuous_acquisition_Q = false;
                    obj.need_restart_for_abnormal_hardwares_Q = true;
                end
            end            
            
            if obj.reimage_current_tile_Q 
                obj.reimage_count = obj.reimage_count + 1;                
                obj.h_logger.write("MESSAGE", "Re-image current stack");
            else
                % Only save the valid tile for acquisition analysis
                obj.acq_count = 0;
                obj.scan_deque.removeFirst();
                if obj.reimage_count >= obj.max_reimage_time
                    obj.h_logger.write("WARNING", "Reach the maximum reimage time. Move on");
                end
                obj.reimage_count = 0;
            end
            
            % Update acquisition queue in the explore mode
            if obj.current_mode == WBIMMicroscopeMode.Explore && ...
                    obj.auto_exploration
                % Analyze buffered image                
                if ~obj.explored_background_tile_Q()
                    % If the current tile is not a background tile
                    tile_neighbor_ind = obj.current_tile_info.neighbor_grid_ind;
                    tile_neighbor_ind = tile_neighbor_ind(tile_neighbor_ind > 0);
                    % Exclude tiles outside the specified bounding box
                    in_e_space_Q = obj.check_grid_center_is_in_e_space(tile_neighbor_ind);
%                     in_e_space_Q = obj.check_tile_overlaps_with_e_space(tile_neighbor_ind);
                    
                    if any(in_e_space_Q) % any([]) = false
                        tile_neighbor_ind = tile_neighbor_ind(in_e_space_Q);
                        obj.hWBIM.tile_manager.add_tile_by_grid_ind(...
                            tile_neighbor_ind, false);
                        msg = sprintf("Detect non-background tile at grid [%d, %d]. Add neighboring tiles to the candidate map", ...
                            obj.current_tile_info.grid_sub);
                        obj.h_logger.write("MESSAGE", msg);
                    end
                end
            end
            
            try
                if obj.continuous_acquisition_Q
                    if ~obj.scan_deque.isEmpty()
                        obj.scan_next_tile();
                        obj.local_acq_post_processing();
                    elseif obj.current_mode == WBIMMicroscopeMode.Explore
                        % Check if any new tile have been added to the
                        % candidate map
                        acq_ind = obj.hWBIM.tile_manager.candidate_map_to_grid_ind(...
                            false, true);
                        if ~isempty(acq_ind) && obj.add_tile_to_acq_queue(acq_ind, true)
                            obj.scan_next_tile();
                            obj.local_acq_post_processing();
                        else
                            obj.hSI.abort();
                        end
                    else
                        obj.hSI.abort();
                    end                    
                else
                    obj.hSI.abort();
                    if obj.mamual_pause_Q
                        obj.h_logger.write("MESSAGE", "Manually pause acquisition.");
                        obj.mamual_pause_Q = false;
                        if strcmp(obj.h_state_timer.Running, 'on')
                            obj.h_state_timer.stop();
                        end
                        if strcmpi(obj.h_delay_timer.timer.Running, 'on')
                            obj.h_delay_timer.timer.stop();
                        end
                    elseif obj.need_restart_for_abnormal_hardwares_Q && ~obj.scan_deque.isEmpty()
                        % Resume acquisition after shut down - to reset
                        % acquistion device state - line shift correction
                        if ~isempty(obj.h_delay_timer)
                            obj.h_delay_timer.delete();
                        end
                        % The timer is started when it is initialized. The
                        % timer only execute once.
                        obj.h_delay_timer = WBIMDelayResumeTimer(...
                            WBIMConfig.SI_ACQ_RESTART_WAITING_TIME_s, obj);
                        obj.need_restart_for_abnormal_hardwares_Q = false;
                        obj.h_logger.write("MESSAGE", sprintf("Acquisition will resume %d seconds later.", ...
                            WBIMConfig.SI_ACQ_RESTART_WAITING_TIME_s));
                    end
                end                
            catch ME
                obj.h_logger.write("ERROR", getReport(ME, 'extended', 'hyperlinks', 'off'));
                obj.hSI.abort();
            end
        end
        
        function bg_tile_Q = explored_background_tile_Q(obj)
            % Background detection does not require post-processing of the
            % acquired tile. Images are retrived from the SI buffered data
            % direclty. 
            
            buffered_data = obj.si_get_buffered_data(obj.bg_detection_channel);   
            % Compute MIP if has more than one layer
            num_z = size(buffered_data{obj.bg_detection_channel(1)}, 3);
            if num_z > 1
                buffered_data = cellfun(@(x) max(x, [], 3), ...
                    buffered_data, 'UniformOutput', false);
            end
            % Might need to change this line for selecting channel
%             im = buffered_data{1}{1};

            if obj.operation_mode == WBIMOperationMode.Manual
                % Visualize and check 
                % Does not work for 4 channels at the moment
                im = WBIMProcessExploreData.merge_channels(buffered_data);
                % Neet to transpose to match the image displayed in the
                % Channel window of SI
                if ismatrix(im)
                    im = im.';
                elseif ndims(im) == 3
                    im = permute(im, [2,1,3]);
                end                
                if isempty(obj.fig_hdl) || ~isvalid(obj.fig_hdl)
                    obj.fig_hdl = figure;
                    obj.ax_hdl = axes(obj.fig_hdl);
                else
                    obj.fig_hdl.Visible = 'on';
                end
                if isempty(obj.im_hdl) || ~isvalid(obj.im_hdl)
                    obj.im_hdl = imagesc(obj.ax_hdl, im);
                else
                    obj.im_hdl.CData = im;
                end
                obj.ax_hdl.DataAspectRatio = [1,1,1];
                response = questdlg('Is this a background tile?', ...
                    obj.COMPONENT_NAME, 'Yes', 'No', 'No');
                switch response
                    case 'Yes'
                        bg_tile_Q = true;
                    case 'No'
                        bg_tile_Q = false;
                end
                obj.fig_hdl.Visible = 'off';
            else
                mip = WBIMExplorationROIClassifier.parse_buffered_image_stack(...
                    buffered_data);
                bg_tile_Q = WBIMExplorationROIClassifier.classify_tiles_by_rules(...
                    mip);
            end
        end
        
        function restartQ = scan_image_quality_check(obj)
            if obj.auto_scan_qc_sync_Q
                buffered_data = obj.si_get_buffered_data(obj.im_qc_channel);
                rs_info = WBIMAcqQC.detect_row_shift_in_scan_stack(buffered_data, ...
                    'max_possible_rs', 30);
                % To be improved
                restartQ = any(abs(rs_info.est_shift_in_range) > 1);
                if restartQ
                    obj.h_logger.write("MESSAGE", sprintf("Detected row shift in range: %s", mat2str(rs_info.est_shift_in_range)));
                elseif any(isfinite(rs_info.est_shift))
                    obj.h_logger.write("MESSAGE", sprintf("Detected row shift out of range: %s", mat2str(rs_info.est_shift)));
                end                
            else
                restartQ = false;
            end
        end
        
        function cleanup_acq_abort(obj, varargin)
            obj.h_logger.write("DEBUG", "acqAbort starts");
            % When acqAbort is fired as part of the end of grab routine, it
            % is preceded by an acqDone event and followed by acqModeDone
            % event. In an intentional abort, acqAbort is the only event
            % fired.
            obj.scanning_in_progress_Q = false;
            % Some of these properties reset themselves when doing
            % looped mode because they are cached and restored
            % before and after the acquisition. This not an issue with
            % grab mode because each is a separate acquisition.
            %             obj.hSI.hRoiManager.mroiEnable = false;
            %             obj.hSI.abort();
            %           The following setting cannot be change during scanimage active state
            obj.hSI.hScan2D.keepResonantScannerOn = false;      
            % Setting peizo park position here has no effect as SI will
            % resotre the acquisition state after this function is run. 
            obj.set_peizo_park_z_r_um(obj.PIEZO_REFERENCE_POS_um);
            
%             if obj.scan_deque.isEmpty()
%                 notify(obj, 'EAcqDone');
%             end            
            % Delete the listeners to avoid accidental triggering
            obj.laser_operation_after_imaging();
            obj.delete_acq_listeners();            
%             if ~isempty(obj.h_laser) && isvalid(obj.h_laser)
%                 obj.h_laser.stop_recording();
%             end
            obj.current_state = WBIMMachineState.Idle;
            obj.local_acq_post_processing(); % Previously cannot read the tiff file
            obj.h_state_timer.UserData.SI_termianted_time = tic;
            obj.h_logger.write("DEBUG", "Terminated");
        end

        function end_of_acq_mode(obj, varargin)
            % Never triggered, as the listener is deleted before the
            % acqModeDone is triggered (after shuting down SI), and the
            % acquisition is terminated by abort. 
            obj.h_logger.write("DEBUG", "Triggered by SI acqModeDone");
        end
        
        function delete_acq_listeners(obj)
            if ~isempty(obj.hScanControlListeners)
                most.idioms.safeDeleteObj(obj.hScanControlListeners);
            end
            if ~isempty(obj.h_listener_laser_monitor)
                delete(obj.h_listener_laser_monitor);
            end
            if ~isempty(obj.h_listener_server)
                delete(obj.h_listener_server);
            end
        end
    end
    %% Next layer   
    methods(Access=protected)
        function update_layer_related_info(obj)
            update_layer_related_info@WBIMControlBase(obj);
            if ~isempty(obj.h_laser) && isvalid(obj.h_laser)
                % Sync laser log file
                if ~strcmpi(obj.h_laser.log_filepath, obj.laser_log_filepath)
                    try
                        % Copy the laser log file to HDD
                        archieve_log_fp = strrep(obj.h_laser.log_filepath, ...
                            obj.file_manager.fp_acquisition_scratch,...
                            obj.file_manager.fp_acquisition_disk);
                        obj.file_manager.copy_file_local_wsl(obj.h_laser.log_filepath, ...
                            archieve_log_fp, false);
                    catch ME
                        obj.h_logger.write("DEBUG", "Failed to move laser log file to local HDD store folder");
                    end
                    % Update log filepath
                    obj.h_laser.set_logger_fp(obj.laser_log_filepath);
                end                
            end
            if obj.run_with_SI_Q
                obj.si_set_save_dir(obj.scratch_root_folder, obj.SI_FILESTEM, 0);
            end
        end
    end
    %% Pipeline control - listener callback
    methods(Hidden)
        function init_timer(obj)
            if isempty(obj.h_state_timer) %|| ~isvalid(obj.h_state_timer)
                t_obj = timer();
                t_obj.Name = 'WBIMImaging State Monitor';
                t_obj.BusyMode = 'drop';
                t_obj.ExecutionMode = 'fixedSpacing';
                t_obj.Period = 5;
                t_obj.TimerFcn = {@obj.timer_state_checker, obj};
                % For debug
                t_obj.StartFcn = {@obj.timer_start, obj};
                t_obj.StopFcn = {@obj.timer_stop, obj}; 
                obj.h_state_timer = t_obj;
            else
                obj.h_logger.write("DEBUG", "The state timer has been initialized");
            end
        end
        
        function reach_end_of_imaging(obj, varargin)
            % This function is triggerd by the timer (h_state_timer)
            % It cannot be triggered using listener. The order
            % of multiple listeners to a single event is undefined and SI
            % has multiple listeners to cleanup the acquisition when
            % SI.abort() is fired. 
            obj.h_logger.write("MESSAGE", ...
                sprintf("Layer %d: Reach the end of imaging task", obj.current_layer));
            if obj.scan_deque.isEmpty()
                obj.local_acq_post_processing();
                % Check tile processing state here. Move on only when all
                % tiles have been processed. Restart an acquisition session
                % if some tiles need to be re-imaged
                if ~isempty(obj.hWBIM) && ~isempty(obj.hWBIM.tile_manager)
                    tile_manager = obj.hWBIM.tile_manager;
                    dispQ = (obj.hWBIM.operation_mode == WBIMOperationMode.Manual);
                    processDoneQ = tile_manager.check_tile_processing_state();
                    try
                        % Save stitched MIP TODO: debug
                        im_fp = obj.file_manager.fp_stitched_merged_mip(obj.experiment_group, ...
                            obj.experiment, obj.current_layer, obj.current_mode);
                        tile_manager.visualize_imaged_tile_merged_mip_sample_yx_um(...
                            'channel', obj.active_channel, 'im_filepath', im_fp, ...
                            'displayQ', dispQ);
                    catch ME
                        obj.h_logger.write("MESSAGE", sprintf("Fail to save the stitched mip\nError message: %s", ...
                            getReport(ME, 'extended', 'hyperlinks', 'off')));
                    end                    
                else
                    processDoneQ = false;
                    obj.h_logger.write("MESSAGE", "No tile was acquired. Terminate continuous operation");
                end
                
                if processDoneQ
                    if obj.current_mode == WBIMMicroscopeMode.Scan
                        msg = sprintf("Layer %d: Reach the end of scaning task", obj.current_layer);
                        obj.h_logger.write("MESSAGE", msg);
                        obj.h_notify.send_email("WBIM MESSAGE", msg);
                        notify(obj, 'EAcqScanDone');
                    elseif obj.current_mode == WBIMMicroscopeMode.Explore
                        msg = sprintf("Layer %d: Reach the end of exploration task", obj.current_layer);
                        obj.h_logger.write("MESSAGE", msg);
                        obj.h_notify.send_email("WBIM MESSAGE", msg);
                        notify(obj, 'EAcqExploreDone');
                    end
                end
            else
                obj.h_logger.write("MESSAGE", "Acquisition terminated. Exist unscanned tiles in the queue");                
            end
        end
        
    end
    methods(Static, Access=private)
        function timer_state_checker(t_obj, evnt, obj)
            % Checking unexpected interruption of the acquisition
            obj.h_logger.write("DEBUG", "Checking acquisition state");            
            if strcmpi(obj.hSI.acqState, 'idle')
                if isfinite(t_obj.UserData.SI_termianted_time)
                    si_settle_time_s = toc(t_obj.UserData.SI_termianted_time);
                else
                    si_settle_time_s = -1;
                end
                if obj.current_state == WBIMMachineState.Idle && ...
                    si_settle_time_s > t_obj.UserData.SI_abort_wait_time_s
                    % Check - probably redundant
                    if obj.normalAcqDone && obj.scan_deque.isEmpty() && ...
                            ~obj.scanning_in_progress_Q && ...
                            ~obj.hSI.active
                        t_obj.stop();
                        obj.reach_end_of_imaging();
                    elseif si_settle_time_s > t_obj.UserData.SI_abnormal_abort_wait_time_s && ...
                            (isempty(obj.h_delay_timer) || strcmpi(obj.h_delay_timer.timer.running, 'off'))
                        % Probably lost frame?
                        obj.h_logger.write("WARNING", "Exist unscanned tiles while the microscope is idle. Resume scanning");
                        assert(~obj.normalAcqDone, 'Reach this point while normalAcqDone = true');
                        if ~obj.normalAcqDone
                            % Delete current SI file
                            si_fp = obj.si_get_current_stack_filepath(obj.hSI.hScan2D.logFileCounter - 1);
                            try
                                delete(si_fp);
                                obj.h_logger.write("INFO", sprintf("Delete previous SI file: %s", ...
                                    si_fp));
                            catch ME
                                obj.h_logger.write("INFO", sprintf("Fail to delete previous SI file: %s\nError message: %s", ...
                                    si_fp, getReport(ME, 'extended', 'hyperlinks', 'off')));
                            end
                        end
                        if obj.continuous_acquisition_Q
                            obj.resume_scanning();
                        end
                    end
                else
                    % Check idle time
                    idle_duration_s = toc(t_obj.UserData.last_busy_time);
                    if idle_duration_s > t_obj.UserData.maximum_idle_time_s
                        warm_msg = sprintf("Idle time exceeds %d seconds", ...
                            t_obj.UserData.maximum_idle_time_s);
                        obj.h_logger.write("WARNING", warm_msg);
                        obj.h_notify.send_email("WBIM WARNING", warm_msg);
                        t_obj.stop();
                        % Step the acquisition?
                        obj.hWBIM.continuous_operation_Q = false;
                        obj.hWBIM.operation_mode = WBIMOperationMode.Manual;
                    end                    
                end
            else
                t_obj.UserData.last_busy_time = tic;                
            end
        end
        
        function timer_start(t_obj, evnt, obj)
            obj.h_logger.write("DEBUG", "State timer starts");
            info = struct;
            info.last_busy_time = tic;
            info.maximum_idle_time_s = WBIMConfig.SI_MAX_IDLE_TIME_s;
            info.SI_termianted_time = nan;
            info.SI_abort_wait_time_s = 3;
            info.SI_abnormal_abort_wait_time_s = WBIMConfig.SI_LOST_FRAME_RESTART_WAITING_TIME_s; % For auto-restart after losing frame
            t_obj.UserData = info;
        end
        
        function timer_stop(t_obj, evnt, obj)
            obj.h_logger.write("DEBUG", "State timer stops");
            t_obj.UserData = [];
        end
    end

    %% Anomaly detection and hanelding
    methods(Hidden)
        function imaging_laser_anomaly_detected(obj, varargin)
            obj.h_logger.write("MESSAGE", "Detect imaging laser anomaly");
            if obj.scanning_in_progress_Q
                obj.reimage_current_tile_Q = true;
                obj.max_reimage_time = WBIMConfig.SI_MAX_REIMAGE_TIME_LASER_ABNORMAL;
                obj.continuous_acquisition_Q = false;
                obj.need_restart_for_abnormal_hardwares_Q = true;
            end
            % TODO: need more handeling here. Long time abnormal state
        end

        function rowshift_detected(obj, varargin)
            % Depreciated. The feedback from the server is too slow. 
%             if obj.auto_scan_qc_async_Q
%                 obj.h_logger.write("MESSAGE", "Detect lineshift artefact");
%                 if obj.current_mode == WBIMMicroscopeMode.Scan && ...
%                         obj.current_state == WBIMMachineState.Busy
%                     % There could be tens of seconds of delay in the
%                     % row shift detection. 
% %                     obj.reimage_current_tile_Q = true;
%                     obj.pause_acquisition(true);
%                 end
%             end
        end
    end
    %% Utility
    methods
        %% Power adjustment
        function set_imaging_power_by_hwp(obj, target_val)
            arguments
                obj WBIMImaging
                target_val (1,1) {double, mustBeNonnegative} = obj.target_imaging_power_mW
            end
            current_val = obj.imaging_power_mW;
            target_ratio = target_val / current_val * obj.h_hwp.fractional_power;
            if target_ratio >= 0 && target_ratio <= 1
                obj.h_hwp.set_fraction_power(target_ratio);
            else
                obj.h_logger.write("WARNING", "Target ratio is not in [0, 1]. Skip");
            end
        end        
        
        function init_laser_power_value(obj, options)
            % Set the laser power 
            arguments
                obj (1, 1) WBIMImaging
                options.stat = 'mean';
            end
            switch options.stat
                case 'current'
                    obj.init_laser_output_power_mW = obj.h_laser.power_mW;
                case 'mean'
                    [obj.init_laser_output_power_mW, ~] = obj.h_laser.get_recent_power_stat();
                    obj.init_laser_output_power_mW = round(obj.init_laser_output_power_mW);
            end            
        end
        
        function val = get_laser_power_amplification_factor(obj, options)
            arguments
                obj (1, 1) WBIMImaging
                options.stat = 'mean';
                options.min_val (1, 1) double = 0.5;
                options.max_val (1, 1) double = 2;
            end
            if (obj.init_laser_output_power_mW > 0)
                switch options.stat
                    case 'current'
                        laser_pwr = obj.h_laser.power_mW;
                    case 'mean'
                        [laser_pwr, ~] = obj.h_laser.get_recent_power_stat();
                end
                val = obj.init_laser_output_power_mW / laser_pwr;
            else
                val = 1;
            end                  
            if val < options.min_val
                obj.h_logger.write("WARNING", "Laser output power has increased by more than 50% from its initial value!");
            elseif val > options.max_val 
                obj.h_logger.write("WARNING", "Laser output power has dropped by more than 50% from its initial value!");
            end            
        end        
        
        %% SI
        function stack_cell = si_get_buffered_data(obj, channel_list)
            if nargin < 2
                channel_list = obj.active_channel;
            end
            % obj.hSI.hDisplay.rollingStripeDataBuffer is a 1 x 1 cell
            % array, inside it, is a 1 x N cell array; 
            % This function will return a 1 x M cell array, where M is the
            % number of the largest displayed channel index
            data = obj.hSI.hDisplay.rollingStripeDataBuffer;
            num_z = numel(data);
            if num_z > 0
                data = cellfun(@(x) obj.si_parse_stripe_img_data(x{1}), data, 'UniformOutput', false);
                num_ch = numel(data{1});
                stack_cell = cell(1, num_ch);
                for i = 1 : num_ch
                    % Need to transpose the image. 
                    tmp_cell = cellfun(@(x) x{i}.', data, 'UniformOutput', false);
                    stack_cell{i} = cat(3, tmp_cell{:});
                end      
                % Check inactive channel
                empty_channel_Q = cellfun(@isempty, stack_cell);
                if any(empty_channel_Q)
                    assert(all(~empty_channel_Q(channel_list)), ...
                        'Missing active channel data');
                    stack_cell = stack_cell(channel_list);
                end
            else
                obj.h_logger.write("DEBUG", "data is empty");
            end
        end           
        
        %% Grid
        function grid_in_c_space_Q = check_grid_center_is_in_c_space(obj, grid_ind)
            tiles_info = obj.current_grid.get_stack_info(grid_ind);
            sample_xy_um = tiles_info.center_xy_um;
            num_tile = size(sample_xy_um, 1);
            tile_sample_xyz_um = cat(2, sample_xy_um, ...
                repelem(obj.sample_xyz_um(3), num_tile, 1));
            grid_in_c_space_Q = obj.sample_xyz_um_is_in_c_space_Q(tile_sample_xyz_um);
        end
        
        function im_grid = construct_2D_grid(obj, grid_name, stack_xyz_size,...
                pxl_xyz_size_um, overlap_xyz_size, z0_offset_um)
            arguments
                obj
                grid_name (1,1) WBIMMicroscopeMode
                stack_xyz_size (1, 3) double
                pxl_xyz_size_um (1, 3) double
                overlap_xyz_size (1, 3) double
                z0_offset_um (1,1) double = 0
            end            
            % TODO: the stage limit should be replaced. This function is
            % for the sample coordinate
            space_xyz_size_pxl = ceil(obj.STAGE_LIMIT_um ./ pxl_xyz_size_um);
            im_grid = WBIMImagingGrid2D(grid_name, space_xyz_size_pxl, ...
                stack_xyz_size, pxl_xyz_size_um, overlap_xyz_size, ...
                obj.IMAGING_2DGRID_AXIS_ORDER, obj.IMAGING_ACQ_AXIS_ORDER);
            % Sample space size in XY (fast, slow) 
            size_um = stack_xyz_size .* pxl_xyz_size_um;
            im_grid.scan_angle_xy = obj.compute_scan_angle(size_um(1:2));
            
            im_grid.piezo_z_list_um = im_grid.z_list_um + obj.PIEZO_REFERENCE_POS_um + ...
                z0_offset_um;
            
            obj.grid(grid_name) = im_grid;
            obj.h_logger.write("DEBUG", sprintf('Finish constructing %s 2D grid', grid_name));
        end
        
        function grid_in_e_space_Q = check_grid_center_is_in_e_space(obj, grid_ind)
            arguments
                obj
                grid_ind (1, :) double
            end
            if isempty(grid_ind)
                grid_in_e_space_Q = [];
                return
            else
                tiles_info = obj.current_grid.get_stack_info(grid_ind);
                % Check if the tile center is in the exploration space
                % bounding box
                sample_xy_um = tiles_info.center_xy_um;
                num_tile = size(sample_xy_um, 1);
                tile_sample_xyz_um = cat(2, sample_xy_um, ...
                    repelem(obj.sample_xyz_um(3), num_tile, 1));
                grid_in_e_space_Q = obj.sample_xyz_um_is_in_e_space_Q(tile_sample_xyz_um);
            end
        end
        
        function grid_in_e_space_Q = check_tile_overlaps_with_e_space(obj, grid_ind)
            arguments
                obj
                grid_ind (1, :) double
            end
            if isempty(grid_ind)
                grid_in_e_space_Q = [];
                return
            else
                tiles_info = obj.current_grid.get_stack_info(grid_ind);
                num_tiles = numel(grid_ind);
                bbox_yx_mmxx_um = tiles_info.tile_mmxx_um;
                bbox_xyz_mmxx_um = cat(2, bbox_yx_mmxx_um(:, [2, 1]), ...
                    repelem(obj.stage_xyz_um(3), num_tiles, 1), ...
                    bbox_yx_mmxx_um(:, [4,3]), ...
                    repelem(obj.stage_xyz_um(3), num_tiles, 1));
                grid_in_e_space_Q = fun_check_bounding_boxes_overlap(...
                    obj.exploration_sample_xyz_mmxx_um, ...
                    bbox_xyz_mmxx_um);
            end            
        end
        
        function validQ = sample_xyz_um_is_in_e_space_Q(obj, xyz_um)
            if iscolumnvector(xyz_um)
                xyz_um = xyz_um.';
            else
                assert(size(xyz_um, 2) == 3, 'xyz_um should have 3 columns');
            end
            validQ = all(bsxfun(@ge, xyz_um, obj.exploration_sample_xyz_mmxx_um(1:3)), 2) & ...
                all(bsxfun(@le, xyz_um, obj.exploration_sample_xyz_mmxx_um(4:6)), 2);
        end
        
        %% Transformation
        function scan_angle = compute_scan_angle(obj, fov_xy_size_um)
            % XY is in the sample coordinate?
            % These two angles do not include the fill factor for the
            % resonant scanner to turn around
            im_x_h_um = fov_xy_size_um(1) / 2;
            im_y_h_um = fov_xy_size_um(2) / 2;
            
            if isvalid(obj.hSI)
                im_corner_point = [-im_x_h_um, -im_y_h_um, 0; ...
                    im_x_h_um, - im_y_h_um, 0; ...
                    im_x_h_um, im_y_h_um, 0; ...
                    -im_x_h_um, im_y_h_um, 0];
                im_cp = scanimage.mroi.coordinates.Points(obj.hSI.hCoordinateSystems.hCSSampleRelative, ...
                    im_corner_point);
                im_cp = im_cp.transform(obj.hSI.hCoordinateSystems.hCSReference);
                im_cp = im_cp.points;
                im_cp(:, 3) = [];
                scan_angle = [abs(im_cp(1,1) - im_cp(2,1)) abs(im_cp(1,2) - im_cp(4,2))];
            else
                % Estimate scan angle according to objective focal length
                scan_angle = [2 * atand(im_x_h_um / obj.OBJECTIVE_F_mm / 1e3) * obj.POST_SCANNER_MAG , ...
                    2 * atand(im_y_h_um / obj.OBJECTIVE_F_mm/ 1e3) * obj.POST_SCANNER_MAG];
            end
            assert(all(scan_angle <= obj.SCANNER_RANGE_DEG), ...
                sprintf('Exceding maximum scan angle of [%d %d] degrees', obj.SCANNER_RANGE_DEG));
        end
        
        function max_fov_um = compute_maximum_fov_size_um(obj)
            if isvalid(obj.hSI)
                % Ref: scanimage.components.tileTools.tileExpanseTool.getAbsMaxScanTileSize();
                fov_cp = obj.hSI.hScan2D.fovCornerPoints;
                fov_cp(:, 3) = 0;
                fov_cp_um = scanimage.mroi.coordinates.Points(obj.hSI.hCoordinateSystems.hCSReference, ...
                    fov_cp);
                % I don't know how does the following transformation work... Where is the
                % matrix from scan angle to distance?
                max_fov_um = fov_cp_um.transform(obj.hSI.hMotors.hCSAlignment);
            else
                max_fov_um = obj.OBJECTIVE_F_mm * tand(obj.SCANNER_RANGE_DEG /...
                    obj.POST_SCANNER_MAG / 2) * 2 * 1000;
            end
        end
        %% Metadata
        function tile_info = get_tile_metadata_str(obj, layer_idx, acq_mode, grid_ind, acq_count)
            arguments
                obj
                layer_idx (1,1) double = obj.current_layer
                acq_mode (1,1) WBIMMicroscopeMode = obj.current_mode
                grid_ind (1,1) double = obj.current_grid_ind
                acq_count (1,1) double = obj.acq_count;
            end
            grid_hdl = obj.grid(acq_mode);
            tile_info = WBIMTileMetadata(obj.experiment_group, obj.experiment, ...
                layer_idx, acq_mode, grid_ind, grid_hdl, acq_count, ...
                obj.active_channel, obj.overwrite_file_Q);   
            tile_info.layer_z_um = obj.layer_abs_z_um;
        end
        %% IO / logging
        %         function save(obj)
        %            file_name = fullfile(obj.data_root_directory, sprintf('%s_imaging_control.mat', ...
        %                obj.data_file_base_name));
        %            % Detach the hSI without deleting the one in the base?
        %            if ~isfolder(obj.data_root_directory)
        %                mkdir(obj.data_root_directory);
        %            end
        %            save(file_name, 'obj');
        %            obj.h_logger.write("MESSAGE", "Finish saving object");
        %         end
        
        % TODO: load function        
    end
    
    methods
       % Utilities 
       function open_current_scratch_folder(obj)
           if ispc
               winopen(obj.current_scratch_folder);
           end
       end
        
       function sync_server_log_file_to_disk(obj)
           % Sync server log file to disk
           try
               server_log_filepath = obj.file_manager.fp_tcp_server_log_file(...
                   obj.experiment_group, obj.experiment);
               log_target_fp = strrep(server_log_filepath, obj.file_manager.SCRATCH_ROOT_PATH, ...
                   obj.file_manager.DATA_ROOT_PATH);
               if isfile(server_log_filepath)
                   target_folder = fileparts(log_target_fp);
                   if ~isfolder(target_folder)
                       mkdir(target_folder);
                   end
                   copyfile(server_log_filepath, log_target_fp);
                   obj.h_logger.write("MESSAGE", "Copy TCP server log to the disk");
               end
           catch ME
               fprintf('Failed to sync the server_log_file to disk\n. Error: %s', ...
                   getReport(ME, 'extended', 'hyperlinks', 'off'));
           end
       end
    end
    
    
    methods(Static)
        %% SI
        function imgData = si_parse_stripe_img_data(stripeData)
            % ScanImage's implementation
            stripeChans = stripeData.channelNumbers;
            if ~isempty(stripeData.roiData)
                imageData = stripeData.roiData{1}.imageData;
            
                chansAvail = 1 : max(stripeData.channelNumbers);

                missingChans = find(ismember(chansAvail, stripeChans)==0);
                % Insert empties for missing channels.
                for i = 1 : numel(missingChans)
                    missingChan = missingChans(i);
                    imageData = {imageData{1:missingChan-1}, {[]}, imageData{missingChan:end}};
                end

                % Remove extra cell layer.
                imgData = cellfun(@(x) x{1}, imageData, 'UniformOutput', false);
            else
                imgData = [];
            end
        end
    end
end