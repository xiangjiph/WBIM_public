classdef WBIMAblation < WBIMControlBase
%% Note: 
% 1. All the positions used in this class are in stage coordinate (um)
% 2. 
    %%
    properties(Hidden, Constant)
       COMPONENT_NAME = 'WBIMAblation';
    end
    
    % Controllable ablation parameters
    properties
        % Ablation parameters
        operation_mode (1,1) WBIMOperationMode = WBIMOperationMode.Manual
        current_state (1,1) WBIMMachineState = WBIMMachineState.Idle
        
        enable_RI_correction_Q (1, 1) logical = false;
    end
    properties(Hidden)
        RI_correction_start_delay_s (1, 1) double = 5;
        RI_correction_max_count (1, 1) double = 1;
    end
    %% 
    
    %% Dependent
    properties(Dependent)
        current_mode (1,1) WBIMMicroscopeMode
    end
    
    properties(SetObservable, SetAccess=protected)
        % Microscope state
        shutter_state (1,1) WBIMDeviceStateOnOff = WBIMDeviceStateOnOff.Unknown
        pump_state (1,1) WBIMDeviceStateOnOff = WBIMDeviceStateOnOff.Unknown
        astrella_gating_state (1,1) WBIMDeviceStateOnOff = WBIMDeviceStateOnOff.Unknown
        diffuser_state (1,1) WBIMDeviceStateOnOff = WBIMDeviceStateOnOff.Unknown        
    end
    
    properties(SetAccess=protected)        
        fractional_power (1,1) double
        % Ablation parameters
        laser_output_power_W (1,1) double
        % Ablation control
        current_num_cut (1,1) double
        target_num_cut (1,1) double    
        
        h_hwp WBIMPowerHWP
        h_refractometer RICorrection
    end 
    
    properties(Hidden, SetAccess=protected)% Development
        % Device handle
        h_astrella_gate_port
        h_pump_port
        % Diffuser 
        h_diffuser_ctrl
        h_diffuser_state
        % Refractometer
%         h_refractometer RICorrection
        % Zaber buffer   
        path_buffer_data (1, :) cell
        
        cached_ablation_paths
        
        target_cut_z_abs_list (1, :) double
        abl_power_list (1, :) double
        diffuser_state_list (1, :) logical
    end
    %% Event
    events(NotifyAccess = protected)
        EAblationLayerDone;
        EAblationModeDone;
        EPumpStateChanged;
    end
    % Event Listeners
    properties(Access=private)
        h_ablation_control_listeners = [];
        h_pump_state_listener = [];
    end
    %% Others
    properties(Access=private)
       Q_finish_init = false; 
       ini_sample_xyz_um (1,3) double
    end
    %% Life cycle
    methods
        function obj = WBIMAblation(exp_group, exp_name, hSI)
            obj = obj @ WBIMControlBase(exp_group, exp_name, hSI);
            if obj.run_with_SI_Q
                obj.init_si_related_handles();
            end
        end
        
        function delete(obj)
            delete@WBIMControlBase(obj)
            obj.cleanup_ablation_control_listeners();
            delete(obj.h_hwp);
            obj.cleanup_refractometer_related_handles();
        end        
    end
    %% Main functions
    methods        
        function execute_ablation_paths(obj, abl_path_objs, confirmQ)
            if nargin < 3
                confirmQ = (obj.operation_mode == WBIMOperationMode.Manual);
            end
            doneQ = obj.setup_volume_ablation(abl_path_objs);
            if doneQ
                obj.start_ablation(confirmQ);
            end
        end        
        
        function start_ablation(obj, ask_for_confirmation_Q)
           if nargin < 2
               ask_for_confirmation_Q = true;
           end
           if isnan(obj.laser_output_power_W) || obj.laser_output_power_W == 0
              error('Set Astrella output power first.'); 
           end           
           obj.ini_sample_xyz_um = obj.sample_xyz_um;
           if ask_for_confirmation_Q
               response = questdlg('Start ablation?', obj.COMPONENT_NAME, ...
                   'OK', 'Cancel', 'Cancel');
               switch response
                   case 'OK'
                       obj.start_next_layer_ablation();
                   case 'Cancel'
                       obj.h_logger.write("MESSAGE", 'Ablation task has been canceled');
                       obj.end_of_ablation_mode(false)
                       return;
               end               
           else
               obj.start_next_layer_ablation();
           end
        end
        
        function obj = stop(obj)
            % Interrupt the planed motion
            stop @ WBIMControlBase(obj)
            obj.end_of_ablation_mode(false);
            obj.current_state = WBIMMachineState.Idle;
            obj.h_logger.write("MESSAGE", "Emergency stop during ablation");
        end
    end    
    %% Initialization 
    methods(Hidden)
        function init_si_related_handles(obj, hSI)
            if nargin < 2
                init_si_related_handles@WBIMControlBase(obj);
            else
                init_si_related_handles@WBIMControlBase(obj, hSI);
            end
            % Digital output port for controlling the Astrella laser gate
            % Astrella gate
            rs = dabs.resources.ResourceStore();
            obj.h_astrella_gate_port = rs.filterByName(obj.ASTRELLA_GATE_PORT_NAME);
            assert(~isempty(obj.h_astrella_gate_port) && isvalid(obj.h_astrella_gate_port), 'The Astrella gating port handle is invalid')
            obj.h_logger.write("MESSAGE", "Finish initializing Astrella gating port");
            % Pump
            obj.h_pump_port = rs.filterByName(obj.ABLATION_PUMP_PORT_CHANNEL_NAME);
            assert(~isempty(obj.h_pump_port) && isvalid(obj.h_pump_port), 'The pump control port handle is invalid')
            obj.h_logger.write("MESSAGE", "Finish initializing pump port");
            % Diffuser control 
            obj.h_diffuser_ctrl = rs.filterByName(obj.ABLATION_DIFFUSER_CONTROL_PORT_NAME);
            assert(~isempty(obj.h_diffuser_ctrl) && isvalid(obj.h_diffuser_ctrl), 'The diffuser control port handle is invalid')
            obj.h_logger.write("MESSAGE", "Finish initializing diffuser control port");
            % Diffuser state
            obj.h_diffuser_state = rs.filterByName(obj.ABLATION_DIFFUSER_STATE_PORT_NAME);
            assert(~isempty(obj.h_diffuser_state) && isvalid(obj.h_diffuser_state), 'The diffuser state port handle is invalid')
            obj.h_logger.write("MESSAGE", "Finish initializing diffuser state port");
            
            % Handle to ablation beam half-wave plate
            try
                obj.init_ablation_hwp();
                obj.h_logger.write("MESSAGE", "Finish initializing ablation beam half-wave plate");
            catch ME
               obj.h_logger.write("ERROR", "Fail to initialize the ablation beam half-wave plate"); 
            end
            
            % Zaber stage monitor
            obj.h_zaber_controller.init_timer();
            obj.init_refractometer();
        end
        
        function obj = init_ablation_hwp(obj)
            obj.h_hwp = WBIMPowerHWP(WBIMConfig.ABLATION_HWP_STAGE_SERIAL_NUMBER, ...
                WBIMConfig.ABLATION_HWP_DELTA_THETA_deg, ...
                WBIMConfig.ABLATION_HWP_DEFAULT_ENERGY_FRACTION);
        end
        
        function init_refractometer(obj, opts)
            arguments
                obj (1, 1) WBIMAblation
                opts.init_listener_Q (1, 1) logical = true;
                opts.target_n (1, 1) double {mustBePositive} = WBIMConfig.RI_SET_POINT
                opts.measurement_period_s (1, 1) {mustBePositive} = WBIMConfig.RI_WAIT_TIME_s
                opts.correction_period_s (1, 1) {mustBePositive} = WBIMConfig.RI_CORRECTION_PERIOD_s
                opts.valve_v_5s_mL (1, 1) {mustBePositive} = WBIMConfig.RI_VALVE_V_5s_mL;
                opts.total_volume_mL (1, 1) {mustBePositive} = WBIMConfig.RI_TOTAL_VOLUME_mL;
            end
            obj.cleanup_refractometer_related_handles()
            obj.h_refractometer = RICorrection('target_n', opts.target_n, ...
                'measurement_period_s', opts.measurement_period_s, ...
                'correction_period_s', opts.correction_period_s, ...
                'valve_v_5s_mL', opts.valve_v_5s_mL, ...
                'total_volume_mL', opts.total_volume_mL, ...
                'logger_handle', obj.h_logger);
            if opts.init_listener_Q
                obj.setup_pump_state_listener();
            end
        end
        
        function setup_pump_state_listener(obj)
            obj.cleanup_pump_state_listener();
            if isempty(obj.h_pump_state_listener)
                obj.h_pump_state_listener = addlistener(...
                    obj, 'EPumpStateChanged', @(~, ~)obj.callback_pump_state_changed);
                obj.h_logger.write("DEBUG", 'Finish setting up pump state listener');
            end
        end
        
    end 
    %% Cleanup
    methods(Hidden)
        function cleanup_refractometer_related_handles(obj)
            if ~isempty(obj.h_refractometer) && isvalid(obj.h_refractometer)
                delete(obj.h_refractometer);                            
            end
            obj.cleanup_pump_state_listener();    
        end
        
        function cleanup_pump_state_listener(obj)
            if ~isempty(obj.h_pump_state_listener)
                delete(obj.h_pump_state_listener);
                obj.h_pump_state_listener = [];
            end
        end
    end
    %% Sets and Gets class properties
    methods                              
        function val = get.current_mode(obj)
            if obj.actuator_at_dichroic_position_Q()
                val = WBIMMicroscopeMode.Imaging;
            elseif obj.actuator_at_mirror_position_Q()
                val = WBIMMicroscopeMode.Ablation;
            else
                val = WBIMMicroscopeMode.Unknown;
            end
        end
        
        function val = get.astrella_gating_state(obj)
            if obj.run_with_SI_Q
                val = WBIMDeviceStateOnOff(obj.h_astrella_gate_port.queryValue());
            else
                val = WBIMDeviceStateOnOff.Unknown;
            end
        end
        
        function set_laser_output_power_using_orthogonal_measurement_mW(obj, oth_mW)
            obj.laser_output_power_W = obj.est_astrella_output_from_orthogonal_path_W(oth_mW/1e3);
        end
        
    end
    %% Ablation laser control
    methods
        function set_ablation_gating_state(obj, openQ, forceOpenQ)
            if nargin < 3
                forceOpenQ = false;
            end
            % Set the gating signal to the Astrella laser
            % openQ: logical scalar; True for setting the gating signal to be high;
            validateattributes(openQ, {'logical'}, {'scalar'});
            % To protect the image acquisition system. The ablation beam can
            % pass the dichroic mirror used for reflecting the excited light
            if openQ && ~forceOpenQ
                assert(obj.actuator_at_mirror_position_Q(),...
                    'The ablation mirror must be in place before turning on the laser output');
            end
            obj.h_astrella_gate_port.setValue(openQ);
        end
        
        function hDigitalTask = setup_ablation_pulses_task(obj, num_pulses, channel_name_cell)
            if nargin < 3
                channel_name_cell = {obj.h_astrella_gate_port.channelName};
            end
            validateattributes(channel_name_cell, {'cell'}, {});
            
            validateattributes(num_pulses, {'numeric'}, {'scalar', 'nonnegative'});
            assert(obj.actuator_at_mirror_position_Q(), 'The ablation mirror must be in place');
            hDigitalTask = obj.si_get_digital_task(channel_name_cell);
            hDigitalTask.sampleMode = 'finite';
%             hDigitalTask.doneCallback = @obj.si_cleanup_pulse_task;
            % Use default sampling rate
%             sample_feq_Hz = 1e6;
%             hDigitalTask.sampleRate = sample_feq_Hz;
            num_high = num_pulses * hDigitalTask.sampleRate / obj.ABLATION_REPETITION_RATE_Hz;
            buffer = cat(1, ones(num_high, 1), 0);
            hDigitalTask.writeOutputBuffer(buffer);
            % The following property determines the number of points is
            % sampled from the buffered waveform. Let it to be the size of
            % the buffer, so the buffered waveform is output once. 
            hDigitalTask.samplesPerTrigger = numel(buffer);
            obj.h_logger.write("MESSAGE", "Finish writing Astrella gate signal buffer");
            %
%             hDigitalTask.start();
%             obj.h_logger.write("MESSAGE", "Start Astrella gate signal");
        end
        
        function enable_laser_output(obj)
            % Enable Astrella output, but do not open the shutter
            % Make sure the shutter is closed
            obj.close_ablation_shutter();
            obj.set_ablation_gating_state(true, true);
            obj.h_zaber_controller.set_ablation_do_value(1);
        end
        
        function disable_laser_output(obj)
            obj.close_ablation_shutter();
            obj.set_ablation_gating_state(false);
            obj.h_zaber_controller.set_ablation_do_value(0);
        end        
    end
    
    methods(Hidden)
        function hDigitalTask = si_get_digital_task(obj, channel_name_cell)
            if nargin < 2
                channel_name_cell = {obj.h_astrella_gate_port.channelName};
            end
            hDigitalTask = dabs.vidrio.ddi.DoTask(obj.h_astrella_gate_port.hDAQ, ...
                'pulse trains');
            hDigitalTask.addChannels(channel_name_cell); 
        end
    end
    
    methods(Static, Hidden)
        function hDigitalTask = si_cleanup_pulse_task(hDigitalTask)
            fprintf('Reach the end of digital output task\n');
            if most.idioms.isValidObj(hDigitalTask)
                hDigitalTask.stop;
                hDigitalTask.abort;
                hDigitalTask.delete;
                clear hDigitalTask
                fprintf("Deleted digital output task\n");
                hDigitalTask = [];
            end
        end
    end
    %% Ablation loop control
    methods        
        function doneQ = prepare_hardware_for_ablation(obj, init_listenerQ)
            if nargin < 2
                init_listenerQ = true;
            end
            try
                if obj.run_with_SI_Q
                    if ~obj.actuator_at_mirror_position_Q
                        obj.actuator_switch_to_mirror();
                    end
                    if obj.pump_state() == WBIMDeviceStateOnOff.Off
                        obj.turn_on_pump();
                        pause(1);
                    end
                    assert(obj.actuator_at_mirror_position_Q, 'Albation mirror is not in place!');
                    obj.open_ablation_shutter();
                    if init_listenerQ 
                        obj.init_listeners();
                    end
                end
                doneQ = true;
            catch ME
                doneQ = false;
                obj.h_logger.write("ERROR", getReport(ME, 'extended', 'hyperlinks', 'off'));
            end
        end
             
        function doneQ = setup_volume_ablation(obj, abl_path_2D)
            % abl_path_2D: WBIMAblationPath2D object(s)
            
            doneQ = false;
            % Remove empty planes
            nonempty_Q = [abl_path_2D.is_not_empty_Q];
            if any(nonempty_Q)
                abl_path_2D = abl_path_2D(nonempty_Q);
                cut_z_fp_abs_um = [abl_path_2D.abl_fp_abs_z_um];
                num_plane_per_path = arrayfun(@(x) numel(x.abl_fp_abs_z_um),...
                    abl_path_2D, 'UniformOutput', true);
                num_plane = numel(cut_z_fp_abs_um);
                assert(sum(num_plane_per_path) == num_plane);
                if any(num_plane_per_path ~= 1)
                    abl_path_2D = repelem(abl_path_2D, num_plane_per_path, 1);
                end
                
                obj.path_buffer_data = cell(1, num_plane);
                obj.target_cut_z_abs_list = cut_z_fp_abs_um;
                obj.abl_power_list = [abl_path_2D.power_W];
                obj.diffuser_state_list = [abl_path_2D.diffuserQ];
                obj.target_num_cut = num_plane;
                obj.current_num_cut = 0;
                try
                    for i = 1 : num_plane
                        path_xy_um_tQ = abl_path_2D(i).trajectory;
                        num_pts = size(path_xy_um_tQ, 1);
                        % Check if the entire trajectory is in the
                        % configuration space
                        tmp_xyz_um = cat(2, path_xy_um_tQ(:, [1,2]), ...
                            repelem(cut_z_fp_abs_um(i), num_pts, 1));
                        in_c_Q = obj.stage_xyz_um_is_in_c_space_Q(tmp_xyz_um);
                        assert(all(in_c_Q), 'Not all the points in the trajectory are in the configuration space');
                        % Generate stream buffer data for later writing to
                        % the controller
                        obj.path_buffer_data{i} = {abl_path_2D(i).fast_v_um_s,...
                            abl_path_2D(i).fast_acc_m_s2, path_xy_um_tQ, 1:2, 1};
                    end
                catch ME
                    obj.turn_off_pump();
                    obj.h_logger.write("ERROR", getReport(ME, 'extended', 'hyperlinks', 'off'));
                    rethrow(ME);
                end
                obj.h_logger.write("DEBUG", 'Finish computing alation path buffers');
                % To be determined: shall we save the ablation trajectory?
%                 obj.save_ablation_parameters(abl_path_2D);
                % Todo: save the path with mask 
                obj.cached_ablation_paths = abl_path_2D;
                obj.save_ablation_parameters();
                doneQ = obj.prepare_hardware_for_ablation();
            else
                obj.h_logger.write("WARNING", 'Empty ablation path object. Terminated');
            end
        end
        
        function obj = start_next_layer_ablation(obj)
            obj.current_state = WBIMMachineState.Busy;
            incrementedQ = false;
            try
                if ~obj.actuator_at_mirror_position_Q
                    error("Ablation mirror is not in place!");
                end
                if obj.pump_state == WBIMDeviceStateOnOff.Off
                    error('Pump is off');
                end
                if obj.current_num_cut >= obj.target_num_cut
                    % This should not happen
                    obj.h_logger.write("WARNING", "Reach the maximum number of ablation planes");
                    obj.end_of_ablation_mode(false);
                else
                    obj.current_num_cut = obj.current_num_cut + 1;
                    incrementedQ = true;
                    % Move z using the piezo
                    obj.move_focal_plane_to_abs_z_um(...
                        obj.target_cut_z_abs_list(obj.current_num_cut));
                    % Write ablation path on the fly 
                    tmp_path_data = obj.path_buffer_data{obj.current_num_cut};
                    % TODO: add error handeling for this function
                    stream_buffer = obj.construct_stream_buffer(tmp_path_data{:});
%                     stream_buffer = obj.h_zaber_controller.construct_stream_buffer_from_trajectory(...
%                         tmp_path_data{:});
                    %
                    % Adjust power and set diffuser state 
                    obj.set_diffuser_state(obj.diffuser_state_list(obj.current_num_cut));
                    obj.set_power_W_after_HWP(obj.abl_power_list(obj.current_num_cut));
                    
                    obj.set_ablation_gating_state(true);
                    obj.h_zaber_controller.execute_zaber_stream_buffer(1:2, ...
                        stream_buffer);
                    obj.h_logger.write("DEBUG", "Succesfully start stream motion");
                end
            catch ME
                err_info = getReport(ME, 'extended', 'hyperlinks', 'off');
                obj.h_logger.write("ERROR", err_info);
                if incrementedQ
                    obj.current_num_cut = obj.current_num_cut - 1;
                end
                obj.end_of_ablation_mode(false);
                rethrow(ME);
            end
        end
        
        function obj = resume_ablation(obj)
            if obj.current_num_cut < obj.target_num_cut
                obj.prepare_hardware_for_ablation();
                obj.start_next_layer_ablation();
                obj.h_logger.write("MESSAGE", "Resume ablation");
            else
                obj.h_logger.write("MESSAGE", "All the buffered ablation streams have been executed");
            end            
        end
    end  
    
    methods(Access=protected)
        function obj = end_of_single_stream_motion(obj, varargin)
            obj.h_logger.write("MESSAGE", ...
                sprintf("Finish ablation plane %d/%d", obj.current_num_cut, obj.target_num_cut));
            obj.set_ablation_gating_state(false);
            if obj.current_num_cut < obj.target_num_cut
               obj.start_next_layer_ablation(); 
            else
                obj.end_of_ablation_mode(true);
                notify(obj, 'EAblationModeDone');
            end
        end
        
        function obj = end_of_ablation_mode(obj, normalQ)
            arguments
                obj WBIMAblation
                normalQ (1,1) logical = false;
            end
            % Setting hardwares
            obj.close_ablation_shutter();
            obj.cleanup_ablation_control_listeners();
            % Shut down the vDAQ port
            obj.set_ablation_gating_state(false);
            % Water pump
            obj.turn_off_pump();
            % Return to the staring position 
            obj.park_piezo();
            obj.current_state = WBIMMachineState.Idle;
            if normalQ
                obj.move_sample_along_axis_um(1:3, obj.ini_sample_xyz_um(1:3));
                obj.h_logger.write("MESSAGE", "Successfully reach the end of ablation mode");
                obj.h_notify.send_email("WBIM MESSAGE", "Reach the end of ablation mode");
            else
                mes = "Abnormal stop of the ablation process";
                obj.h_logger.write("WARNING", mes);
                obj.h_notify.send_email("WBIM WARNING", mes);                
            end
        end
        
        function init_listeners(obj)
            obj.cleanup_ablation_control_listeners();
            if isempty(obj.h_ablation_control_listeners)
                obj.h_ablation_control_listeners = addlistener(...
                    obj.h_zaber_controller, 'EMotionDone', @obj.end_of_single_stream_motion);
                obj.h_logger.write("DEBUG", 'Finish setting up ablation control listeners');
            end
        end
        
        function cleanup_ablation_control_listeners(obj)
            if ~isempty(obj.h_ablation_control_listeners)
               delete(obj.h_ablation_control_listeners); 
               obj.h_ablation_control_listeners = [];
            end
        end
    end
    
    methods
        function stream_buffer = construct_stream_buffer(obj, ...
                max_speed_um_s, max_acc_m_s2, abl_trajectory, stream_axis, buffer_id)
            try
                stream_buffer = obj.h_zaber_controller.construct_stream_buffer_from_trajectory(...
                    max_speed_um_s, max_acc_m_s2, abl_trajectory, stream_axis, buffer_id);
            catch ME
                obj.h_logger.write("WARNING", getReport(ME, 'extended', 'hyperlinks', 'off'));
                if strcmp(ME.identifier, 'MATLAB:Java:GenericException')
                    % if is out of memory error
                    try
                        obj.h_logger.write("MESSAGE", "Try to reset the controller");
                        obj.recover_from_stage_controller_error();
                        obj.h_logger.write("MESSAGE", "Try to regenerate the stream buffer");
                        stream_buffer = obj.h_zaber_controller.construct_stream_buffer_from_trajectory(...
                            max_speed_um_s, max_acc_m_s2, abl_trajectory, stream_axis, buffer_id);
                    catch ME2
                        obj.h_logger.write("WARNING", sprintf("Fail to generate the stream buffer again.\n%s", ...
                            getReport(ME2, 'extended', 'hyperlinks', 'off')));
                        rethrow(ME2);
                    end
                end
            end
        end
    end
    %% Other Hardware 
    methods
        %% Shutter control
        function val = get.shutter_state(obj)
            if obj.run_with_SI_Q
                val = WBIMDeviceStateOnOff(obj.hSI.hShutters.hShutters{obj.ABLATION_SHUTTER_IDX}.isOpen);
            else
                val = WBIMDeviceStateOnOff.Unknown;
            end                    
        end
        
        function obj = open_ablation_shutter(obj)
            if ~obj.hSI.hShutters.hShutters{obj.ABLATION_SHUTTER_IDX}.isOpen
                obj.hSI.hShutters.hShutters{obj.ABLATION_SHUTTER_IDX}.open;
                obj.h_logger.write("MESSAGE", 'Ablation beam external shutter opened');
            end
        end
        
        function obj = close_ablation_shutter(obj)
            if obj.hSI.hShutters.hShutters{obj.ABLATION_SHUTTER_IDX}.isOpen
                obj.hSI.hShutters.hShutters{obj.ABLATION_SHUTTER_IDX}.close;
                obj.h_logger.write("MESSAGE", 'Ablation beam external shutter closed'); 
            end
        end
        %% Ablation laser HWP
        function val = get.fractional_power(obj)
            if obj.run_with_SI_Q
                val = obj.h_hwp.fractional_power;
            else
                val = nan;
            end
        end
        
        function obj = set_ablation_fractional_power(obj, frac)
           assert(frac >= 0 && frac <= 1, 'Fractional power should be between 0 and 1');
           if abs(frac - obj.h_hwp.fractional_power) > 1e-6
               obj.h_hwp.set_fraction_power(frac);
               obj.h_logger.write("MESSAGE", sprintf("Set ablation beam fractional power to be %.4f%%", frac * 100));
           end
        end
        
        function val = est_astrella_output_from_orthogonal_path_W(obj, power_W)
            persistent oth_data
            if isempty(oth_data)
                oth_data = load(obj.ABLATION_OTH_DATA_FP);
            end
            if ~isempty(obj.h_hwp) && isvalid(obj.h_hwp)      
                hwp_pos = obj.h_hwp.getPosition();
                assert(hwp_pos >= oth_data.angle(1) && ...
                    hwp_pos <= oth_data.angle(end), 'Angle out of interpolation range');
                % Todo: What if the angle goes out of the range of the
                % measurement?                
                oth_power_W = oth_data.angle2power(hwp_pos);                    
                val = (power_W / oth_power_W) * obj.ABLATION_HWP_PMax_W;
            else
                val = nan;
            end
        end
        
        function val = est_ablation_power_from_orthogonal_path_W(obj, power_W)
            if ~isempty(obj.h_hwp) && isvalid(obj.h_hwp)
                est_astrella_power = obj.est_astrella_output_from_orthogonal_path_W(power_W);
                val = est_astrella_power * obj.h_hwp.fractional_power;
            else
                val = nan;
            end
        end
        
        function val = est_peak_fluence_J_per_cm2(obj, power_mw, options)
            arguments
                obj WBIMAblation
                power_mw (1,:) double = obj.laser_output_power_W
                options.diffuser_state (1,1) logical = obj.diffuser_state;
            end
            transmitted_W = power_mw * obj.ABLATION_TRANSMISSION_FRACTION * ...
                obj.h_hwp.fractional_power;            
            % Need two state variables - depends on the diffuser state
            spot_shape = WBIMSPAblation.estimate_ablation_area(options.diffuser_state);            
            val = (transmitted_W / obj.ABLATION_REPETITION_RATE_Hz)...
                / (spot_shape.deff_um_fast * spot_shape.deff_um_slow * pi / 4 / 1e8);
        end
        
        function set_fluence_J_per_cm2(obj, target_fluence, type)
            if nargin < 3
                type = 'peak';
            end
            % Get current ablation fluence
            switch type
                case 'peak'
                    current_fluence = obj.est_peak_fluence_J_per_cm2();
            end
            assert(current_fluence > 0, ...
                'Estimated current fluence is not finite. Please check if the laser output power has been set.');
            multiplier = target_fluence / current_fluence;
            target_fractional_power = obj.fractional_power * multiplier;
            assert(target_fractional_power >= 0 && target_fractional_power <=1, ...
                'Estimated target fractional power is out of range');
            obj.set_ablation_fractional_power(target_fractional_power);            
        end
        
        function set_power_W_after_HWP(obj, power_W)
            frac = power_W / obj.laser_output_power_W;
            obj.set_ablation_fractional_power(frac);            
        end
        %% Pump
        function val = get.pump_state(obj)
            if obj.run_with_SI_Q
                val = WBIMDeviceStateOnOff(obj.h_pump_port.queryValue());
            else
                val = WBIMDeviceStateOnOff.Unknown;
            end
        end
        
        function set_pump_state(obj, turnOnQ)
            validateattributes(turnOnQ, {'logical'}, {'scalar'});
            obj.h_pump_port.setValue(turnOnQ);
        end
        
        function turn_on_pump(obj)
            prior_state = obj.pump_state;
            obj.set_pump_state(true);
            obj.h_logger.write("DEBUG", "Turn on pump");
            if prior_state == WBIMDeviceStateOnOff.Off
                notify(obj, 'EPumpStateChanged');
            end
        end
        
        function turn_off_pump(obj)
            prior_state = obj.pump_state;
            obj.set_pump_state(false);
            obj.h_logger.write("DEBUG", "Turn off pump");
            if prior_state == WBIMDeviceStateOnOff.On
                notify(obj, 'EPumpStateChanged');
            end
        end
        %% Diffuser 
        function set_diffuser_state(obj, val)
            validateattributes(val, {'logical'}, {'scalar'});                
                        
            switch val
                case WBIMDeviceStateOnOff.On
                    obj.use_diffuser();
                case WBIMDeviceStateOnOff.Off
                    obj.remove_diffuser();
                case WBIMDeviceStateOnOff.Unknown
                    % Do nothing
            end
        end
        
        function val = get.diffuser_state(obj)
            if obj.run_with_SI_Q
                val = WBIMDeviceStateOnOff(obj.h_diffuser_state.queryValue());
            else
                val = WBIMDeviceStateOnOff.Unknown;
            end            
        end
                
        function use_diffuser(obj)
            if obj.diffuser_state == WBIMDeviceStateOnOff.Off
                obj.h_diffuser_ctrl.setValue(true);
                pause(1);
                assert(obj.diffuser_state == WBIMDeviceStateOnOff.On, 'Diffuser is not in place!');
                obj.h_logger.write("DEBUG", "Add diffuser");                
            end
        end
        
        function remove_diffuser(obj)
            if obj.diffuser_state == WBIMDeviceStateOnOff.On
                obj.h_diffuser_ctrl.setValue(false);
                pause(1);
                assert(obj.diffuser_state == WBIMDeviceStateOnOff.Off, 'Diffuser is still in place!');
                obj.h_logger.write("DEBUG", "Remove diffuser");
            end            
        end
        %% Refractometer
        function measure_solution_refractive_index(obj)
            if ~isempty(obj.h_refractometer) && isvalid(obj.h_refractometer)
                obj.h_refractometer.measure_RI();
            end
        end
        
        function callback_pump_state_changed(obj)
            obj.h_logger.write("DEBUG", sprintf("Reach the callback_pump_state_changed"));
            if obj.enable_RI_correction_Q && ~isempty(obj.h_refractometer) &&...
                    isvalid(obj.h_refractometer)
                if obj.pump_state == WBIMDeviceStateOnOff.On
                    doneQ = obj.h_refractometer.measure_and_correct_index_once(...
                        'start_delay_s', obj.RI_correction_start_delay_s, ...
                        'force_restart_Q', false);
                    if doneQ
                        obj.h_logger.write("DEBUG", sprintf("Start measuring and correcting for solution refractive index"));
                    end
                else
                    % Turn off the timer? 
                    doneQ = obj.h_refractometer.stop();
                    if doneQ
                        obj.h_logger.write("DEBUG", sprintf("Pump is off. Stop RICorrection"));
                    end                    
                end                
            end
        end
    end
    %% Utilities
    methods        
        function record = record_ablation_control_inputs(obj)
            record = struct;
            log_fields = {'path_buffer_data', 'target_cut_z_abs_list',...
                'abl_power_list'};
            num_field = numel(log_fields);
            for i = 1 : num_field
                tmp_fn = log_fields{i};
                record.(tmp_fn) = obj.(tmp_fn);
            end            
        end
        
        function save_ablation_parameters(obj, abl_trajectory)
            if nargin < 2
                abl_trajectory = [];
            end            
            record = obj.record_ablation_control_inputs();
            record.filepath = obj.file_manager.fp_ablation_instruction(...
                obj.experiment_group, obj.experiment, obj.current_layer);
            record.abl_trajectory = abl_trajectory;
            save_folder = fileparts(record.filepath);
            if ~isfolder(save_folder)
                mkdir(save_folder);
            end
            save(record.filepath, '-struct', 'record');
            obj.h_logger.write("DEBUG", sprintf("Finish saving ablation parameters to %s", ...
                record.filepath));
        end
        
        function abl_path_2D = construct_path_for_sample_space_bbox(obj, bbox_xyz_mmxx_um, ...
                z_step_um, ssf_fast, ssf_slow, peak_fluence_J_cm2, diffuserQ)
            arguments
                obj (1,1) WBIMAblation
                bbox_xyz_mmxx_um (1,6) double
                z_step_um (1,1) double 
                ssf_fast (1,:) double 
                ssf_slow (1,:) double 
                peak_fluence_J_cm2 (1,:) double 
                diffuserQ (1,:) logical
            end
            assert(all(bbox_xyz_mmxx_um(1:3) <= bbox_xyz_mmxx_um(4:6)));
            % Get stage space bounding box
            stage_xyz_mmxx_um = obj.sample_to_stage_bbox_xyz_um(bbox_xyz_mmxx_um);
            assert(all(obj.stage_xyz_um_is_in_c_space_Q(bbox_xyz_mmxx_um([1:3;4:6]))), ...
                sprintf("Stage bbox (%d, %d, %d, %d, %d, %d) is outside the configuration space", ...
                round(stage_xyz_mmxx_um)));
            
            z_abs_um_list = (stage_xyz_mmxx_um(3) : z_step_um : stage_xyz_mmxx_um(6))...
                + obj.piezo_z_r_um;
            abl_para = WBIMSPAblation.compute_ablation_parameters_for_multiple_planes(ssf_fast, ssf_slow, ...
                diffuserQ, peak_fluence_J_cm2);
                        
            abl_vol_3D = WBIMAblationPath3D(abl_para.fast_axis_speed_um_s, ...
                abl_para.fast_acceleration_length_um, abl_para.slow_axis_step_um, ...
                z_abs_um_list, abl_para.power_W, diffuserQ, stage_xyz_mmxx_um([1,2,4,5]));
            abl_vol_3D.construct_path_from_mask();
            abl_path_2D = abl_vol_3D.path;
        end        

        function abl_path_2D = construct_path_for_sample_space_mask(obj, ...
                mask_yx_um, mask_center_yx_um, abl_fp_z_abs_um, ...
                ssf_fast, ssf_slow, peak_fluence_J_cm2, diffuserQ, options)
            % Convert roi mask and bounding box in SAMPLE coordinate to
            % an ablation path object with coordinates in STAGE coordinate
            arguments
                obj (1,1) WBIMAblation
                mask_yx_um logical % Local mask of the ablation region in sample coordinate
                mask_center_yx_um (1, 2) double % Position of the ablation mask center in sample coordinate
                abl_fp_z_abs_um (1, :) double % Absolute z position of the ablation plane(s)
                ssf_fast (1, :) double 
                ssf_slow (1, :) double 
                peak_fluence_J_cm2 (1,:) double 
                diffuserQ (1,:) logical
                options.visQ (1,1) logical = false;
            end
            if ~isscalar(abl_fp_z_abs_um) && isscalar(peak_fluence_J_cm2)
                peak_fluence_J_cm2 = repelem(peak_fluence_J_cm2, 1, ...
                    numel(abl_fp_z_abs_um));
            end
            % Convert the sample coordinate ablation mask and mask center
            % to the stage coordinate.
            mask_xy_um = permute(obj.sample_to_stage_im_yx_um(mask_yx_um), [2,1,3]);
            size_xy_um = size(mask_xy_um, [1,2]);
            % Pad 0 in z - in our setup we assume the sample to stage
            % transformation only changes xy coordinates
            center_stage_xy_um = obj.sample_to_stage_xyz_um([mask_center_yx_um(2), ...
                mask_center_yx_um(1), 0]).';
            center_stage_xy_um = center_stage_xy_um(1:2);
            mask_xy_mm_um = round(center_stage_xy_um - size_xy_um/2 + 1);
            mask_xy_xx_um = mask_xy_mm_um + size_xy_um - 1;
            stage_bbox_xy_um = [mask_xy_mm_um, mask_xy_xx_um];
            
            abl_para = WBIMSPAblation.compute_ablation_parameters_for_multiple_planes(...
                ssf_fast, ssf_slow, diffuserQ, peak_fluence_J_cm2);
            % Construct the 2D ablation path object(s)
            abl_vol_3D = WBIMAblationPath3D(abl_para.fast_axis_speed_um_s, ...
                abl_para.fast_acceleration_length_um, abl_para.slow_axis_step_um, ...
                abl_fp_z_abs_um, abl_para.power_W, abl_para.diffuserQ, stage_bbox_xy_um, mask_xy_um);
            abl_vol_3D.construct_path_from_mask();
            abl_path_2D = abl_vol_3D.path;
            if options.visQ
                for i = 1 : numel(abl_path_2D)
                    abl_path_2D(i).visualize_mask_with_trajectory(); 
                end
            end
        end
                
        function abl_path_obj = convert_sample_roi_to_ablation_path_xy(obj,...
                abl_fp_z_abs_um, center_yx_um, mask_yx_um)
%             Deprecated
            % Convert roi mask and bounding box in SAMPLE coordinate to
            % an ablation path object with coordinates in STAGE coordinate
            % center_yx_um: 1-by-2 vector. Center of the mask in sample
            %   coordinate
            % mask_yx_um: logical matrix, local mask of the ROI in sample
            %   coordinate
            
            % Convert the sample coordinate ablation mask and mask center
            % to the stage coordinate. 
            num_mask_sec = size(mask_yx_um, 3);
            mask_xy_um = permute(obj.sample_to_stage_im_yx_um(mask_yx_um), [2,1,3]);
            size_xy_um = size(mask_xy_um, [1,2]);
            % Pad 0 in z - in our setup we assume the sample to stage
            % transformation only changes xy coordinates
            center_stage_xy_um = obj.sample_to_stage_xyz_um([center_yx_um(2), ...
                center_yx_um(1), 0]).';
            center_stage_xy_um = center_stage_xy_um(1:2);
            mask_xy_mm_um = round(center_stage_xy_um - size_xy_um/2 + 1);
            mask_xy_xx_um = mask_xy_mm_um + size_xy_um - 1;
            stage_bbox_xy_um = [mask_xy_mm_um, mask_xy_xx_um];
            % Construct the 2D ablation path object(s)
            if num_mask_sec == 1            
                abl_path_obj = WBIMAblationPath2D(obj.fast_axis_speed_um_s, ...
                    obj.fast_acceleration_length_um, obj.slow_axis_step_um, ...
                    abl_fp_z_abs_um, stage_bbox_xy_um, mask_xy_um);
            else
                assert(numel(abl_fp_z_abs_um) == num_mask_sec, ...
                    'number of focal plane z position should be the same as number of sections in the mask stack');
                abl_path_obj = cell(num_mask_sec, 1);
                for i = 1 : num_mask_sec
                    abl_path_obj{i} = WBIMAblationPath2D(obj.fast_axis_speed_um_s, ...
                        obj.fast_acceleration_length_um, obj.slow_axis_step_um, ...
                        abl_fp_z_abs_um(i), stage_bbox_xy_um, mask_xy_um(:, :, i));
                end
                abl_path_obj = cat(1, abl_path_obj{:});
            end
        end
        
        
        function abl_str = construct_ablation_path_from_AVD(obj, abl_vol_dect)
            arguments
                obj (1,1) WBIMAblation
                abl_vol_dect (:, 1) WBIMAVD
            end
            num_abl_paras = numel(abl_vol_dect.ablation_parameter);
            abl_str = cell(num_abl_paras, 1);
            for i = 1 : num_abl_paras
                tmp_abl_para = abl_vol_dect.ablation_parameter(i);
                tmp_abl_mask = abl_vol_dect.ablation_mask{i};
                tmp_abl_z_r_um = abl_vol_dect.ablation_z_r_um{i};
                tmp_abl_fp_z_um = abl_vol_dect.abs_fp_z_um_mm + tmp_abl_z_r_um;
                abl_str{i} = obj.construct_path_for_sample_space_mask(...
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
        
        % Set up 
        function exit_code = setup_flow_system_with_periodic_flow(obj, ...
                total_duration_s, pause_time_s, on_time_s)
            arguments
                obj WBIMAblation
                total_duration_s (1,1) double {mustBeNonnegative} = 1800
                pause_time_s (1,1) double {mustBeNonnegative} = 5
                on_time_s (1,1) double {mustBeNonnegative} = 10
            end
            obj.h_logger.write("MESSAGE", "Start setting up flow system with periodic flow");
            num_cycles = total_duration_s / (pause_time_s + on_time_s);
            for i = 1 : num_cycles
                obj.turn_on_pump()
                pause(on_time_s);
                obj.turn_off_pump();
                pause(pause_time_s);
            end
            obj.h_logger.write("MESSAGE", "Finish periodic flow task");
            exit_code = 0;
        end
    end
    %% Error handleing
    methods
        function recover_from_stage_controller_error(obj)
            obj.close_ablation_shutter();
            obj.turn_off_pump();
            obj.cleanup_ablation_control_listeners();
            obj.hWBIM.reset_stage_controller();
            obj.init_listeners();
            obj.turn_on_pump();
            pause(3);
            obj.open_ablation_shutter();
            obj.h_notify.send_email("WBIM MESSAGE", "Finish resetting the stage controller");
        end        
    end
end