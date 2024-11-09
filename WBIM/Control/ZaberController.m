classdef ZaberController < handle
    %% Handles
    properties
        COMPONENT_NAME = 'ZaberController';        
        h_device
        h_axes
        h_all_axes
    end
    properties(Hidden, SetAccess=private)
        connections = []
    end
    %%
    properties(Constant=true)
%         num_axes = 4;
%         axis_slow_id = 1;
%         axis_fast_id = 2;
%         axis_z_id = 3;
%         axis_actuator_id = 4;
        num_axes = WBIMConfig.STAGE_NUM_AXIS;
        axis_slow_id = WBIMConfig.STAGE_AXIS_SLOW_ID;
        axis_fast_id = WBIMConfig.STAGE_AXIS_FAST_ID;
        axis_z_id = WBIMConfig.STAGE_AXIS_Z_ID;
        axis_actuator_id = WBIMConfig.STAGE_AXIS_ACUATOR_ID;
        axis_id2name = 'xyz'
        axis_name2id = containers.Map({'x', 'X', 'y', 'Y', 'z', 'Z'}, ...
            {1, 1, 2, 2, 3, 3});
    end
    %% Options
    properties
        ablation_do_port (1,1) double = 1;
        % X-MCC4 can have up to 4 stream beging active simultaneously and
        % has 100 stream buffers  
        stream_id (1,1) double = 1;
        ablation_active_region_trigger_Q (1,1) double = true;
    end
    %% Stage state monitor
    events(NotifyAccess = protected)
       EMotionDone 
    end
    properties(Hidden)
       h_timer 
       verbose_Q (1,1) logical = false;
       monitoring_axis_id = 1 : 4;
    end
    
    %% Ablation trigger
    properties(Constant=true, Hidden)
        val_do_ablation_start logical = true;
        val_do_ablation_end logical = false;
    end
    %% Stream
    properties(Access=protected)
        measurement_array
    end        
    %% Utilities
    properties
        unit_um_per_s
        unit_um
        unit_m_per_s2
    end
    %% Life cycle
    methods
        function obj = ZaberController(axes_array)
            arguments
                axes_array = [];
            end
            % axes_array can be obtained from hSI.hMotors.hMotors{1}.hAxes;
            % Initialization
            if ~isempty(axes_array)
                obj.init_handles(axes_array);
            else
               obj.init_by_serial_port(); 
            end
            obj.init_utils();
%             obj.init_timer();
        end
        
        function delete(obj)
            delete(obj.h_timer);
            if ~isempty(obj.connections)
                try
                    obj.connections.close();
                catch ME
                    fprintf('%s', getReport(ME, 'extended', 'hyperlinks', 'off'));
                end
            end
           % Need to debug the following part - delete cell array; 
           % Maybe not...
%            delete(obj.h_all_axes);
%            delete(obj.h_device);
%            delete(obj.h_axes);
        end
    end
    %% Initialization
    methods(Hidden)
        function obj = init_handles(obj, axes_array)
            obj.h_axes = axes_array;
            obj.h_device = obj.h_axes{1}.getDevice();
            obj.h_all_axes = obj.h_device.getAllAxes();
        end
        
        function obj = init_by_serial_port(obj, port_id)
            % TODO: To be tested
            arguments
               obj (1,1) ZaberController
               port_id = WBIMConfig.STAGE_SERIAL_PORT_ID;
            end
            import zaber.motion.ascii.Connection;
            obj.connections = Connection.openSerialPort(port_id);
            try
                deviceList = obj.connections.detectDevices();
                device = deviceList(1);
                num_stage_axis = device.getAxisCount;
                stage_axis = cell(num_stage_axis, 1);
                for iter_axis = 1 : num_stage_axis
                    stage_axis{iter_axis} = device.getAxis(iter_axis);
                end
                obj.init_handles(stage_axis);
            catch ME
                obj.connections.close();
                rethrow(ME);
            end            
        end
        
        
        function obj = init_utils(obj)
            import zaber.motion.Units;
            import zaber.motion.Measurement;   
            
            obj.measurement_array = javaArray('zaber.motion.Measurement', 1);
            obj.unit_um = Units.LENGTH_MICROMETRES;
            obj.unit_um_per_s = Units.VELOCITY_MICROMETRES_PER_SECOND;
            obj.unit_m_per_s2 = Units.ACCELERATION_METRES_PER_SECOND_SQUARED;
            % TODO
        end
        
        function obj = init_timer(obj)
            t_obj = timer();
            t_obj.Name = 'Zaber Stage State Monitor';
            t_obj.BusyMode = 'drop';
            t_obj.ExecutionMode = 'fixedRate';
            t_obj.Period = 1; % second
            t_obj.TimerFcn = {@obj.monitor_check_state, obj};
            
            % For debug purpose
            t_obj.StartFcn = {@obj.monitor_start, obj};
            t_obj.StopFcn = {@obj.monitor_end, obj};
            
            obj.h_timer = t_obj;
        end
        
        function init_actuator_parameters(obj)
            % To limit the maximum current and position (?) in the actuator
        end
    end
    %% Utility
    methods
        function set_max_speed_um_s(obj, axis_id, speed_um_s)
           obj.h_axes{axis_id}.getSettings().set('maxspeed', ...
               speed_um_s, obj.unit_um_per_s);
        end
        
        % How to set the acceleration? Use COM directly? 
        function set_max_acceleration_m_per_s2(obj, axis_id, a_m_per_s2)
           obj.h_axes{axis_id}.getSettings().set('accel', ...
               a_m_per_s2, obj.unit_m_per_s2);
        end
        
        function stop(obj)
           obj.h_all_axes.stop(); 
           for i = 1 : obj.num_axes
               obj.h_axes{i}.stop();
           end
           obj.end_of_planed_motion();
        end
        
        function obj = set_port_output(obj, output_type, port_id, val)
            switch output_type
                case {'digital', 'Digital', 'D', 'd'}
                    obj.h_device.getIO().setDigitalOutput(port_id, val);
                case {'analog', 'Analog', 'A', 'a'}
                    obj.h_device.getIO().setAnalogOutput(port_id, val);
                otherwise
                    error('Unrecognize output type');
            end
        end
        
        function obj = set_ablation_do_value(obj, val)
            obj.set_port_output('digital', obj.ablation_do_port, val);
        end
        
        function val = get_ablation_do_value(obj)
           val = obj.get_port_output_value('digital', obj.ablation_do_port); 
        end
        
        function val = get_port_output_value(obj, output_type, port_id)
            switch output_type
                case {'digital', 'Digital', 'D', 'd'}
                    val = obj.h_device.getIO().getDigitalOutput(port_id);
                case {'analog', 'Analog', 'A', 'a'}
                    val = obj.h_device.getIO().getAnalogOutput(port_id);
                otherwise
                    error('Unrecognize output type');
            end            
        end
        
        function val = get_port_input_value(obj, input_type, port_id)
            switch input_type
                case {'digital', 'Digital', 'D', 'd'}
                    val = obj.h_device.getIO().getDigitalInput(port_id);
                case {'analog', 'Analog', 'A', 'a'}
                    val = obj.h_device.getIO().getAnalogInput(port_id);
                otherwise
                    error('Unrecognize input type');
            end
        end
        
        function val = get_axes_pos_um(obj, axis_id)
            if isscalar(axis_id)
                val = obj.h_axes{axis_id}.getPosition(obj.unit_um);
            else
                val = cellfun(@(x) x.getPosition(obj.unit_um), obj.h_axes(axis_id));
            end
        end
        
        function move_to_absolute_pos_um(obj, axis_id, pos_um, waitQ)
            if nargin < 4
                waitQ = false;
            end
            num_axis = numel(axis_id);
            assert(num_axis == numel(pos_um));
            if isscalar(axis_id)
                obj.h_axes{axis_id}.moveAbsolute(pos_um, obj.unit_um, waitQ)
            else
                for iter_axis = 1 : num_axis
                    tmp_id = axis_id(iter_axis);
                    obj.h_axes{tmp_id}.moveAbsolute(pos_um(iter_axis), ...
                        obj.unit_um, waitQ);
                end
            end
        end
        
        function home(obj, axis_id, waitQ)
            if nargin < 2
                axis_id = 1 : 4;
            end
            if nargin < 3
                waitQ = true;
            end            
            is_z_Q = (axis_id == obj.axis_z_id);
            non_z_ind = axis_id(~is_z_Q);           
            % Always home the z stage first
            if any(is_z_Q)
                obj.h_axes{obj.axis_z_id}.home(true);
            end 
            for iter_ax = 1 : numel(non_z_ind)
                obj.h_axes{non_z_ind(iter_ax)}.home(waitQ);
            end
            obj.h_all_axes.waitUntilIdle();
            % Move the x y stages back to the mid point: 
            obj.move_to_absolute_pos_um([obj.axis_fast_id, obj.axis_slow_id], ...
                [37500, 37500], true);
        end
        
        function [homeState] = isHomed(obj, axis_id)
            if nargin < 2
                axis_id = 1 : 4;
            end
            homeState = false(size(axis_id));
            for i = 1 : numel(axis_id)
                homeState(i) = obj.h_axes{i}.isHomed;
            end            
        end
    end
    %% Stream
    methods
        function obj = execute_line_sets_immediately(obj, line_str)
            try
                % Does the following part need to be run every time the
                % function is called?
                stage_stream = obj.h_device.getStream(1);
                stage_stream.setupLive(line_str.stream_axis_id_array_2d);
                stage_stream.setMaxSpeed(line_str.fast_speed_um_s, obj.unit_um_per_s);
                
                stage_stream = obj.write_lines_to_stream(line_str, stage_stream);
                stage_stream.disable();
                obj.h_timer.start();                
            catch ME
                stage_stream.disable();
                rethrow(ME);
            end
        end
        
        function [stream_buffer] = construct_steam_buffer_from_lines(obj, line_str,...
                buffer_id)
            import zaber.motion.Measurement;
            import zaber.motion.ascii.StreamBuffer;
            if nargin < 3
                buffer_id = 1;
            end
            
            try
                % Get zaber handles
                stream_buffer = obj.h_device.getStreamBuffer(buffer_id);
                stage_stream = obj.h_device.getStream(obj.stream_id);
                
                stream_buffer.erase(); % Make sure the buffer is empty before storing to it
                stage_stream.setupStore(stream_buffer,...
                    line_str.stream_axis_id_array_2d);
                
                stage_stream.setMaxSpeed(line_str.fast_speed_um_s, ...
                    obj.unit_um_per_s);
                stage_stream.setMaxTangentialAcceleration(line_str.fast_acceleration_m_s2,...
                    obj.unit_m_per_s2);
                
                % Write streams
                stage_stream.setDigitalOutput(obj.ablation_do_port, ...
                    obj.val_do_ablation_end);
                
                stage_stream = obj.write_lines_to_stream(line_str, stage_stream);
                
                stage_stream.setDigitalOutput(obj.ablation_do_port, ...
                    obj.val_do_ablation_end);
                
                stage_stream.disable();
            catch ME
                stage_stream.disable();
                rethrow(ME)
            end
        end
        
        function [stream_buffer] = construct_stream_buffer_from_trajectory(obj, ...
                max_speed_um_s, max_acc_m_s2, abl_trajectory, stream_axis, buffer_id)
            % abl_trajectory is for xy movement. The XY stages move much
            % faster than what the z-stage can move. Zaber current does not
            % support streaming stepping motor and servo motors with
            % different maximum speed. 
            import zaber.motion.Measurement;
            import zaber.motion.ascii.StreamBuffer;
            if nargin < 6
                buffer_id = 1;
            end
            if iscolumnvector(stream_axis)
                stream_axis = stream_axis.';
            end
            try               
                % Get zaber handles
                stream_buffer = obj.h_device.getStreamBuffer(buffer_id);
                stage_stream = obj.h_device.getStream(obj.stream_id);
                
                stream_buffer.erase(); % Make sure the buffer is empty before storing to it
                stage_stream.setupStore(stream_buffer, stream_axis);
                
                stage_stream.setMaxSpeed(max_speed_um_s, ...
                    obj.unit_um_per_s);
                stage_stream.setMaxTangentialAcceleration(max_acc_m_s2,...
                    obj.unit_m_per_s2);
                
                stage_stream.setDigitalOutput(obj.ablation_do_port, ...
                    obj.val_do_ablation_end);
                % Write stream
                abl_trajectory = abl_trajectory.';
                num_pts = size(abl_trajectory, 2);
%                 zaber_pos_array = javaArray('zaber.motion.Measurement', 3);
                previous_pos_um = nan(3, 1);
                previous_do_stage = nan;
                for iter_pts = 1 : num_pts
                    pt_xyz_um = abl_trajectory(stream_axis, iter_pts);
                    pt_do_val = abl_trajectory(end, iter_pts);
                    for i_ax = stream_axis
                        if pt_xyz_um(i_ax) ~= previous_pos_um(i_ax)
                            obj.measurement_array(1) = Measurement(pt_xyz_um(i_ax), ...
                                obj.unit_um);
                            stage_stream.lineAbsoluteOn(i_ax - 1, obj.measurement_array);
                            previous_pos_um(i_ax) = pt_xyz_um(i_ax);
                        end
%                         zaber_pos_array(i_ax) = Measurement(pt_xyz_um(i_ax), ...
%                             obj.unit_um);
                    end
%                     stage_stream.lineAbsoluteOn(0:2, zaber_pos_array);
                    if pt_do_val == 1 
                        if isnan(previous_do_stage)
                            stage_stream.setDigitalOutput(obj.ablation_do_port, ...
                                obj.val_do_ablation_start);
                            previous_do_stage = 1;
                        else
                            stage_stream.toggleDigitalOutput(obj.ablation_do_port);
                        end
                    elseif pt_do_val == -1
                        if isnan(previous_do_stage)
                            stage_stream.setDigitalOutput(obj.ablation_do_port, ...
                                obj.val_do_ablation_end);
                            previous_do_stage = -1;
                        else
                            stage_stream.toggleDigitalOutput(obj.ablation_do_port);
                        end
                    end
                end                
                % For safe, set trigger output again
                stage_stream.setDigitalOutput(obj.ablation_do_port, ...
                    obj.val_do_ablation_end);                
                stage_stream.disable();
            catch ME
                stage_stream.disable();
                rethrow(ME)
            end
        end
        
        function [stream_buffer] = construct_stream_buffer_one_step(obj, ...
                step_size_um, stream_axis, buffer_id)
            import zaber.motion.Measurement;
            import zaber.motion.ascii.StreamBuffer;
            if nargin < 4
                buffer_id = 1;
            end
            validateattributes(step_size_um, {'numeric'}, {'scalar'});
            validateattributes(stream_axis, {'numeric'}, {'scalar', 'positive'});
            try
                % Get zaber handles
                stream_buffer = obj.h_device.getStreamBuffer(buffer_id);
                stage_stream = obj.h_device.getStream(obj.stream_id);
                stream_buffer.erase(); % Make sure the buffer is empty before storing to it
                stage_stream.setupStore(stream_buffer, stream_axis);
                obj.measurement_array(1) = Measurement(step_size_um, obj.unit_um);
                % There is only one axis in this stream
                stage_stream.lineRelativeOn(0, obj.measurement_array);
                stage_stream.disable();
            catch ME
                stage_stream.disable();
                rethrow(ME)
            end
        end
                
        function stage_stream = write_lines_to_stream(obj, line_str, stage_stream)
            % Deprecated
            % The first 0 is the index for the live axes and is 0-based
            slow_ax_idx = line_str.slow_axis_id_0;
            fast_ax_idx = line_str.fast_axis_id_0;
            
            import zaber.motion.Measurement;
            previous_slow_pos = -1;
            for iter_line = 1 : line_str.num_line
                % Slow axis
                tmp_slow_pos_um = line_str.slow_ep_um_seq(iter_line);
                if tmp_slow_pos_um ~= previous_slow_pos
                    obj.measurement_array(1) = Measurement(tmp_slow_pos_um, ...
                        obj.unit_um);
                    previous_slow_pos = tmp_slow_pos_um;
                end
                stage_stream.lineAbsoluteOn(slow_ax_idx, obj.measurement_array);
                % Fast axis
                if iter_line == 1
                    obj.measurement_array(1) = Measurement(line_str.fast_ep_um_seq(1, iter_line), ...
                        obj.unit_um);
                    stage_stream.lineAbsoluteOn(fast_ax_idx, obj.measurement_array);
                end
                obj.measurement_array(1) = Measurement(line_str.fast_ep_um_seq(2, iter_line), ...
                    obj.unit_um);
                stage_stream.lineAbsoluteOn(fast_ax_idx, obj.measurement_array);
                % Enter the constant speed zone
                if obj.ablation_active_region_trigger_Q
                    stage_stream.setDigitalOutput(obj.ablation_do_port, ...
                        obj.val_do_ablation_start);
                end
                obj.measurement_array(1) = Measurement(line_str.fast_ep_um_seq(3, iter_line), ...
                    obj.unit_um);
                stage_stream.lineAbsoluteOn(fast_ax_idx, obj.measurement_array);
                % Exit the constant speed zone
                if obj.ablation_active_region_trigger_Q
                    stage_stream.setDigitalOutput(obj.ablation_do_port, ...
                        obj.val_do_ablation_end);
                end
                obj.measurement_array(1) = Measurement(line_str.fast_ep_um_seq(4, iter_line), ...
                    obj.unit_um);
                stage_stream.lineAbsoluteOn(fast_ax_idx, obj.measurement_array);
            end
        end
        
        function obj = execute_zaber_stream_buffer(obj, stream_axis_id_array, ...
                stream_buffer)
            try
                stage_stream = obj.h_device.getStream(obj.stream_id);
                if all(stage_stream.getMode == 'LIVE')
                    % Got an error here for some reason after iterrupting
                    % an albation process
                    stage_stream.disable();
                end
                stage_stream.setupLive(stream_axis_id_array);
                % Check if the initial digital output is at high:                
                if isa(stream_buffer, 'zaber.motion.ascii.StreamBuffer')
                    stage_stream.call(stream_buffer);
                elseif iscell(stream_buffer)
                    % Append all the buffers
                    for iter_cell = 1 : numel(stream_buffer)
                        stage_stream.call(stream_buffer{iter_cell});
                    end
                end
                obj.monitoring_axis_id = stream_axis_id_array;
                if isempty(obj.h_timer) || ~isvalid(obj.h_timer)
                    obj.init_timer();
                end
                obj.h_timer.start();
            catch ME
                stage_stream.disable();
                rethrow(ME);wbi
            end
        end
        
        function obj = end_of_planed_motion(obj)
            if obj.verbose_Q
                fprintf('Reach the end of the planed motion\n');
            end
            % TODO
            stage_stream = obj.h_device.getStream(obj.stream_id);
            stage_stream.disable();
            if strcmpi(obj.h_timer.Running, 'on')
                obj.h_timer.stop();
            end
            % Shut down the vDAQ port first
            %             obj.zaber_set_output('digital', obj.zaber_trigger_output_port, ...
            %                 0);
            % Determine whethere it reaches the end of a single plane or the
            % end of the ablation volume
        end
        
        function clear_stream_buffer(obj)
            for buffer_id = 1 : 20
                stream_buffer = obj.h_device.getStreamBuffer(buffer_id);
                stage_stream = obj.h_device.getStream(obj.stream_id);
                stream_buffer.erase(); % Make sure the buffer is empty before storing to it
                stage_stream.disable();
            end
        end
    end
    %% Timer functions
    methods(Static, Hidden)
        % This function has to be static due to the default input format to
        % the timer callback function
        function monitor_check_state(t_obj, evnt, obj)
            h_axes = obj.h_axes(obj.monitoring_axis_id);
            % t_data should be the handle of the axes
            if isscalar(h_axes)
                is_busy_Q = h_axes.isBusy();
            elseif iscell(h_axes)
                is_busy_Q = any(cellfun(@(x) x.isBusy(), h_axes));
            end
            if ~is_busy_Q
                if obj.verbose_Q
                    fprintf('All the stages are idle\n');
                end
                obj.end_of_planed_motion();
                t_obj.stop();
                obj.monitoring_axis_id = 1:obj.num_axes;
                notify(obj, 'EMotionDone');
            end
        end
        
        function monitor_start(t_obj, evnt, obj)
            if obj.verbose_Q
                fprintf("Stage state monitor start\n");
            end
        end
        
        function monitor_end(t_obj, evnt, obj)
            if obj.verbose_Q
                fprintf("Stage state monitor end\n");
            end
        end
    end
end