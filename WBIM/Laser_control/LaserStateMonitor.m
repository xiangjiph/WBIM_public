classdef LaserStateMonitor < CoherentChameleon
    properties
        COMPONENT_NAME = 'ImagingLaserMonitor';
    end
    properties(Dependent)
        laser_connected (1,1) logical
        time_since_last_tuning_s (1, 1) double
    end
    properties
        previous_tuning_time char
        % Control properties
        target_wavelength_nm (1,1) double
        log_filepath string
        
        notify_enabled_Q (1,1) logical = true;
        
        h_logger mlog.Logger
        % Abnormal detection
        normal_min_power_mW (1, 1) double = 0;
        normal_std_power_mW (1, 1) double = nan;
        min_post_tuning_time_s (1, 1) double = 5;
    end
    
    properties(Hidden, SetAccess=protected)
        power_history_buffer_size = 120;
        power_history_pointer (1, 1) {mustBeNonnegative} = 1;
        power_history (1, :) double = [];
        previous_tuning_t = [];
        latest_state_stat_str = []
        h_timer timer
        log_period_s (1,1) double;
        log_buffer_size (1,1) double;
    end
    
    properties(Constant, Hidden)
        LOGGER_NAME = 'Imaging Laser';
    end
    
    %%
    events(NotifyAccess=protected)
        ELaserAbnormal
    end
    %% Life cycle
    methods
        function obj = LaserStateMonitor(port_id, baud_rate, wavelength_nm, ...
                time_step_s, record_period_s, log_fp)
            obj @ CoherentChameleon(port_id, baud_rate);
            obj.verbose_Q = false;
            
            obj.target_wavelength_nm = wavelength_nm;
            obj.log_period_s = time_step_s;
            obj.log_buffer_size = record_period_s;
            obj.log_filepath = log_fp;
            obj.initialization();
        end
        
        function delete(obj)
            delete@CoherentChameleon(obj)
            if ~isempty(obj.h_timer)
                if strcmpi(obj.h_timer.Running, 'on')
                    obj.h_timer.stop();
                end
                obj.h_timer.delete();
            end
        end
    end
    %% Sets and Gets
    methods
        function val = get.laser_connected(obj)
            if isempty(obj.port_handle)
                val = false;
            else
                val = true;
            end
        end
    end
    %% Set new values
    methods
        function set_wavelength_nm(obj, val)
            set_wavelength_nm@CoherentChameleon(obj, val);
            obj.target_wavelength_nm = val;
            obj.update_state();
            obj.h_logger.write("INFO", sprintf("Set laser wavelength to %d nm", val));
        end
        
        function set_logger_fp(obj, fp)
            obj.log_filepath = fp;
            obj.h_logger.LogFile = fp;
        end
        
        function set_GDD_fs2(obj, val)
            reply = set_GDD_fs2@CoherentChameleon(obj, val);
            obj.h_logger.write("INFO", sprintf("Set laser GDD to %d fs2", val));
        end
    end
    %% Initialization
    methods
        function initialization(obj)
            obj.init_logger();
            obj.init_timer();
            obj.setup_laser();
            obj.init_power_statistics();
        end
        
        function inititialize_and_luanch(obj)
            obj.initialization();
            obj.start_recording();
        end
        
        function setup_laser(obj)
            try
                obj.connect_to_laser();
            catch ME
                obj.h_logger.write("ERROR", getReport(ME, 'extended', 'hyperlinks', 'off'));
                rethrow(ME);
            end
            try
                obj.turn_on_laser();
                obj.set_wavelength_nm(obj.target_wavelength_nm);
            catch ME
                obj.close_port();
                obj.h_logger.write("ERROR", getReport(ME, 'extended', 'hyperlinks', 'off'));
                rethrow(ME)
            end
        end
    end
    %% Operations
    methods
        function init_logger(obj)
            [folder, ~] = fileparts(obj.log_filepath);
            if ~isfolder(folder)
                mkdir(folder);
            end
            obj.h_logger = mlog.Logger(obj.LOGGER_NAME, obj.log_filepath);
            obj.h_logger.FileThreshold = mlog.Level.DEBUG;
            obj.h_logger.CommandWindowThreshold = mlog.Level.MESSAGE;
            obj.h_logger.MessageReceivedEventThreshold = mlog.Level.WARNING;
        end
        
        function connect_to_laser(obj)
            try
                obj.open_port();
            catch ME
                obj.h_logger.write("ERROR", getReport(ME, 'extended', 'hyperlinks', 'off'));
            end
        end
        
        function disconnect_from_laser(obj)
            % Only close the port. Do not turn off the laser
            % Note that the laser might automatically shut down in the
            % abscence of RS232 activity
            obj.h_timer.stop();
            obj.close_port();
        end
        
        function turn_on_laser(obj)
            if ~obj.soft_key()
                obj.set_soft_key(true);
                obj.h_logger.write("INFO", "Turning on the laser soft key");
            else
                obj.h_logger.write("INFO", "The laser is already on");
            end
            obj.update_state();
        end
        
        function turn_off_laser(obj)
            obj.h_timer.stop();
            if obj.tunable_shutter_state
                obj.close_tunable_output_shutter();
                obj.h_logger.write("INFO", "Close the tunable laser shutter");
            end
            if obj.soft_key()
                obj.set_soft_key(0);
                obj.h_logger.write("INFO", "Imaging laser off");
            end
            obj.update_state();
        end
        
        function obj = update_state(obj)
            if obj.tuning_state
                obj.previous_tuning_time = datestr(now);
                obj.previous_tuning_t = tic;
            end
        end
        
        function str = get_laser_state_string(obj)
            str = sprintf("SoftKey: %s\tTuning: %s\tTunableLaserShutter: %s\tWavelength: %d\tTunableLaserPower: %d", ...
                obj.soft_key, obj.tuning_state, obj.tunable_shutter_state,...
                obj.wavelength_nm, obj.power_mW);
        end
        
        function val = get.time_since_last_tuning_s(obj)
            if ~isempty(obj.previous_tuning_t)
               val = toc(obj.previous_tuning_t); 
            end
        end
    end
    %%
    methods
        function start_recording(obj)
            if isempty(obj.h_timer)
                obj.init_timer();
            end
            if obj.laser_connected
                if isequal(obj.h_timer.Running, 'off')
                    info_str = sprintf("Start recording imaging laser state");
                    obj.h_timer.start();
                    obj.h_logger.write("DEBUG", info_str);
                end
            else
                obj.h_logger.write("WARNING", "Not connected to the laser. Skip recording");
            end
        end
        
        function stop_recording(obj)
            if isequal(obj.h_timer.Running, 'on')
                stop(obj.h_timer);
                info_str = sprintf("Stop recording imaging laser state");
                obj.h_logger.write("DEBUG", info_str);
            end
        end
    end
    %% Timer
    methods
        function init_timer(obj)
            t_obj = timer();
            t_obj.Name = 'Imaging Laser State Monitor';
            t_obj.BusyMode = 'drop';
            t_obj.ExecutionMode = 'fixedSpacing';
            t_obj.Period = obj.log_period_s; % second
            t_obj.TimerFcn = {@obj.timer_fun, obj};
            
            % For debug purpose
            t_obj.StartFcn = {@obj.timer_start, obj};
            t_obj.StopFcn = {@obj.timer_end, obj};
            
            obj.h_timer = t_obj;
        end
    end
    
    methods(Static, Hidden)
        % This function has to be static due to the default input format to
        % the timer callback function
        function timer_fun(tobj, evnt, obj)
            try
                r_str = tobj.UserData;
                
                r_str.num_record = r_str.num_record + 1;
                obj.update_state();
                state_string = obj.get_laser_state_string();
                obj.h_logger.write("DEBUG", state_string);
                
                r_str.time(r_str.num_record) = now;
                r_str.soft_key_on_Q(r_str.num_record) = obj.soft_key;
                r_str.is_tuning_Q(r_str.num_record) = obj.tuning_state;
                r_str.shutter_opened_Q(r_str.num_record) = obj.tunable_shutter_state;
                r_str.wavelength_nm(r_str.num_record) = obj.wavelength_nm;
                r_str.power_mW(r_str.num_record) = obj.power_mW;
                if obj.is_tuning_done_Q()
                    obj.update_power_history(r_str.power_mW(r_str.num_record))
                end
                
                if r_str.num_record == r_str.buffer_size
                    power_stat = obj.compute_basic_statistics(r_str.power_mW);
                    stat_str = sprintf('Basic statistics about the power in the previous record section:\n\tMean: %.2f mW\tSTD: %.2f mW\tMin: %.2f mW\tMedian: %.2f\tMax: %.2f mW', ...
                        power_stat.mean, power_stat.std, power_stat.prctile_val(1), ...
                        power_stat.prctile_val(5), power_stat.prctile_val(end));
                    obj.h_logger.write("DEBUG", stat_str);
                    obj.latest_state_stat_str = stat_str;
                    % Re-initialize the log structure
                    r_str = obj.init_log_str(obj.log_filepath, ...
                        obj.log_buffer_size);
                end
                
                if obj.notify_enabled_Q
                    if ~obj.soft_key
                        notify(obj, 'ELaserAbnormal', WBIMEvent(ChameleonState.SoftKeyOff));
                        obj.h_logger.write("WARNING", "Laser OFF");
                    elseif obj.tuning_state
                        notify(obj, 'ELaserAbnormal', WBIMEvent(ChameleonState.Tuning));
                        obj.h_logger.write("WARNING", "Laser tuning");
                    elseif obj.is_low_power_Q()
                        notify(obj, 'ELaserAbnormal', WBIMEvent(ChameleonState.LowPower));
                        obj.h_logger.write("WARNING", "Laser power low");
                        % TODO: Check power fluctuation
                    end
                end
                tobj.UserData = r_str;
            catch ME
                rethrow(ME);
            end
            % Determine whether the laser power is abnormal
        end
        
        function timer_start(tobj, evnt, obj)
            % Initialize buffer
            r = LaserStateMonitor.init_log_str(obj.log_filepath, ...
                obj.log_buffer_size);
            [log_folder, ~] = fileparts(r.filepath);
            if ~isfolder(log_folder)
                mkdir(log_folder);
            end
            tobj.UserData = r;
        end
        
        function timer_end(tobj, evnt, obj)
            %             tobj.stop();
            % Write the remaining part in the recorded buffer
            % Close the laser port
        end
        
        function r = init_log_str(filepath, buffer_size)
            r = struct;
            r.filepath = filepath;
            r.buffer_size = buffer_size;
            r.num_record = 0;
            [r.time, r.wavelength_nm, r.power_mW] = deal(nan(r.buffer_size, 1));
            [r.shutter_opened_Q, r.is_tuning_Q, r.soft_key_on_Q] = ...
                deal(zeros(r.buffer_size, 1, 'int8'));
        end
        
        function exit_code = write_state_log(filename, record_str)
            record_str.time = datestr(record_str.time, 'yyyymmdd hh:MM:ss').';
            
            [fp, ~, ext] = fileparts(filename);
            exit_code = 1;
            if ~isfolder(fp)
                mkdir(fp);
            end
            assert(any(strcmpi(ext, {'.CSV', '.TXT'})), 'Only support writing csv file at the moment');
            f_h = fopen(filename, 'a');
            try
                if ~isfile(filename)
                    fprintf(f_h, 'Time, SoftKeyState, TuningState, ShutterState, Wavelength_nm, Power_mW\n');
                end
                
                num_record = numel(record_str.power_mW);
                for iter_t = 1 : num_record
                    fprintf(f_h, '%s, %d, %d, %d, %d, %d\n', ...
                        record_str.time(:, iter_t), record_str.soft_key_on_Q(iter_t), ...
                        record_str.is_tuning_Q(iter_t), record_str.shutter_opened_Q(iter_t), ...
                        record_str.wavelength_nm(iter_t), record_str.power_mW(iter_t));
                end
                fclose(f_h);
                exit_code = 0;
            catch ME
                fclose(f_h);
                rethrow(ME);
            end
        end
    end
    %% Statistics
    % Abnormal detection
    methods
        function init_power_statistics(obj)
            obj.power_history = nan(obj.power_history_buffer_size, 1);
            obj.power_history_pointer = 1;                        
        end
        
        function update_power_history(obj, new_pwr_mW)
            if obj.power_history_pointer <= obj.power_history_buffer_size
                obj.power_history(obj.power_history_pointer) = new_pwr_mW;
                obj.power_history_pointer = obj.power_history_pointer + 1;
            else
                obj.power_history = [obj.power_history(2:end), new_pwr_mW];
            end
        end
        
        function [pm, pstd] = get_recent_power_stat(obj)
            pm = mean(obj.power_history, 'omitnan');
            pstd = std(obj.power_history, 'omitnan');            
        end
        
        function normal_Q = is_power_normal_Q(obj)
            normal_Q = ~obj.is_low_power_Q();
            % Tuning has been done for xx seconds. 
            normal_Q = normal_Q && obj.is_tuning_done_Q();
            
            % TODO: Check fluctuation; autoregression to compute the mean
            % and variance. 
        end
        
        function doneQ = is_tuning_done_Q(obj)
            doneQ = (obj.time_since_last_tuning_s > obj.min_post_tuning_time_s);
        end        
        
        function low_Q = is_low_power_Q(obj)
            low_Q = (obj.power_mW < obj.normal_min_power_mW);
        end        
    end
    
    methods(Static)
        function stat = compute_basic_statistics(data)
            validateattributes(data, {'numeric'}, {'vector'});
            data = double(data);
            stat.num_data = numel(data);
            stat.sum = sum(data);
            stat.mean = stat.sum ./ stat.num_data;
            stat.std = sqrt(sum(data .^ 2) / stat.num_data - stat.mean.^2);
            stat.cv = stat.std ./ stat.mean;
            stat.prctile_th = [0, 1, 5, 25, 50, 75, 95, 99, 100];
            stat.prctile_val = prctile(data, stat.prctile_th);
        end
        
        function result = analyze_laser_state(r_str)
            result = struct;
            result.start_time = r_str.time(1);
            result.end_time = r_str.time(end);
            
            soft_key_u = unique(r_str.soft_key_on_Q);
            result.softkey_normal_Q = isscalar(soft_key_u);
            
            tuning_state_u = unique(r_str.is_tuning_Q);
            result.no_tuning_Q = all(tuning_state_u == 0);
            
            tunable_laser_shutter_u = unique(r_str.shutter_opened_Q);
            result.tunable_shutter_normal_Q = isscalar(tunable_laser_shutter_u);
            
            laser_wavelength_u = unique(r_str.wavelength_nm);
            result.wavelength_nm = laser_wavelength_u;
            
            result.power_stat = LaserStateMonitor.compute_basic_statistics(r_str.power_mW);
            
        end
        
    end
    %% Utilities
    % To do:
    % 1. GUI for interacting with discovery
    % 2.
end