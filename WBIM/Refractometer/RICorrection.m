classdef RICorrection < handle
    properties
        hRefractometer InlineRefractometer
        hValve SolenoidValve
        h_timer_state timer
        target_n {mustBePositive}
        measurement_period_s {mustBePositive} %s
        correction_period_s {mustBePositive} 
        total_volume_mL {mustBePositive} %mL
        fineVolume {mustBePositive} = 0.5;
        
        logQ (1, 1) logical = true;
    end

    properties (Access = private)
        running_RI_measurement_Q logical
        state {mustBeInteger, mustBeNonnegative}
        time_last_measurement {mustBeNonnegative}
        
        correction_count (1, 1) double = 0;
        max_correction_count (1, 1) double = inf;
    end
    
    properties(Hidden)
        default_timer_start_delay (1, 1) double = 0;
    end
    
    properties (Transient)
        hLogger mlog.Logger
    end

    properties (Constant, Access = private)
        RI_ELEM1 = 1.333 %ddH2O
        RI_ELEM2 = 1.4793 %DMSO
    end
    
    methods
        function obj = RICorrection(opts)
            arguments
                opts.target_n (1, 1) double {mustBePositive} = 1.4286
                opts.measurement_period_s (1, 1) {mustBePositive} = 15 * 60;
                opts.correction_period_s = 60;
                opts.valve_v_5s_mL = 8.25;   
                opts.total_volume_mL = 1500;
                opts.log_fp = [];
                opts.logger_handle = [];
            end
            try
                obj.target_n = opts.target_n;
                obj.measurement_period_s = opts.measurement_period_s;
                obj.correction_period_s = opts.correction_period_s;
                obj.total_volume_mL = opts.total_volume_mL;
                
                obj.hRefractometer = InlineRefractometer();
%                 obj.hRefractometer.logQ = false;
                obj.hRefractometer.hRIC = obj;
                obj.hValve = SolenoidValve();
                obj.hValve.v_5s_mL = opts.valve_v_5s_mL;
                
                obj.init_timer();
                if ~isempty(opts.log_fp)
                    obj.init_logger(opts.log_fp);
                elseif ~isempty(opts.logger_handle)
                    obj.hLogger = opts.logger_handle;
                end
                obj.state = 3;
            catch ME
                obj.delete();
                rethrow(ME);
            end
        end
        
        
        function val = get.running_RI_measurement_Q(obj)
            if isequal(obj.hRefractometer.timerMeasure.Running, "on")
                val = true;
            else
                val = false;
            end
        end
%% Timer
        function init_timer(obj)
            arguments
                obj (1, 1) RICorrection
            end
            obj.h_timer_state = timer();
            obj.h_timer_state.ExecutionMode = 'fixedSpacing';
            obj.h_timer_state.Period = 0.5; % second
            obj.h_timer_state.StartDelay = obj.default_timer_start_delay;
            obj.h_timer_state.TimerFcn = @(~, ~)obj.stateFunction();
            obj.h_timer_state.StopFcn = @(~, ~)obj.timer_end();
        end
%%
        function startCorrection(obj)
            if isempty(obj.hValve.v_5s_mL)
                warning("Must set hValve 5000ms calibration volume first!")
            else
                wait(obj.hRefractometer.timerMeasure)
                obj.state = 3;
                start(obj.h_timer_state)
            end
        end
        
        function done_Q = measure_and_correct_index_once(obj, opt)
            arguments
               obj (1, 1) RICorrection
               opt.start_delay_s (1, 1) double {mustBeNonnegative} = obj.default_timer_start_delay;
               opt.force_restart_Q (1, 1) logical = false;
            end
            done_Q = false;
            if isequal(obj.h_timer_state.Running, 'off') || opt.force_restart_Q
                obj.h_timer_state.stop();
                obj.h_timer_state.StartDelay = opt.start_delay_s;
                obj.correction_count = 0;
                obj.max_correction_count = 1;
                obj.startCorrection();
                done_Q = true;
            end
        end
        
        function exeQ = stop(obj)
            exeQ = false;
            if ~isempty(obj.h_timer_state) && isvalid(obj.h_timer_state) && ...
                    isequal(obj.h_timer_state.Running, 'on')
                obj.h_timer_state.stop();
                exeQ = true;
            end
        end
        
        function state = get_state(obj)
            % What does this function do???
            state = obj.state;
        end
        
        function init_logger(obj, fp)
            if isempty(which('WBIMConfig.m'))
                logger_name = 'RICorrection';
            else
                logger_name = WBIMConfig.RI_LOGGER_NANME;
            end
            obj.hLogger = mlog.Logger(logger_name, fp);
            obj.hLogger.FileThreshold = WBIMConfig.LOGGER_FILE_THRESHOLD;
            obj.hLogger.CommandWindowThreshold = WBIMConfig.LOGGER_COMMAND_WINDOW_THRESHOLD;
            obj.hLogger.MessageReceivedEventThreshold = WBIMConfig.LOGGER_EVENT_THRESHOLD;
            obj.log_info("DEBUG", "Initialized Logger");
            obj.logQ = true;
        end

        function delete(obj)
            % Don't delete the logger handle!
            if ~isempty(obj.h_timer_state)
                stop(obj.h_timer_state);
            end
            if ~isempty(obj.hRefractometer)
                obj.hRefractometer.delete;
            end
            if ~isempty(obj.hValve)
                obj.hValve.delete;
            end
            if ~isempty(obj.h_timer_state)
                obj.h_timer_state.delete;
            end
        end

        function measure_RI(obj)
            obj.hRefractometer.measure_RI();
        end
        
        function add_water_mL(obj, add_vol)
            arguments
                obj (1,1) RICorrection
                add_vol (1,1) double {mustBeNonnegative}
            end
            if add_vol > 0
                obj.hValve.valve_volume(add_vol);
            end
            % Not sure if this is necessary. Estimating the volume is hard
            % given the temperature fluctuation. 
%             obj.total_volume_mL = obj.total_volume_mL + add_vol;
        end
        
    end

    methods(Access = private)
        function stateFunction(obj)
            switch obj.state
                case 0 % intermeasurement idle
                    if toc(obj.time_last_measurement) >= obj.measurement_period_s
                        obj.state = 3; % start measurement again
                    else
                        obj.state = 0; % keep waiting
                    end
                case 1 % measurement idle
                    if obj.running_RI_measurement_Q
                        obj.state = 1; % keep waiting
                    else
                        obj.log_info("DEBUG", sprintf("RI: %0.4f\tBRIX: %0.1f",...
                            obj.hRefractometer.RI, obj.hRefractometer.BRIX))
                        obj.state = 4; % start corrections
                    end
                case 2 % intercorrection idle
                    if toc(obj.time_last_measurement) >= obj.correction_period_s
                        obj.state = 3; % start measurement again
                    else
                        obj.state = 2; % keep waiting
                    end
                case 3 % measurement
                    obj.time_last_measurement = tic;
                    obj.log_info("MESSAGE", sprintf("Start measuring RI"));
                    obj.hRefractometer.measure_RI();
                    obj.state = 1; % wait for measurement to end
                case 4 % correction
                    if obj.hRefractometer.RI > obj.target_n
                        obj.log_info("MESSAGE", sprintf("Measured refractive index is higher than the target value of %.4f.", ...
                            obj.target_n));
                        obj.adjust_solution_index();
                        obj.state = 2; % start waiting
                    else
                        obj.state = 0; % idle
                        obj.log_info("MESSAGE", sprintf("Measured refractive index is lower than the target value of %.4f.", ...
                            obj.target_n));
                    end
                    obj.correction_count = obj.correction_count + 1;
                    
                    if obj.correction_count >= obj.max_correction_count
                        obj.h_timer_state.stop();
                    end
            end
        end
        
        function timer_end(obj)
            obj.log_info("DEBUG", "RICorrection timer ends. Reset values");
            obj.correction_count = 0;
            obj.max_correction_count = inf;
            obj.h_timer_state.StartDelay = obj.default_timer_start_delay;
        end
        
    end
    
    methods
        function add_vol = adjust_solution_index(obj)
            if obj.hRefractometer.RI > obj.target_n
                est_add_vol_mL = obj.estimate_add_water_volume();
                if est_add_vol_mL > 10
                    est_add_vol_mL = 10;
                    obj.log_info("WARNING", "Unreasonably large estimated volume. Cap to the maximum value");
                end
                
                if est_add_vol_mL > (4 * obj.fineVolume)
                    add_vol = est_add_vol_mL * 0.5;
                else
                    add_vol = obj.fineVolume;
                end
                obj.add_water_mL(add_vol);
                obj.log_info("MESSAGE", sprintf('Add %.2f mL of water', add_vol));
            end
        end
        
        function add_vol = estimate_add_water_volume(obj, target_n)
            arguments
                obj (1,1) RICorrection
                target_n (1,1) double {mustBePositive} = obj.target_n;
            end
            assert(obj.total_volume_mL > 0, 'total_volume_mL needs to be initialized');
            add_vol = obj.total_volume_mL * ...
                (obj.hRefractometer.RI - target_n) / (target_n - obj.RI_ELEM1);
            assert(isscalar(add_vol) && isfinite(add_vol));
        end
        
        function turn_on_camera_view(obj)
            obj.hRefractometer.cam.preview();
        end
        
        function turn_off_camera_view(obj)
            obj.hRefractometer.cam.closePreview();
        end
        %% Logging
        function log_info(obj, msg_type, msg)
            if obj.logQ
                msg_type = upper(msg_type);
                if ~isempty(obj.hLogger) && isvalid(obj.hLogger)
                    obj.hLogger.write(msg_type, msg);
                else
                    fprintf('RICorrection %s: %s\n', msg_type, msg);
                end
            end
        end
    end

end