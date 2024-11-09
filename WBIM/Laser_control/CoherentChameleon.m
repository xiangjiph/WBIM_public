classdef CoherentChameleon < handle
    %% TODO: 
    % 1. Run laser state monitor asynchronously using parallel workder 
    % 2. Write laser state directly to the SSD
    % 3. Send laser tuning information back to the parallel client - does
    % afterEach blocks execuition? 
    % 4. Imaging control writes file to SSD to notify the laser monitor -
    % for e.g. tuning - the laser state monitor needs to continuously
    % checking if an instruction file exist. 
    % 5. Add dispersion-compensation related API
    
    %%    
    properties(Access=private)
        port_id char
        baud_rate (1,1) double                
    end
    properties(Hidden)
        max_query_attempt = 5
        verbose_Q = false
    end
    
    properties(Access=protected)
       port_handle 
    end
    
    properties(SetAccess=protected)
        % State properties
       soft_key (1,1) ChameleonSoftKey = -1;
       tunable_shutter_state (1,1) ChameleonTunableShutter = -1;
       tuning_state (1,1) ChameleonTuningState = -1;
       wavelength_nm (1,1) double = 0;
       power_mW (1,1) double = 0;        
       GDD_fs2 (1,1) double = nan;
    end
    
    methods
        function obj = CoherentChameleon(port_id, baud_rate)
            if nargin < 1
                port_id = "COM7";
            end
            if nargin < 2
                baud_rate = 19200;
            end
            
            obj.port_id = port_id;
            obj.baud_rate = baud_rate;
        end
        
        function delete(obj)
            % Do not turn off the soft key
            if ~isempty(obj.port_handle)
                obj.close_tunable_output_shutter();
                obj.close_port();
            end
        end
        
        function exit_code = open_port(obj)
            exit_code = 1;
            if ~isempty(obj.port_handle)
                if obj.verbose_Q
                   fprintf('The serial port has been opened. Close and then re-open.\n'); 
                end
                obj.port_handle = [];
            end            
            obj.port_handle = serialport(obj.port_id, obj.baud_rate);
            configureTerminator(obj.port_handle, "CR/LF");
            obj.port_handle.readline();
            obj.port_handle.writeread("ECHO=0");
            obj.port_handle.writeread("Prompt=0");
            exit_code = 0;
            if obj.verbose_Q 
                fprintf('Successfully open the serial port\n');
            end
        end
        
        function obj = close_port(obj)
            obj.port_handle = [];
            if obj.verbose_Q
               fprintf('Serial port closed\n'); 
            end
        end
        
        function flush_port(obj)
            obj.port_handle.flush();
        end
    end
    %% To do
    % 1. Error message handeling
    % 2. Check states before operating the shutter
    %% Operations
    methods                
        function reply = set_soft_key(obj, n)
            assert(n == 0 || n == 1, 'key state should be either 0 or 1');
            reply = obj.port_handle.writeread(sprintf("LASER=%d", n));
        end
        
        function reply = clear_inactive_falut(obj)
            reply = obj.port_handle.writeread("FC");
        end
        
        function reply = set_tunable_output_shutter(obj, n)
            assert(n == 0 || n == 1, 'key state should be either 0 or 1');
            reply = obj.port_handle.writeread(sprintf("SHUTTER=%d", n));
        end
        
        function reply = open_tunable_output_shutter(obj)
            reply = obj.set_tunable_output_shutter(1);
        end
        
        function reply = close_tunable_output_shutter(obj)
            reply = obj.set_tunable_output_shutter(0);
        end
        
        function reply = set_wavelength_nm(obj, lambda)
            assert(lambda >= 650 && lambda <= 1300, 'The input wavelength is out of range');
            reply = obj.port_handle.writeread(sprintf("WAVELENGTH=%d", round(lambda)));
        end
        
        function reply = set_heart_beat_rate(obj, val)
           assert(val >= 0 && val <= 100, 'Heart beat rate should be between [0, 100] seconds');
           reply = obj.port_handle.writeread(sprintf("HBR=%d", round(val)));
        end
        
        function reply = set_GDD_fs2(obj, gdd_fs2)
%             assert(gdd_fs2 >= 0 && gdd_fs2 <= 17482, "GDD should be in [0, 17482] fs2");
            reply = obj.port_handle.writeread(sprintf("GDD=%d", round(gdd_fs2)));
        end
    end
    %% Query Laser State
    %     methods(Access = private)
    methods
        function info = query_info(obj, query_string, verbose_Q)
            if nargin < 3
                verbose_Q = obj.verbose_Q;
            end
            info = [];
            num_attempt = 0;
            while isempty(info) && num_attempt <= obj.max_query_attempt
                num_attempt = num_attempt + 1;
                try
                    obj.flush_port();
                    info = obj.port_handle.writeread(query_string);
                    if strcmp(info, 'COMMAND NOT DEFINED')
                        info = [];
                    end
                catch ME
                    info = [];
                    warning('Fail to query %s after %d attempts', query_string, num_attempt);
                end
            end
            if isempty(info)
               warning('Fail to query %s after %d attempts', query_string, num_attempt);
            end
            if verbose_Q
                fprintf("%s %s\n", query_string, info);
            end
        end
    end    
    %% Unit operations
    methods
        function info = get.soft_key(obj)
            info = obj.query_info("?LASER");
            try
                info = str2double(info);
                if ismember(info, [0, 1])
                    info = ChameleonSoftKey(info);
                else
                    info = ChameleonSoftKey(-1);
                end
            catch ME
                fprintf('Current info:\n');
                disp(info);
                warning('Fail to get the soft key state. Error message: %s', ...
                    getReport(ME, 'extended'));
                info = ChameleonSoftKey(-1);
            end
        end
        
        function info = get.tuning_state(obj)
            info = obj.query_info("?TS");
            try 
                info = ChameleonTuningState(str2double(info));
            catch ME
                fprintf('Current info:\n');
                disp(info);
                warning('Fail to get the tuning laser state. Error message: %s', ...
                    getReport(ME, 'extended'));
                info = ChameleonTuningState(-1);
            end                
        end
        
        function info = get.tunable_shutter_state(obj)
            info = obj.query_info("?SVAR");
            try
                info = ChameleonTunableShutter(str2double(info));
            catch ME
                fprintf('Current info:\n');
                disp(info);
                warning('Fail to get the tunable laser shutter state. Error message: %s', ...
                    getReport(ME, 'extended'));
                info = ChameleonTunableShutter(-1);
            end
        end
        
        function info = get.wavelength_nm(obj)
            info = obj.query_info("?VW");
            try
                info = str2double(info);
            catch ME
               warning('Fail to get the wavelength of the laser. Error message: %s', ...
                   getReport(ME, 'extended'));
               info = nan;
            end
        end
        
        function info = get.power_mW(obj)
            info = obj.query_info("?PVAR");
            try
                info = str2double(info);
            catch ME
               warning('Fail to get the tunable laser power. Error message: %s', ...
                   getReport(ME, 'extended'));
               info = nan;
            end
        end
        
        function info = get_active_fault(obj)
            info_1 = obj.query_info("?F");
            info_2 = obj.query_info("?FT");
            info = {info_1, info_2};
        end
        
        function info = get.GDD_fs2(obj)
           info = obj.query_info("?GDD");
           try
              info = str2double(info); 
           catch ME
               warning('Fail to get the wavelength of the laser. Error message: %s', ...
                   getReport(ME, 'extended'));
               info = nan;
           end
        end
    end
    %% Query info
    methods
        function [soft_key, is_tuning, is_opened, wavelength, power_mw] = ...
                get_tunable_laser_state(obj, verboseQ)
            if nargin < 2
                verboseQ = obj.verbose_Q;
            end
            soft_key = obj.soft_key();
            is_tuning = obj.tuning_state();
            is_opened = obj.tunable_shutter_state();
            wavelength = obj.wavelength_nm();
            power_mw = obj.power_mW();
            if verboseQ
                fprintf('Soft key: %d\tTuning state: %d\tTunable shutter state: %d\tWavelength: %d nm\tPower: %d mW\n', ...
                    soft_key, is_tuning, is_opened, wavelength, power_mw);
            end
        end
    end
end