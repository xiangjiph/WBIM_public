classdef SolenoidValve < handle
% Class for controlling and logging data from an inline refractometer

    properties
        port_DEV %(1,1) internal.Serialport
        v_5s_mL %volume for 5000ms in mL
    end
    
    properties (Access=private)
        use_Arduino_Q % True for Arduino, false for vDAQ
    end
    
    properties
        h_sidport_valve_on
        h_sidtask_LCD_on
    end
    
    methods
        function obj = SolenoidValve(dev_port)
            arguments
                dev_port = [];
            end
            [~, hostname] = system('hostname');
            hostname = hostname(1:end-1);
            if isempty(dev_port)
                if ~strcmpi(hostname, 'pia')
                    try
                        % Begin serial communications with Arduino
                        obj.port_DEV  = serialport('/dev/ttyACM0',9600);
                        obj.use_Arduino_Q = true;
                    catch ME
                        error("Arduino port not accessible or currently in use.")
                    end
                else
                    try 
                        obj.init_si_do_task();
                    catch ME
                        obj.clean_up_si_related_handles();
                        rethrow(ME);
                    end                                        
                end     
            else
                obj.port_DEV = dev_port;
                % Need to check if dev_port is a serial port object
                obj.use_Arduino_Q = true;
            end
        end
%% Initialization
        function init_si_do_task(obj)
            rs = dabs.resources.ResourceStore();
            obj.h_sidport_valve_on = rs.filterByName(WBIMConfig.RI_VALVE_PORT_NAME);
        end     

        function clean_up_si_related_handles(obj)
            InlineRefractometer.si_cleanup_pulse_task(obj.h_sidtask_LCD_on);
            obj.h_sidport_valve_on = [];
        end
        
%%
        function open_valve_for_t_ms(obj, time_ms)
            % Open for valve for "time" ms
            arguments
                obj
                time_ms {mustBeNonnegative,mustBeInteger}
            end
            
            if obj.use_Arduino_Q
                obj.port_DEV.writeline(strcat("VALVE ",string(time_ms)));
            else
                obj.h_sidtask_LCD_on = InlineRefractometer.si_get_do_pulse_task(...
                    obj.h_sidport_valve_on, time_ms/1000, true);
                obj.h_sidtask_LCD_on.start()
            end
        end

        function valve_volume(obj, volume)
            % Open valve for "volume" mL
            arguments
                obj
                volume {mustBeNonnegative}
            end
            obj.open_valve_for_t_ms(floor(volume/obj.v_5s_mL * 5000))
        end
        
        function delete(obj)
            obj.clean_up_si_related_handles();
        end
        
    end

end
