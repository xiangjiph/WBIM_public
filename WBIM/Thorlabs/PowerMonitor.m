classdef PowerMonitor < handle
    % API for communicating with Thorlabs PM100D
    properties(Constant)
        TYPE_ID = "0x8078"
    end
    properties
        % The serial number is written on the backside of the power meter.
        serial_number string = 'P0027382';
        resource_name string
        h_device
        channel = 1;
    end
    %% Life cycle
    methods
        function obj = PowerMonitor(serial_number)
            if nargin == 1
                obj.serial_number = serial_number;
            end
            obj.resource_name = sprintf("USB0::0x1313::%s::%s::INSTR", ...
                obj.TYPE_ID, obj.serial_number);
            obj.h_device = visa('NI', obj.resource_name);
        end

        function obj = open_port(obj)
            fopen(obj.h_device);
        end

        function obj = close_port(obj)
            fclose(obj.h_device);
        end

        function pwr = get_current_power(obj)
            pwr = str2double(query(obj.h_device, sprintf(':FETCH:POW:%d:VAL?', obj.channel)));
        end
    end
    %% Utility
    methods(Static)
        function list = list_dev_list()
            list = visadevlist();
        end
    end
end
