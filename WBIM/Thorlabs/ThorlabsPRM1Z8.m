classdef ThorlabsPRM1Z8 < handle
    properties
        serial_number char
        time_out_ms (1,1) double = 1000 * 50; % In millisecond
        device 
%         LastKnownPosition (1,1) double
    end

    properties(Constant, Hidden)
        ThorlabsLibraryRoot = 'C:\Program Files\Thorlabs\Kinesis'
        NET_Assembly_name_list = {'Thorlabs.MotionControl.DeviceManagerCLI', ...
            'Thorlabs.MotionControl.GenericMotorCLI', ...
            'Thorlabs.MotionControl.KCube.DCServoCLI'};
        InitSettingTime_ms = 5000;
    end
    
    properties(Access=protected)
    end

    methods
        function obj = ThorlabsPRM1Z8(serial_num)
            % TODO: error handeling - during initialization
            if isnumeric(serial_num)
                obj.serial_number = sprintf('%d', serial_num);
            else 
                validateattributes(serial_num, {'char', 'string'}, {'nonempty'});
                obj.serial_number = serial_num;
            end
            try
                obj.init_net();
                obj.init_dev_handle();
                obj.initDev();
            catch ME
                obj.disconnect();
                rethrow(ME);
            end
        end

        function init_net(obj)
%             net_domain = System.AppDomain.CurrentDomain;
%             ass_list = net_domain.GetAssemblies;
%             ass_list.Get
%             for i = 1 : numel(obj.NET_Assembly_name_list)
%                ass_fp = fullfile(obj.ThorlabsLibraryRoot, obj.NET_Assembly_name_list{i});
%                if ~exist(ass_fp, 'class')
%                    NET.addAssembly(sprintf('%s.dll', ass_fp));
%                end
%             end
            NET.addAssembly(fullfile(obj.ThorlabsLibraryRoot, ...
                'Thorlabs.MotionControl.DeviceManagerCLI.dll'));
            NET.addAssembly(fullfile(obj.ThorlabsLibraryRoot, ...
                'Thorlabs.MotionControl.GenericMotorCLI.dll'));
            NET.addAssembly(fullfile(obj.ThorlabsLibraryRoot, ...
                'Thorlabs.MotionControl.KCube.DCServoCLI.dll'));
        end

        function init_dev_handle(obj)
            import Thorlabs.MotionControl.DeviceManagerCLI.*
            import Thorlabs.MotionControl.GenericMotorCLI.*
            import Thorlabs.MotionControl.KCube.DCServoCLI.*

            %Initialize Device List
            DeviceManagerCLI.BuildDeviceList();
            DeviceManagerCLI.GetDeviceListSize();
            % Set up device and configuration
            obj.device = KCubeDCServo.CreateKCubeDCServo(obj.serial_number);
            obj.device.Connect(obj.serial_number);
            try
                if ~obj.device.IsSettingsInitialized()
                    obj.device.WaitForSettingsInitialized(obj.InitSettingTime_ms);
                end
                % configure the stage
                motorSettings = obj.device.LoadMotorConfiguration(obj.serial_number);
                motorSettings.DeviceSettingsName = 'PRM1-Z8';
                % update the RealToDeviceUnit converter
                motorSettings.UpdateCurrentConfiguration();
                
                % push the settings down to the device
                MotorDeviceSettings = obj.device.MotorDeviceSettings;
                % Settings
%                 mot_dir = Thorlabs.MotionControl.GenericMotorCLI.Settings.RotationDirection.Quickest;
%                 MotorDeviceSettings.Rotation.RotationDirection = mot_dir;
%                 rot_range = Thorlabs.MotionControl.GenericMotorCLI.Settings.RotationMode.RotationalUnlimited;
%                 MotorDeviceSettings.Rotation.RotationMode = rot_range;
                obj.device.SetSettings(MotorDeviceSettings, true, false);
            catch ME
                obj.device.Disconnect();
                rethrow(ME);
            end
        end

        function initDev(obj)
            import Thorlabs.MotionControl.DeviceManagerCLI.*
            import Thorlabs.MotionControl.GenericMotorCLI.*
            import Thorlabs.MotionControl.KCube.DCServoCLI.*
            obj.device.StartPolling(250);
            pause(1.0); %wait to make sure device is enabled
            try
                if obj.device.NeedsHoming
                    obj.device.Home(obj.time_out_ms);
                end                
            catch ME
                obj.disconnect();
                rethrow(ME);
            end
        end

        function disconnect(obj)
            obj.device.StopPolling();
            obj.device.Disconnect();
        end

        function delete(obj)
            obj.disconnect();
        end
    end
    %%
    methods
        function home(obj)
            obj.device.Home(obj.time_out_ms);
        end

        function moveTo(obj, angle_deg)
            if ~obj.device.IsEnabled
                obj.device.EnableDevice;
            end            
            obj.device.MoveTo(angle_deg, obj.time_out_ms);
        end
        
        function moveRelative(obj, delta_angle_deg)
            if ~obj.device.IsEnabled
                obj.device.EnableDevice;
            end
            delta_angle_deg = System.Decimal(delta_angle_deg);
            obj.device.SetMoveRelativeDistance(delta_angle_deg);
            obj.device.MoveRelative(obj.time_out_ms);
%             obj.device.MoveRelative(delta_angle_deg, obj.time_out_ms);
        end

        function pos = getPosition(obj)
            pos = System.Decimal.ToDouble(obj.device.Position);
%             obj.LastKnownPosition = pos;
        end
    end
end