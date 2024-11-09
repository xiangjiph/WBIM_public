classdef WBIMPowerMeter < handle

properties
    meter ThorlabsPowerMeter
end

properties
    device_name = ''
    avg_time_s = 0.1
    timeout_time_s = 1e3
    meter_attenuation = 0;
    wavelength_nm = 800;
    measurement_unit = 'W';
end

methods
    function obj = WBIMPowerMeter(device_name)
        obj.device_name = device_name;
        meterLst = ThorlabsPowerMeter;
        meterIdx = find(strcmp(meterLst.modelName, device_name));
        if ~isempty(meterIdx)
            if meterLst.DeviceAvailable(meterIdx)
                obj.meter = meterLst.connect(meterLst.listdevices, meterIdx);
                obj.init_meter();
            else
                fprintf('%s device can be found but is not available.\n', device_name);
            end
        else
            fprintf('%s device cannot be found.\n', device_name);
        end
        meterLst.delete();
    end
    
    function init_meter(obj)
        obj.meter.setPowerAutoRange(true);
        obj.meter.setAverageTime(obj.avg_time_s);
        obj.meter.setTimeout(obj.timeout_time_s);
        obj.meter.setWaveLength(obj.wavelength_nm);
        % set power measurement unit to W
        obj.set_measurement_unit('W');
    end
    
    function delete(obj)
        obj.meter.disconnect;
        obj.meter.delete;
    end
    
    function power_W = measure_power(obj, num_measurements, wait_time_s)
        arguments
            obj WBIMPowerMeter
            num_measurements {mustBePositive, mustBeInteger} = 1
            wait_time_s {mustBeNonnegative, mustBeNumeric} = 0
        end
        obj.set_measurement_unit(obj.measurement_unit);
        power_W = zeros(num_measurements, 1);
        for i = 1 : num_measurements
            % The response time is slightly longer than the averaging time
            % for the measurement. But what's the sampling rate? 
            [~, power_W(i)] = obj.meter.deviceNET.measPower;
            if wait_time_s
                pause(wait_time_s);
            end
        end        
    end    
end

methods
    function set.wavelength_nm(obj, val)
        obj.meter.setWaveLength(val); %#ok<MCSUP>
        obj.wavelength_nm = val;
    end
    
    function set.avg_time_s(obj, val)
        obj.meter.setAverageTime(val); %#ok<MCSUP>
        obj.avg_time_s = val;
    end
    
    function set_measurement_unit(obj, val)
        arguments
            obj
            val = 'W'
        end
        switch val
            case 'W'
                obj.meter.deviceNET.setPowerUnit(0);
            otherwise 
                error('To be implemented');
        end
    end
end


methods(Static)
    
    
end
    
end