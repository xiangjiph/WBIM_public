% clear all;
NET.addAssembly('C:\Program Files (x86)\Microsoft.NET\Primary Interop Assemblies\Thorlabs.TLPM_64.Interop.dll');

import Thorlabs.TLPM_64.Interop.*;

%Create a dummy TLPM object to check for compatible devices.
handle = System.IntPtr(0);
device = TLPM(handle);

%Look for connected devices
try
    [~,deviceCount] = device.findRsrc();
catch
    print 'Unable to find compatible connected devices. Is the device connected, on, and using the TLPM driver? This example will not work with the legacy IVI drivers.';
    return
end

if deviceCount == 0
    print 'Unable to find compatible connected devices. Is the device connected, on, and using the TLPM driver? This example will not work with the legacy IVI drivers.';
    return
end

%if only one device is connected, connect to this one.
deviceName=System.Text.StringBuilder(256);
if deviceCount == 1
    device.getRsrcName(0, deviceName);
%if multiple are connected, ask which to use
else
    for i = 0:deviceCount-1
        device.getRsrcName(i, deviceName);
        disp(' ');
        disp(['Device #', num2str(i)]);
        disp(deviceName.ToString);
    end
    val=input('Select a device by the number from the above detected devices: ');
    if (floor(val) < deviceCount) && (floor(val) > -1)
        device.getRsrcName(floor(val), deviceName);
    else
        print 'Invalid selection';
        return;
    end
end

deviceName.ToString()

%reassign device to the selected power meter console
device = TLPM(deviceName.ToString(), true, false);

%set wavelength
wavelength= 800;
try
    device.setWavelength(wavelength);
    disp('Wavelength Setting (nm):')
    disp(wavelength);
catch ME
    disp('Error Setting Wavelength');
    device.Dispose()
%     rethrow(ME);
end
    
%measure power
try
    [~, power]=device.measPower();
    disp('Measured Power:');
    disp(power);
catch ME
    disp('Error Measuring Power');
    device.Dispose()
%     rethrow(ME);
end

device.Dispose()
