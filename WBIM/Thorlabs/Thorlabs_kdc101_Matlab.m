%Example for programming the Thorlabs KDC101 with Kinesis in MATLAB, with PRM1-Z8 stage.


%Load assemblies
NET.addAssembly('C:\Program Files\Thorlabs\Kinesis\Thorlabs.MotionControl.DeviceManagerCLI.dll');
NET.addAssembly('C:\Program Files\Thorlabs\Kinesis\Thorlabs.MotionControl.GenericMotorCLI.dll');
NET.addAssembly('C:\Program Files\Thorlabs\Kinesis\Thorlabs.MotionControl.KCube.DCServoCLI.dll');

%Initialize Device List
import Thorlabs.MotionControl.DeviceManagerCLI.*
import Thorlabs.MotionControl.GenericMotorCLI.*
import Thorlabs.MotionControl.KCube.DCServoCLI.*

%Initialize Device List
DeviceManagerCLI.BuildDeviceList();
DeviceManagerCLI.GetDeviceListSize();

%Should change the serial number below to the one being used.
serial_num='27255786';
timeout_val=60000;

%Set up device and configuration
device = KCubeDCServo.CreateKCubeDCServo(serial_num);
device.Connect(serial_num);
device.WaitForSettingsInitialized(5000);


% configure the stage
motorSettings = device.LoadMotorConfiguration(serial_num);
motorSettings.DeviceSettingsName = 'PRM1-Z8';
% update the RealToDeviceUnit converter
motorSettings.UpdateCurrentConfiguration();
 
% push the settings down to the device
MotorDeviceSettings = device.MotorDeviceSettings;
device.SetSettings(MotorDeviceSettings, true, false);


device.StartPolling(250);

pause(1); %wait to make sure device is enabled

%Home
device.Home(timeout_val);
fprintf('Motor homed.\n');

%Move to unit 100
device.MoveTo(100, timeout_val);

%Check Position
pos = System.Decimal.ToDouble(device.Position);
fprintf('The motor position is %d.\n',pos);

device.StopPolling()
device.Disconnect()
%% Test API
h_hwp_abl = WBIMPowerHWP(WBIMConfig.ABLATION_HWP_STAGE_SERIAL_NUMBER, ...
    WBIMConfig.ABLATION_HWP_DELTA_THETA_deg);

theta_list = -30 : 1 : 90;
f_list = h_hwp_abl.get_fractional_power_at_angle(theta_list);
theta_inv_list = h_hwp_abl.get_angle_at_fractional_power(f_list);
fig_hdl = figure;
ax_hdl_1 = subplot(1,2,1);
plot(ax_hdl_1, theta_list, f_list);
grid(ax_hdl_1, 'on');
ax_hdl_1.XLabel.String = '\theta (deg)';
ax_hdl_1.YLabel.String = 'Fraction energy';
ax_hdl_2 = subplot(1,2,2);
plot(ax_hdl_2, theta_list, theta_inv_list);
grid(ax_hdl_2, 'on');