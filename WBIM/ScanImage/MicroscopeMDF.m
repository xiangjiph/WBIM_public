% Most Software Machine Data File

%% scanimage.SI (ScanImage)

% Global microscope properties
objectiveResolution = 40.48;     % Resolution of the objective in microns/degree of scan angle

% Data file location

% Custom Scripts
startUpScript = '';     % Name of script that is executed in workspace 'base' after scanimage initializes
shutDownScript = '';     % Name of script that is executed in workspace 'base' after scanimage exits

fieldCurvatureZs = [];     % Field curvature for mesoscope
fieldCurvatureRxs = [];     % Field curvature for mesoscope
fieldCurvatureRys = [];     % Field curvature for mesoscope
useJsonHeaderFormat = false;     % Use JSON format for TIFF file header

fieldCurvatureTip = 0;
fieldCurvatureTilt = 0;

%% scanimage.components.Motors (SI Motors)
% SI Stage/Motor Component.
motorXYZ = {'Stages' 'Stages' 'Stages'};     % Defines the motor for ScanImage axes X Y Z.
motorAxisXYZ = [1 2 3];     % Defines the motor axis used for Scanimage axes X Y Z.
scaleXYZ = [1 1 1];     % Defines scaling factors for axes.
backlashCompensation = [0 0 0];     % Backlash compensation in um (positive or negative)

%% scanimage.components.Photostim (SI Photostim)
photostimScannerName = '';     % Name of scanner (from first MDF section) to use for photostimulation. Must be a linear scanner

% Monitoring DAQ AI channels
BeamAiId = [];     % AI channel to be used for monitoring the Pockels cell output

loggingStartTrigger = '';     % PFI line to which start trigger for logging is wired to photostim board. Leave empty for automatic routing via PXI bus

stimActiveOutputChannel = '';     % Digital terminal on stim board to output stim active signal. (e.g. on vDAQ: 'D2.6' on NI-DAQ hardware: '/port0/line0'
beamActiveOutputChannel = '';     % Digital terminal on stim board to output beam active signal. (e.g. on vDAQ: 'D2.7' on NI-DAQ hardware: '/port0/line1'
slmTriggerOutputChannel = '';     % Digital terminal on stim board to trigger SLM frame flip. (e.g. on vDAQ: 'D2.5' on NI-DAQ hardware: '/port0/line2'

%% dabs.generic.GalvoPureAnalog (Galvo)
AOControl = '/vDAQ0/AO2';     % control terminal  e.g. '/vDAQ0/AO0'
AOOffset = '';     % control terminal  e.g. '/vDAQ0/AO0'
AIFeedback = '/vDAQ0/AI0';     % feedback terminal e.g. '/vDAQ0/AI0'

angularRange = 20;     % total angular range in optical degrees (e.g. for a galvo with -20..+20 optical degrees, enter 40)
voltsPerOpticalDegrees = 1.0239;     % volts per optical degrees for the control signal
parkPosition = 0;     % park position in optical degrees
slewRateLimit = Inf;     % Slew rate limit of the analog output in Volts per second

% Calibration settings
feedbackVoltLUT = [-2.4858 -10.239;-1.97414 -7.96367;-1.42371 -5.68833;-0.855936 -3.413;-0.280566 -1.13767;0.29201 1.13767;0.866119 3.413;1.42621 5.68833;2.01035 7.96367;2.5183 10.239];     % [Nx2] lut translating feedback volts into position volts
offsetVoltScaling = 1;     % scalar factor for offset volts

voltsOffset = 0;

%% dabs.zaber.ZaberMultiDevice (Stages)
comPort = 'COM4';     % Serial port the stage is connected to (e.g. 'COM3')
baudRate = 115200;     % Baudrate for serial communication
communicationProtocol = 'ASCII';     % Communication protocol ('ASCII' or 'Binary')
deviceLibraryPath = '';     % Path to '.sqlite' device library. Only required for offline use.
homingTimeout_s = 20;     % Timeout for homing move in seconds

%% dabs.generic.ResonantScannerAnalog (CRS8K)
AOZoom = '/vDAQ0/AO0';     % zoom control terminal  e.g. '/vDAQ0/AO0'
DOEnable = '/vDAQ0/D0.1';     % digital enable terminal e.g. '/vDAQ0/D0.1'
DISync = '/vDAQ0/D0.0';     % digital sync terminal e.g. '/vDAQ0/D0.0'

nominalFrequency = 7910;     % nominal resonant frequency in Hz
angularRange = 26;     % total angular range in optical degrees (e.g. for a resonant scanner with -13..+13 optical degrees, enter 26)
voltsPerOpticalDegrees = 0.1923;     % volts per optical degrees for the control signal
settleTime = 0.5;     % settle time in seconds to allow the resonant scanner to turn on

% Calibration Settings
amplitudeToLinePhaseMap = [0.198 1.5375e-06;0.2 1.775e-06;0.281 1.1125e-06;0.286 1.2625e-06;0.327 1.075e-06;0.392 2.75e-07;0.487 7.5e-07;0.497 8.125e-07;0.667 6.25e-07;1.136 4e-07;1.666 2.88e-07;2 3.125e-07;2.066 3.425e-06;2.261 2.25e-07;2.499 2.48e-07;2.857 2.875e-07;3.004 2e-07;3.333 2.25e-07;3.999 2.5e-07;4.885 2.5e-07;4.999 2.24e-07;6.665 2.375e-07;8.403 2.375e-07;8.638 2.625e-07;10 2.16e-07;10.203 1.92e-07;14.203 2.16e-07;17.066 2.16e-07;17.839 2.375e-07;18.033 2.375e-07;18.226 2.08e-07;18.357 2.125e-07;18.676 2.5e-07;18.706 2.125e-07;19.139 2.25e-07;19.229 1.92e-07;19.258 2.125e-07;19.967 2.625e-07;19.999 2e-07;20.82 2.08e-07;23.154 2.08e-07];     % translates an amplitude (degrees) to a line phase (seconds)
amplitudeToFrequencyMap = [0.279 7924.59;0.645 7936.43;2.001 7936.48;2.501 7937.3;2.81 7936.63;3.004 7928.68;3.006 7937.14;3.333 7928.28;3.606 7938.3;3.607 7937.63;4 7923.21;4.312 7936.74;4.507 7936.49;4.508 7936.83;4.523 7936.42;4.524 7936.48;4.547 7933.06;4.548 7935.76;4.558 7933.39;4.588 7937.19;4.589 7937.1;4.785 7936.52;4.902 7935.83;4.916 7935.18;4.918 7935.96;4.926 7935.8;4.927 7936.52;4.934 7936.98;4.958 7936.53;4.959 7936.57;4.999 7932.35;5.003 7933.11;5.875 7934.92;6.011 7936.33;6.143 7937.12;6.445 7937.06;6.446 7937.47;6.483 7934.98;6.665 7935.52;6.667 7935.76;7.998 7917.5;8 7921.43;8.265 7935.85;8.319 7921.85;8.32 7922.36;8.403 7936.48;8.404 7935.96;8.421 7936.26;8.433 7920.8;8.434 7915.62;8.449 7935.61;8.451 7935.7;8.606 7934.66;8.608 7934.48;8.613 7934.89;8.624 7934.32;8.638 7936.56;8.639 7936.32;8.644 7935.32;8.645 7936.06;8.666 7936.1;8.667 7935.67;8.974 7935.82;8.976 7935.53;8.985 7933.89;8.986 7935.09;9.016 7935.42;9.017 7935.9;9.047 7935.78;9.048 7935.55;9.065 7936.7;9.136 7931.55;9.182 7936.3;9.184 7934.63;9.216 7934.39;9.239 7936.8;9.241 7936.78;9.654 7932.2;10 7920.23;10.203 7921.55;10.204 7921.55;10.255 7921.62;10.256 7922.21;10.953 7935.23;10.954 7934.71;11.049 7935.25;11.093 7920.38;11.11 7922.73;11.25 7936.48;11.748 7936.22;11.749 7933.42;11.975 7921.98;11.976 7919.35;12.891 7934.5;12.893 7934.82;13.332 7918.41;13.426 7918.71;13.428 7919.64;13.871 7919.35;14.203 7917.01;14.204 7918.82;14.276 7921.26;14.277 7917.35;14.48 7922.25;14.545 7919.86;14.927 7921.72;15.583 7919.56;15.999 7921.33;16 7918.64;16.722 7920.6;16.723 7918.75;17.066 7920.58;17.067 7916.15;17.76 7931.63;17.762 7928;17.839 7931.23;17.84 7932.12;17.938 7934.57;17.939 7932.24;18.033 7933.82;18.034 7928.6;18.095 7930.94;18.097 7934.91;18.166 7930.76;18.167 7927.18;18.18 7929.55;18.182 7929.35;18.19 7933.58;18.226 7919.45;18.228 7921;18.357 7931.87;18.358 7930.39;18.676 7927.93;18.692 7928.85;18.706 7933.63;18.714 7934.5;18.741 7917.8;18.779 7934.87;19.139 7934.37;19.141 7933.37;19.164 7930.3;19.165 7935.25;19.196 7934.33;19.198 7930.89;19.212 7932.63;19.229 7920.36;19.23 7917.44;19.259 7933.38;19.341 7918.36;19.418 7915.22;19.419 7916.3;19.519 7919.54;19.779 7920.73;19.945 7932.41;19.946 7935.19;19.967 7934.5;19.968 7934.44;19.999 7925.37;20 7918.32;20.009 7926.73;20.01 7924.29;20.066 7917.93;20.067 7915.5;20.139 7927.52;20.141 7934.84;20.491 7920.62;20.493 7915.67;20.575 7920.59;20.576 7912.27;20.82 7919.95;20.821 7910.15;21.988 7913.66;21.989 7920.05;22.221 7920.76;22.222 7918.18;23.154 7921.02;23.155 7917.07];     % translates an amplitude (degrees) to a resonant frequency (Hz)
amplitudeLUT = zeros(0,2);     % translates a nominal amplitude (degrees) to an output amplitude (degrees)

%% dabs.npoint.LC40x (Piezo)
AOControl = '/vDAQ0/AO1';     % control terminal  e.g. '/vDAQ0/AO0'
AIFeedback = '/vDAQ0/AI1';     % feedback terminal e.g. '/vDAQ0/AI0'
FrameClockIn = '';     % frame clock input terminal e.g. '/Dev1/PFI0'

parkPositionUm = -50;     % park position in micron
travelRangeUm = [-200 200];     % travel range in micron

voltsPerUm = 0.05;     % volts per micron
voltsOffset = 0;     % volts that sets actuator to zero position

% Calibration Data
positionLUT = zeros(0,2);     % Position LUT
feedbackVoltLUT = [-9.86648 -10;-7.67913 -7.77778;-5.49206 -5.55556;-3.30263 -3.33333;-1.10831 -1.11111;1.11912 1.11111;3.30639 3.33333;5.49224 5.55556;7.68315 7.77778;9.86663 10];     % [Nx2] lut translating feedback volts into position volts
comPort = 'COM6';     % Name of the COM port associated with the controller (e.g. 'COM1')
channel = 1;     % Number of the LC40x channel which the device is connected (e.g. '1')

%% scanimage.components.scan2d.RggScan (Microscope)

acquisitionDeviceId = 'vDAQ0';     % RDI Device ID
acquisitionEngineIdx = 1;

resonantScanner = 'CRS8K';     % Name of the resonant scanner
xGalvo = '';     % Name of the x galvo scanner
yGalvo = 'Galvo';     % Name of the y galvo scanner
beams = {'Imaging Pockels Cell'};     % beam device names
fastZs = {'Piezo'};     % fastZ device names
shutters = {'Imaging Shutter'};     % shutter device names

channelsInvert = [false false false false];     % Logical: Specifies if the input signal is inverted (i.e., more negative for increased light signal)
keepResonantScannerOn = false;     % Always keep resonant scanner on to avoid drift and settling time issues

externalSampleClock = false;     % Logical: use external sample clock connected to the CLK IN terminal of the FlexRIO digitizer module
externalSampleClockRate = 1.25e+08;     % [Hz]: nominal frequency of the external sample clock connected to the CLK IN terminal (e.g. 80e6); actual rate is measured on FPGA
externalSampleClockMultiplier = 1;     % Multiplier to apply to external sample clock

extendedRggFov = 0;     % If true and x galvo is present, addressable FOV is combination of resonant FOV and x galvo FOV.

% Advanced/Optional
PeriodClockDebounceTime = 1e-07;     % [s] time the period clock has to be stable before a change is registered
TriggerDebounceTime = 5e-07;     % [s] time acquisition, stop and next trigger to be stable before a change is registered
reverseLineRead = 0;     % flips the image in the resonant scan axis
defaultFlybackTimePerFrame = 0.001;     % [s] default time to allow galvos to fly back after one frame is complete. overridden by cfg file
defaultFlytoTimePerScanfield = 0.001;     % [s] time to allow galvos to fly from one scanfield to the next. overridden by cfg file

% Aux Trigger Recording, Photon Counting, and I2C are mutually exclusive

% Aux Trigger Recording
auxTriggersTimeDebounce = 1e-07;     % [s] time after an edge where subsequent edges are ignored
auxTriggerLinesInvert = [false false false false];     % [logical] 1x4 vector specifying polarity of aux trigger inputs
auxTrigger1In = '';     % Digital input lines for aux trigger 1
auxTrigger2In = '';     % Digital input lines for aux trigger 2
auxTrigger3In = '';     % Digital input lines for aux trigger 3
auxTrigger4In = '';     % Digital input lines for aux trigger 4

% Signal Conditioning
disableMaskDivide = [false false false false];     % disable averaging of samples into pixels; instead accumulate samples
photonDiscriminatorThresholds = [500 500];
photonDiscriminatorModes = {'threshold crossing' 'threshold crossing'};
photonDiscriminatorDifferentiateWidths = [4 4];

% I2C
i2cEnable = false;
i2cSdaPort = '';
i2cSclPort = '';
i2cAddress = 0;     % [byte] I2C address of the FPGA
i2cDebounce = 1e-07;     % [s] time the I2C signal has to be stable high before a change is registered
i2cStoreAsChar = false;     % if false, the I2C packet bytes are stored as a uint8 array. if true, the I2C packet bytes are stored as a string. Note: a Null byte in the packet terminates the string
i2cSendAck = true;     % When enabled FPGA confirms each packet with an ACK bit by actively pulling down the SDA line

% Laser Trigger
LaserTriggerPort = '/vDAQ0/CLK_IN';     % Digital input where laser trigger is connected.

% Trigger Outputs
frameClockOut = '';     % Output line for the frame clock
lineClockOut = '';     % Output line for the line clock
beamModifiedLineClockOut = '';     % Output line for beam clock
volumeTriggerOut = '';     % Output line for the volume clock

% Calibration data
scannerToRefTransform = [1 0 0;0 1 0;0 0 1];
LaserTriggerDebounceTicks = 1;

virtualChannelsSource = {'AI0' 'AI1' 'AI2' 'AI3'};
virtualChannelsMode = {'analog' 'analog' 'analog' 'analog'};
virtualChannelsThreshold = [false false false false];
virtualChannelsBinarize = [false false false false];
virtualChannelsEdgeDetect = [false false false false];
virtualChannelsLaserGate = [false false false false];
virtualChannelsDisableDivide = [false false false false];
virtualChannelsThresholdValue = [100 100 100 100];
virtualChannelsLaserFilterWindow = {[0 1] [0 1] [0 1] [0 1]};

enableHostPixelCorrection = false;
hostPixelCorrectionMultiplier = 500;
sampleClockPhase = [];

desiredSampleRateCtl = 1e+06;
useCustomFilterClock = false;
customFilterClockPeriod = 1;

%% dabs.generic.BeamModulatorFastAnalog (Imaging Pockels Cell)
AOControl = '/vDAQ0/AO4';     % control terminal  e.g. '/vDAQ0/AO0'
AIFeedback = '';     % feedback terminal e.g. '/vDAQ0/AI0'

outputRange_V = [0 1.075];     % Control output range in Volts
feedbackUsesRejectedLight = false;     % Indicates if photodiode is in rejected path of beams modulator.
calibrationOpenShutters = {'Imaging Shutter'};     % List of shutters to open during the calibration. (e.g. {'Shutter1' 'Shutter2'}

powerFractionLimit = 0.8;     % Maximum allowed power fraction (between 0 and 1)

% Calibration data
powerFraction2ModulationVoltLut = [0.00242215 0;0.00426759 0.025;0.00882353 0.05;0.0156286 0.075;0.0243368 0.1;0.0369666 0.125;0.0511534 0.15;0.0677047 0.175;0.0866782 0.2;0.107555 0.225;0.130623 0.25;0.155421 0.275;0.181776 0.3;0.211765 0.325;0.242272 0.35;0.274221 0.375;0.308535 0.4;0.341984 0.425;0.376586 0.45;0.411765 0.475;0.44752 0.5;0.483852 0.525;0.520761 0.55;0.557093 0.575;0.592849 0.6;0.627451 0.625;0.662053 0.65;0.695502 0.675;0.72895 0.7;0.760669 0.725;0.790657 0.75;0.818339 0.775;0.846021 0.8;0.871972 0.825;0.895617 0.85;0.917532 0.875;0.936563 0.9;0.95271 0.925;0.967128 0.95;0.979815 0.975;0.988466 1;0.995386 1.025;0.998847 1.05;1 1.075];
powerFraction2PowerWattLut = [0 0.0042;1 1.734];
powerFraction2FeedbackVoltLut = [0 0.0460709;1 9.99975];
feedbackOffset_V = 0;

% Calibration settings
calibrationNumPoints = 200;     % number of equidistant points to measure within the analog ouptut range
calibrationAverageSamples = 5;     % per analog output voltage, average N analog input samples. This helps to reduce noise
calibrationNumRepeats = 5;     % number of times to repeat the calibration routine. the end result is the average of all calibration runs
calibrationSettlingTime_s = 0.01;     % pause between measurement points. this allows the beam modulation to settle
calibrationFlybackTime_s = 0.2;     % pause between calibration runs

% Advanced Settings. Note: these settings are unused for vDAQ based systems
modifiedLineClockIn = '';     % Terminal to which external beam trigger is connected. Leave empty for automatic routing via PXI/RTSI bus
frameClockIn = '';     % Terminal to which external frame clock is connected. Leave empty for automatic routing via PXI/RTSI bus
referenceClockIn = '';     % Terminal to which external reference clock is connected. Leave empty for automatic routing via PXI/RTSI bus
referenceClockRate = 1e+07;     % if referenceClockIn is used, referenceClockRate defines the rate of the reference clock in Hz. Default: 10e6Hz

%% dabs.generic.DigitalShutter (Imaging Shutter)
DOControl = '/vDAQ0/D0.2';     % control terminal  e.g. '/vDAQ0/DIO0'
invertOutput = true;     % invert output drive signal to shutter
openTime_s = 0.1;     % settling time for shutter in seconds

%% dabs.generic.DigitalShutter (Ablation Shutter)
DOControl = '/vDAQ0/D0.3';     % control terminal  e.g. '/vDAQ0/DIO0'
invertOutput = true;     % invert output drive signal to shutter
openTime_s = 0.1;     % settling time for shutter in seconds

%% scanimage.components.CoordinateSystems (SI CoordinateSystems)
% SI Coordinate System Component.
classDataFileName = 'default-CoordinateSystems_classData.mat';     % File containing the previously generated alignment data corresponding to the currently installed objective, SLM, scanners, etc.

