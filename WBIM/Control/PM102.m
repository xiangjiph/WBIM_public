%% Parameters

%Imaging laser parameters
wavelenArray = [810 950 1029]; %nm

%Meter parameters
meterName = 'PM102';
meterAttenuation = 0;
avgTime = 0.1; %seconds
timeoutTime = 1000; %milliseconds

%Pockel cell parameters
voltLow = 0;
voltHigh = 1.5;
dVolt = 0.05;
voltArray = [voltLow:dVolt:voltHigh]';

% Measurement times
updateTime = 0.2;   %seconds
nSamples = 5;      %number of samples to collect
powerOutVals = zeros(1,nSamples);   %output power in watts
powerInVals = zeros(1,nSamples);    %input power in milliwatts


%% Initialize calibration arrays
calibrationPts = struct;
calibrationPts.biasVoltage = -40;
calibrationPts.wavelengths = wavelenArray;
calibrationPts.voltages = voltArray;
nPoints = length(calibrationPts.voltages);
calibrationPts.measuredOutputPower = zeros(nPoints,length(wavelenArray));
calibrationPts.measuredOutputPowerError = zeros(nPoints,length(wavelenArray));
calibrationPts.measuredInputPower = zeros(nPoints,length(wavelenArray));

%% Connect to pockel cell
rs = dabs.resources.ResourceStore();
pockel_ao = rs.filterByName('/vDAQ0/AO4');

%% Connect to power meter
meterLst = ThorlabsPowerMeter;
meterAddress = meterLst.listdevices;
meterIdx = find(strcmp(meterLst.modelName, meterName));
if ~isempty(meterIdx)
%     meterAddress = meterLst.listdevices{meterIdx};
else
    warning(strcat(meterName, ' could not be found!'))
    return
end

if meterLst.DeviceAvailable(meterIdx)
    meter = meterLst.connect(meterAddress,meterIdx);
else
    warning(strcat(meterName,' is not available!'))
    return
end

%% Initialize power meter parameters
meter.setPowerAutoRange(true);
pause(5)    % Pause the program a bit to allow the power meter to autoadjust
meter.setAverageTime(avgTime);
meter.setTimeout(timeoutTime);

%% Iterate through all wavelengths
for i = 1:length(wavelenArray)

    % Set wavelength parameter
    wavelen = wavelenArray(i);
    wbim.imaging.h_laser.set_wavelength_nm(wavelen)
    meter.setWaveLength(wavelen);
    pause(10)   % give time for the wavelength to adjust
    
    tic
    for j = 1:nPoints
        fprintf('Measuring %1.0fnm at %0.2fV\r',wavelen,voltArray(j))
        pockel_ao.setValue(voltArray(j));
        pause(2)    %give time for the voltage to adjust

        % trackPower = figure;
        wbim.imaging.h_laser.open_tunable_output_shutter;
        pause(10)    %Give time for output shutter to open and meter to reach equilibrium
        % tic
        for k=1:nSamples
            powerInVals(k) = wbim.imaging.h_laser.power_mW;
            meter.updateReading(updateTime);
            powerOutVals(k) = meter.meterPowerReading;            
        %     fprintf('%.10f%c\r',powerOutVals(k),meter.meterPowerUnit);

        %     figure(trackPower)
        %     plot([0:(k-1)]*updateTime,powerOutVals(1:k))
        %     drawnow
        end
        % toc
        % fprintf('%.10f%c\r',mean(powerOutVals),meter.meterPowerUnit);
        % plot([0:(nSamples-1)]*updateTime,powerOutVals)
        calibrationPts.measuredOutputPower(j,i) = mean(powerOutVals);
        calibrationPts.measuredOutputPowerError(j,i) = std(powerOutVals);
        calibrationPts.measuredInputPower(j,i) = mean(powerInVals)/1000;
    end
    toc
    wbim.imaging.h_laser.close_tunable_output_shutter
    pockel_ao.setValue(0);
    if i~=length(wavelenArray)
        pause(300)  %relax for 5 minutes
    end
end

save('calibrationPts_810_950_1029nm_20231104_0V_1d5V_fast.mat','-struct','calibrationPts')

%% Plotting
figure;
hold on
% for i = 1:length(calibrationPts.wavelengths)
%     tempFun = 100*calibrationPts.measuredOutputPower(:,i)./calibrationPts.measuredInputPower(:,i);
%     plot(calibrationPts.voltages,100*tempFun/max(tempFun(:)))
% end

for i = 1:length(calibrationPts.wavelengths)
    tempFun = 100*calibrationPts.measuredOutputPower(:,i)./calibrationPts.measuredInputPower(:,i);
    plot(calibrationPts.voltages,tempFun)
end

lgd = legend(num2str(calibrationPts.wavelengths'));
lgd.Location = 'northwest';

xlabel('Applied Voltages (V)')
ylabel('Intensity Percentages (%)')
title('Calibration Curves')

%% Fitting

mdl = fittype('a*sin(b*x+c)+d','indep','x');
fittedmdl = fit(x,100*(y/max(y(:))),mdl,'start',[max(y(:))/2,2*pi/(2.4),0,max(y(:))/2]);

%% Disconnect
wbim.imaging.h_laser.close_tunable_output_shutter
meter.disconnect
meter.delete
meterLst.delete