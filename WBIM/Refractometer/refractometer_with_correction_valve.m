% File logging
date_str = string(datetime('today','Format','yyyy_MM_dd'));

log_folder = fullfile(pwd, 'logTests', date_str,filesep);

% create log folder if it does not exist yet
if ~exist(log_folder,'dir')
    mkdir(log_folder)
end

% create log file and open for appending
log_fp = fullfile(log_folder,strcat(date_str,".log"));

fid = fopen(log_fp,'a+');

%% Initialize objects
% Begin serial communications with Arduino
% myArduino = serialport('/dev/ttyACM0',9600);
rs = dabs.resources.ResourceStore();

digitalPins = struct;
digitalPins.LCD_ON = rs.filterByName('/vDAQ0/D1.7');
digitalPins.BRIX_READ = rs.filterByName('/vDAQ0/D3.6');
digitalPins.RED_LED = rs.filterByName('/vDAQ0/D1.6');
digitalPins.LCD_STATE = rs.filterByName('/vDAQ0/D0.7');
digitalPins.VALVE = rs.filterByName('/vDAQ0/D3.7');

% Initialize refractometer object with Arduino
% myRefractometer = InlineRefractometer(myArduino);
myRefractometer = InlineRefractometer(digitalPins);

% Initialize solenoid valve object with Arduino
% myValve = SolenoidValve(myArduino);
myValve = SolenoidValve(digitalPins);
myValve.valve_timer(5000); % open valve for 5s to get initial volume
myValve.v_5s_mL = 7.5; % set initial volume in mL

%%
set_RI = 1.4201;

tmax = 6; % duration of run in hours
tmax = tmax*60*60; % duration of run in seconds
wait_time = 600; % wait time in seconds

time_total = 0;
tic
i = 0;
while i<=(tmax/wait_time)

    if (i==0 || toc>wait_time)
        time_total = time_total+toc;
        tic

        myRefractometer.measure_RI;
        wait(myRefractometer.timerMeasure)
        
        time_total_temp = time_total+toc;
        
        fprintf(fid, '%s ..... t: %s ..... RI: %s ..... BRIX: %s\n', ...
            string(datetime('now','Format','yyyy-MM-dd HH:mm:ss')), ...
            num2str(time_total_temp,'%0.3f'), ...
            num2str(myRefractometer.RI,'%.4f'), ...
            num2str(myRefractometer.BRIX,'%0.1f'));
        i=i+1;

        while myRefractometer.RI > set_RI
            % Add volume
            myValve.valve_volume(0.5); % add volume in mL 
            pause(30);
            try
                myRefractometer.measure_RI;
                wait(myRefractometer.timerMeasure);
                time_total_temp = time_total+toc;
        
                fprintf(fid, '%s ..... t: %s ..... RI: %s ..... BRIX: %s\n', ...
                    string(datetime('now','Format','yyyy-MM-dd HH:mm:ss')), ...
                    num2str(time_total_temp,'%0.3f'), ...
                    num2str(myRefractometer.RI,'%.4f'), ...
                    num2str(myRefractometer.BRIX,'%0.1f'));
            catch ME
            end
        end
    end

end