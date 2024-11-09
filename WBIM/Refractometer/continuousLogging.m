date_str = string(datetime('today','Format','yyyy_MM_dd'));

log_folder = fullfile(pwd, 'logTests', date_str,filesep);

% create log folder if it does not exist yet
if ~exist(log_folder,'dir')
    mkdir(log_folder)
end

% create log file and open for appending
log_fp = fullfile(log_folder,strcat(date_str,".log"));
%
% if isfile(log_fp)
%     delete(log_fp)
% end

fid = fopen(log_fp,'a+');

myRefractometer = InlineRefractometer;

tmax = 0.05; % duration of run in hours
tmax = tmax*60*60; % duration of run in seconds
wait_time = 10; % wait time in seconds

time_total = 0;
tic
i = 0;
while i<=(tmax/wait_time)

    if (i==0 || toc>wait_time)
        time_total = time_total+toc;
        tic

        myRefractometer.measure_RI;
        
        time_total_temp = time_total+toc;
        
        fprintf(fid, '%s ..... t: %s ..... RI: %s ..... BRIX: %s\n', ...
            string(datetime('now','Format','yyyy-MM-dd HH:mm:ss')), ...
            num2str(time_total_temp,'%0.3f'), ...
            num2str(myRefractometer.RI,'%.4f'), ...
            num2str(myRefractometer.BRIX,'%0.1f'));
        i=i+1;
    end

end

myRefractometer.delete