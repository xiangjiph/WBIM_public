set_env;
%% Long-term stability test
% Setting
test_duration_s = 3600 * 72;
record_period_s = 60;
time_step_s = 1;
task_name = sprintf('Stability_test_%s', datestr(now, 'yyyymmdd'));
%%
fprintf('%s Connecting to the laser\n', datestr(now, 'yyyymmdd hh:MM:ss'));
im_laser_h = CoherentChameleon("COM7", 19200);
im_laser_h.set_logger(task_name, record_period_s, time_step_s, test_duration_s);
matlab_window_log_fp = fullfile(im_laser_h.state_logger.root_folder, 'matlab_log.txt');
diary matlab_window_log_fp
diary on
im_laser_h.verbose_Q = false;
im_laser_h.open_port();
if im_laser_h.get_soft_key_state()
    im_laser_h.set_soft_key(0);
    pause(1);
end
%%
im_laser_h.set_soft_key(1);
im_laser_h.set_wavelength_nm(900);
if ~im_laser_h.get_tunable_shutter_state()
    im_laser_h.open_tunable_output_shutter();
    pause(1);
end

try
    im_laser_h.run_logger(true);
catch ME
    warning('%s Unexpected failure. Error message: %s\n', ...
        datestr(now, 'yyyymmdd hh:MM:ss'), getReport(ME, 'extended', 'hyperlinks', 'off'));
end
im_laser_h.close_tunable_output_shutter();
im_laser_h.set_soft_key(0);
im_laser_h.close_port();
fprintf('%s Finish recording laser states\n', datestr(now, 'yyyymmdd hh:MM:ss'));
diary off

%% View history
% im_laser_h.write_state_log(fp, im_laser_record);
% tmp = readtable(fp);
%% Test timmer
