%% Laser stablity test log visualization
task_name = 'Stability_test_20210903';
DataManager = WBIMFileManager;
log_folder = fullfile(DataManager.SCRATCH_ROOT_PATH, ...
    'Imaging_laser_state', task_name);
assert(isfolder(log_folder), 'Folder does not exist');
%% Get all the 
log_list = dir(fullfile(log_folder, sprintf('%s_*.txt', task_name)));
log_list_name = {log_list.name};
log_time = {log_list.date};
log_time = datetime(log_time);
[log_time, log_time_idx] = sortrows(log_time.');
log_list_name = log_list_name(log_time_idx);
num_log_files = numel(log_list);
log_table_cell = cell(num_log_files, 1);
tic_read = tic;
for iter_cell = 1 : num_log_files
   tmp_fp = fullfile(log_folder, log_list_name{iter_cell});
   tmp_table = readtable(tmp_fp); 
   log_table_cell{iter_cell} = tmp_table;
end
fprintf('Finish reading all the log files in the folder. Elapsed time is %.2f seconds\n', ...
    toc(tic_read));
log_table = cat(1, log_table_cell{:});
log_table.Time = datenum(log_table.Time, 'yyyymmdd hh:MM:ss');
log_table.Time = (log_table.Time - log_table.Time(1));
%%
if all(log_table.SoftKeyState)
    fprintf('Soft key state was always on.\n');
else
    warning('Soft key was transiently turned off.');
end
if all(log_table.ShutterState)
    fprintf('Shutter was always opened.\n');
else
    warning('Shutter was transiently closed.');
end
%%
fig_hdl = figure;
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* [2,1.5];
ax_hdl_1 = subplot(3,1,1);
plot(ax_hdl_1, log_table.Time, log_table.Power_mW);
ax_hdl_1.XLabel.String = 'Time (day)';
ax_hdl_1.YLabel.String = 'Power (mW)';
grid(ax_hdl_1, 'on');

ax_hdl_2 = subplot(3,1,2);
plot(ax_hdl_2, log_table.Time, log_table.TunningState);
ax_hdl_2.YLabel.String = 'Tuning state';
ax_hdl_2.XLabel.String = 'Time (day)';
linkaxes([ax_hdl_1, ax_hdl_2], 'x');
ax_hdl_1.XLim(2) = log_table.Time(end);
grid(ax_hdl_2, 'on');

ax_hdl_3 = subplot(3,1,3);
t_s = log_table.Time * 24 * 3600;
plot(ax_hdl_3, t_s, log_table.Power_mW);
ax_hdl_3.XScale = 'log';
ax_hdl_3.XLabel.String = 'Time (second)';
grid(ax_hdl_3, 'on');
ax_hdl_3.XLim(2) = t_s(end);
ax_hdl_3.YLabel.String = deal('Power (mW)');

vis_fp = fullfile(log_folder, sprintf('%s_log_vis.png', task_name));
fun_print_image_in_several_formats(fig_hdl, vis_fp );