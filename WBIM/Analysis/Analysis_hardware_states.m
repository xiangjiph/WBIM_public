clc;clear;close all;
DataManager = WBIMFileManager;
exp_group = 'Ablation_test';
exp_name = 'WBIM20240305002_20240801';
acq_root_folder = DataManager.fp_experiment(exp_group, exp_name);
vis_folder = fullfile(acq_root_folder, 'visualization', 'log');
sub_folders = dir(acq_root_folder);
sub_folders = arrayfun(@(x) x.name, sub_folders, 'UniformOutput', false);
is_layer_folder_Q = ~cellfun(@isempty, regexp(sub_folders, '[0-9]{5}'));
sub_folders = sort(cellfun(@str2double, sub_folders(is_layer_folder_Q)));
% Load setting file 
acq_info = load(DataManager.fp_acq_control_parameter(exp_group, exp_name, 0));
%%
num_layer = numel(sub_folders);
laser_layer_log = cell(num_layer, 1);
for i = 1 : num_layer
    tmp_layer = sub_folders(i);
    tmp_layer_fp = fullfile(DataManager.fp_layer(exp_group, exp_name, tmp_layer), ...
        'log');
    tmp_laser_log = fullfile(tmp_layer_fp, DataManager.fn_imaging_laser_log_file(...
        exp_group, exp_name, tmp_layer));
    tmp_text = readlines(tmp_laser_log);
    selected_Q = ~contains(tmp_text, {'WARNING', 'INFO', 'Mean', 'Basic statistics', 'Start', 'Stop'}) & ...
        cellfun(@numel, tmp_text) > 1;
    tmp_text = tmp_text(selected_Q);
    tmp_log = arrayfun(@(x) strsplit(x, {' ', '\t'}), tmp_text, 'UniformOutput', false);
    tmp_log = cat(1, tmp_log{:});
    laser_layer_log{i} = tmp_log;
end
%%
laser_layer_log = cat(1, laser_layer_log{:});
laser_log = struct;
laser_log.power_mW = str2double(cat(1, laser_layer_log(:, 13)));
laser_log.time_s = cellfun(@(x1, x2) [x1, ' ', x2], laser_layer_log(:, 1), laser_layer_log(:, 2), ...
    'UniformOutput', false);
laser_log.time_s = datetime(laser_log.time_s);
laser_log_t0 = laser_log.time_s(1);
laser_log.time_s = laser_log.time_s - laser_log_t0;

t_max = laser_log.time_s(end);
fig_hdl = figure;
ax = axes(fig_hdl);
plot(ax, laser_log.time_s, laser_log.power_mW);
ax.YLim = [2700, 3000];
ax.XLabel.String = 'Time';
ax.YLabel.String = 'Imaging laser output power (mW)';
fig_fp = fullfile(vis_folder, sprintf('%s_%s_laser_power_log.png', exp_group, ...
    exp_name));
fun_print_image_in_several_formats(fig_hdl, fig_fp);
%% 
num_layer = numel(sub_folders);
layer_log = cell(num_layer, 1);
for i = 1 : num_layer
    tmp_layer = sub_folders(i);
    tmp_layer_fp = fullfile(DataManager.fp_layer(exp_group, exp_name, tmp_layer), ...
        'log');
    tmp_machine_log = fullfile(tmp_layer_fp, DataManager.fn_microscope_log_file(...
        exp_group, exp_name, tmp_layer));
    tmp_text = readlines(tmp_machine_log);

    tmp_info = struct;
    % Refractive index measurement
    tmp_info.RI = tmp_text(contains(tmp_text, {'MESSAGE RI'}));
    tmp_info.add_water = tmp_text(contains(tmp_text, {'mL of water'}));
    tmp_info.pump_on = tmp_text(contains(tmp_text, 'Turn on pump'));
    tmp_info.pump_off = tmp_text(contains(tmp_text, 'Turn off pump'));
    tmp_info.set_ablation_power = tmp_text(contains(tmp_text, 'Set ablation beam fractional power to be'));
    tmp_info.ablation_done = tmp_text(contains(tmp_text, 'Ablation beam external shutter closed'));
    tmp_info.ablation_start = tmp_text(contains(tmp_text, 'Ablation beam external shutter opened'));
    tmp_info.imaging_laser = tmp_text(contains(tmp_text, 'laser'));
    tmp_info.row_shift = tmp_text(contains(tmp_text, 'Detected row shift in range'));
    tmp_info.reimage = tmp_text(contains(tmp_text, 'Re-image'));
    layer_log{i} = tmp_info;
end
layer_log = cat(1, layer_log{:});
%% Parse log
machine_log = struct;
fn = fieldnames(layer_log);
for i = 1 : numel(fn)
    tmp_fn = fn{i};
    tmp_log = cat(1, layer_log.(tmp_fn));
    tmp_t = WBIMUtil.get_log_timestamp(tmp_log, laser_log_t0);
    tmp_log = arrayfun(@(x) strsplit(x, {' ', '\t'}), tmp_log, 'UniformOutput', false);
    switch tmp_fn
        case 'set_ablation_power'
            abl_pwr = cellfun(@(x) str2double(x{11}(1:end-1)), tmp_log);
            machine_log.abl_pwr_f = abl_pwr / 100;
        case 'RI'
            ri_val = cellfun(@(x) str2double(x{5}), tmp_log);
            [tmp_t, tmp_idx] = sort(tmp_t);
            machine_log.RI = ri_val(tmp_idx);
        case 'add_water'
            machine_log.water_vol = cellfun(@(x) str2double(x{5}), tmp_log);

    end
    machine_log.(sprintf('t_%s', tmp_fn)) = tmp_t;
end
%%
fig_hdl = figure;
a1 = subplot(3, 1, 1);
plot(a1, laser_log.time_s, laser_log.power_mW);
a1.YLim = [2700, 3000];
a1.XLabel.String = 'Time';
a1.YLabel.String = 'Imaging laser output power (mW)';
yyaxis(a1, 'right');
hold(a1, 'on');
scatter(a1, machine_log.t_reimage, ones(size(machine_log.t_reimage)), 'rx');
a1.YLabel.String = 'Re-image';
a1.Title.String = sprintf("Experiment: %s", exp_name);
a1.Title.Interpreter = 'none';

a2 = subplot(3,1,2);
[ri_t, sort_idx] = sort(machine_log.t_RI);
plot(a2, ri_t, machine_log.RI(sort_idx));
a2.YLim = [1.4280, 1.43];
a2.YLabel.String = 'Refractive index';
hold(a2, 'on');
yyaxis(a2, 'right');
scatter(a2, machine_log.t_add_water, machine_log.water_vol , 'o', 'filled');
a2.YLabel.String = 'Added water volume (mL)';

abl_t = cat(1, machine_log.t_ablation_start, machine_log.t_set_ablation_power, ...
    machine_log.t_ablation_done);
abl_pwr = cat(1, zeros(size(machine_log.t_ablation_start)), machine_log.abl_pwr_f, ...
    zeros(size(machine_log.t_ablation_done)));
[abl_t, tmp_idx] = sort(abl_t);
abl_pwr = abl_pwr(tmp_idx);

abl_t = cat(1, abl_t, abl_t(2:end) - 10 * eps);
abl_pwr = cat(1, abl_pwr, abl_pwr(1:end-1));
[abl_t, tmp_idx] = sort(abl_t);
abl_pwr = abl_pwr(tmp_idx);
spot_shape = WBIMSPAblation.estimate_ablation_area(false);
transmitted_W = acq_info.ablation.laser_output_power_W * WBIMConfig.ABLATION_TRANSMISSION_FRACTION * ...
                abl_pwr;  
abl_fluence = (transmitted_W / WBIMConfig.ABLATION_REPETITION_RATE_Hz)...
    / (spot_shape.deff_um_fast * spot_shape.deff_um_slow * pi / 4 / 1e8);

a3 = subplot(3,1,3);
plot(a3, abl_t, abl_fluence);
% scatter(a3, abl_t, abl_pwr * acq_info.ablation.laser_output_power_W * 1e3, 'r.');
a3.YLabel.String = sprintf("Ablation fluence (J/cm^2)");

fig_fp = fullfile(vis_folder, sprintf('%s_%s_exp_log.png', exp_group, ...
    exp_name));
fun_print_image_in_several_formats(fig_hdl, fig_fp);