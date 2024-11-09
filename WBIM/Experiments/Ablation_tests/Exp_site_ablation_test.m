%% Single site ablation test

%% Initialization 
% Both imaging and ablation 
% 1. Take stack of local region
% 2. Start ablation experiments
% 3. Take stacks
%% Power
wbim.ablation.set_ablation_fractional_power(0.08);
    %% Switch to ablation mode
wbim.ablation.turn_on_pump();
wbim.ablation.actuator_switch_to_mirror();
wbim.ablation.h_zaber_controller.set_ablation_do_value(1);
wbim.ablation.open_ablation_shutter();
%% Parameters
num_pulse_list = [1, 2, 3, 4, 5];
num_trial = 5;
trial_spacing_um = 100;
%% Number of pulses
test_parameter = num_pulse_list;
z_idx = 1;
pulse_energy_idx = 1;
num_p = numel(test_parameter);
total_num_trial = num_p * num_trial;

x_um_vec = (0 : 1 : (num_p - 1)) * trial_spacing_um;
x_um_vec = x_um_vec - x_um_vec(ceil(num_trial/2));
y_um_vec = (0 : 1 : (num_trial - 1)) * trial_spacing_um;
y_um_vec = y_um_vec - y_um_vec(ceil(num_trial/2));
init_pos_um = wbim.sample_xyz_um;
x_um_vec = x_um_vec + init_pos_um(1);
y_um_vec = y_um_vec + init_pos_um(2);
for iter_np = 1 : num_p
    % Generate pulse task 
    tmp_num_pulse = test_parameter(iter_np);
    tmp_output_task = wbim.ablation.setup_ablation_pulses_task(tmp_num_pulse);
    wbim.move_sample_along_axis_um(1, x_um_vec(iter_np), false);
    for iter_t = 1 : num_trial
        wbim.move_sample_along_axis_um(2, y_um_vec(iter_t), false);
        fprintf('%d pulses at (%d, %d, %d) um\n', tmp_num_pulse, ...
            round(wbim.sample_xyz_um));
        tmp_output_task.start();
        % Does DIO task has a waitUntilIdle function? Use pause at the
        % moment 
        pause(1);
    end
    wbim.ablation.si_cleanup_pulse_task(tmp_output_task);
end
%% Ablation density
tmp_num_pulse = 1;
tmp_output_task = wbim.ablation.setup_ablation_pulses_task(tmp_num_pulse);
trial_spacing_um = 50;
num_p = 10;
num_trial = 10;
x_um_vec = (0 : 1 : (num_p - 1)) * trial_spacing_um;
x_um_vec = x_um_vec - x_um_vec(ceil(num_trial/2));
y_um_vec = (0 : 1 : (num_trial - 1)) * trial_spacing_um;
y_um_vec = y_um_vec - y_um_vec(ceil(num_trial/2));
init_pos_um = wbim.sample_xyz_um;
x_um_vec = x_um_vec + init_pos_um(1);
y_um_vec = y_um_vec + init_pos_um(2);
for iter_np = 1 : num_p
    % Generate pulse task     
    wbim.move_sample_along_axis_um(1, x_um_vec(iter_np), false);
    for iter_t = 1 : num_trial
        wbim.move_sample_along_axis_um(2, y_um_vec(iter_t), false);
        fprintf('%d pulses at (%d, %d, %d) um\n', tmp_num_pulse, ...
            round(wbim.sample_xyz_um));
        tmp_output_task.start();
        % Does DIO task has a waitUntilIdle function? Use pause at the
        % moment 
    end
end
wbim.ablation.si_cleanup_pulse_task(tmp_output_task);
%% Swtich back to imaging mode
wbim.ablation.close_ablation_shutter();
wbim.ablation.h_zaber_controller.set_ablation_do_value(0);
wbim.ablation.actuator_switch_to_dichroic();
wbim.ablation.turn_off_pump();
wbim.imaging.move_sample_along_axis_um(1:3, init_pos_um(1:3))
%% Image ablation site
% Take stacks in scan mode 
[exp_x_um, exp_y_um] = meshgrid(x_um_vec, y_um_vec);
exp_x_um = exp_x_um(:);
exp_y_um = exp_y_um(:);
site_xy_um = cat(2, exp_x_um, exp_y_um);

acq_mode = WBIMMicroscopeMode.Scan;
wbim.imaging.set_acq_mode(acq_mode);
wbim.imaging.h_tile_manager.add_tile_at_xy(site_xy_um, true);
%%
continuous_acq_Q = true;
wbim.imaging.init_scanning(continuous_acq_Q);
wbim.imaging.scan_next_tile();
