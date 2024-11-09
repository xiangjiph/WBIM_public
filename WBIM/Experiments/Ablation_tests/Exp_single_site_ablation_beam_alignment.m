%%
wbim.ablation.set_diffuser_state(false);
wbim.ablation.set_fluence_J_per_cm2(16);
% wbim.ablation.set_ablation_fractional_power(0.01);
num_pulse_per_site = 5;
tmp_output_task = wbim.ablation.setup_ablation_pulses_task(num_pulse_per_site);
%% Single site ablation
wbim.ablation.turn_on_pump();
wbim.ablation.actuator_switch_to_mirror();
wbim.ablation.h_zaber_controller.set_ablation_do_value(1);
wbim.ablation.open_ablation_shutter();
pause(2);
tmp_output_task.start();
t_tic = tic;
wbim.ablation.close_ablation_shutter();
wbim.ablation.turn_off_pump();
wbim.ablation.h_zaber_controller.set_ablation_do_value(0);
wbim.ablation.actuator_switch_to_dichroic();
hSI.startFocus
% hSI.startGrab;
fprintf("Time between ablation and start acquisition: %.2f seconds\n", ...
    toc(t_tic));
%%
wbim.ablation.si_cleanup_pulse_task(tmp_output_task);
%%
wbim.ablation.turn_on_pump();
wbim.ablation.turn_off_pump();
