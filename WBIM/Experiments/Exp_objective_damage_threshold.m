%% Reduce the ablation speed
% ablation_ctrl.actuator_switch_to_dichroic();
wbim.ablation.h_hwp.home();
wbim.ablation.set_ablation_fractional_power(1);
%%
wbim.ablation.actuator_switch_to_mirror();
wbim.ablation.open_ablation_shutter();
abl_trigger_task = wbim.ablation.setup_ablation_pulses_task(1000 * 60 * 15);
zaber_port_id = wbim.ablation.h_zaber_controller.ablation_do_port;
%%
wbim.ablation.h_zaber_controller.set_port_output('digital', zaber_port_id, 1);
abl_trigger_task.start();
fprintf('Ablation test start: %s\n', datestr(now));
%% Clean up
wbim.ablation.h_zaber_controller.set_port_output('digital', zaber_port_id, 0);
wbim.ablation.close_ablation_shutter();
%%
abl_trigger_task.stop();
abl_trigger_task.abort();
abl_trigger_task.delete();
clear abl_trigger_task

%% Note: 
% 1. 20% during ROI ablation - damaged (second line) 
% 2. 4% for 1 minute - appear to be safe
% 3. 4% for 5 minute - appear to be safe (inspected under microscope)
% 4. 10% for 5 minute - visible damage (by eye)
% 5. 7% for 5 mintue - appear to be safe
% 6. 7% for 15 minute - visible damaged
%%
wbim.ablation.set_ablation_fractional_power(0.25);
wbim.ablation.open_ablation_shutter();
fprintf('Ablation test starts at %s\n', datestr(now));
pause(60 * 15);
wbim.ablation.close_ablation_shutter();
fprintf('Ablation test ends at %s\n', datestr(now));
%%
