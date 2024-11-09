% Run this after initializing wbim 
myRI = RICorrection();
% myRI.hValve.v_5s_mL = 8.25; %mL % What's this?
% myRI.init_logger(wbim.fp_RI_log);
% myRI.total_volume_mL = 1500;
%%
rs = dabs.resources.ResourceStore();
h_pump_port = rs.filterByName(WBIMConfig.ABLATION_PUMP_PORT_CHANNEL_NAME);
h_pump_port.setValue(0);
%%
% myRI.total_volume_mL = 1500;
% myRI.intercorrectionTime = 300;
% % pause(5000)
% myRI.startCorrection
%%
myRI.hRefractometer.port_DEV.LCD_ON.setValue(0);
for i = 1 : 100
    lcd_state = myRI.hRefractometer.port_DEV.LCD_STATE.readValue();
    fprintf('%s Current LDC state: %d\n', datestr(now), lcd_state);
    pause(1);
end
%%
rt_lcd_on_port = myRI.hRefractometer.port_DEV.LCD_ON;

pause_duration_s = 0.1;
h_si_task = dabs.vidrio.ddi.DoTask(rt_lcd_on_port.hDAQ, 'pulse trains');
h_si_task.addChannel(rt_lcd_on_port.name);
h_si_task.sampleMode = 'finite';
num_high = pause_duration_s * h_si_task.sampleRate;
buffer = cat(1, ones(num_high, 1), 0);
h_si_task.writeOutputBuffer(buffer);
h_si_task.samplesPerTrigger = numel(buffer);
h_si_task.start();
for i = 1 : 10
    pause(0.1);
    fprintf('%s: Current LCD state: %d\n', datestr(now, 'HH:MM:SS.FFF'),  myRI.hRefractometer.port_DEV.LCD_STATE.readValue());
end
WBIMAblation.si_cleanup_pulse_task(h_si_task);
%%
wbim.ablation.init_refractometer();