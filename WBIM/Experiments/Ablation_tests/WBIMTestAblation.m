classdef WBIMTestAblation < handle
    
properties
    
    
end


methods(Static)
    function single_plane_site_ablation(wbim, opt)
        arguments
            wbim WBIMControl
            opt.fluence_J_cm2 {mustBeNonnegative}
            opt.num_pulse_per_site {mustBeNonnegative}
            opt.site_sample_pos_x_um
            opt.site_sample_pos_y_um
        end
        
        
        
        
        array_size = size(opt.site_sample_pos_x_um);
        assert(all(array_size == size(opt.site_sample_pos_y_um)));
        
        wbim.h_logger.write("MESSAGE", "Start single site ablation test");
        wbim.switch_to_ablation_mode();
        wbim.ablation.h_zaber_controller.set_ablation_do_value(1);
        wbim.ablation.open_ablation_shutter();
        init_pos_um = wbim.sample_xyz_um;                
        for iter_2 = 1 : array_size(2)
            for iter_1 = 1 : array_size(1)
                % Generate pulse task
                tmp_num_pulse = opt.num_pulse_per_site(iter_1, iter_2);
                tmp_fluence_J_cm2 = opt.fluence_J_cm2(iter_1, iter_2);
                
                tmp_output_task = wbim.ablation.setup_ablation_pulses_task(tmp_num_pulse);
                wbim.h_logger.write("MESSAGE", sprintf('Target fluence %.2f J/cm2.',... 
                    tmp_fluence_J_cm2));
                wbim.ablation.set_fluence_J_per_cm2(tmp_fluence_J_cm2 / ND_decay_ratio);
                
                wbim.move_sample_along_axis_um(1:2, ...
                    [opt.site_sample_pos_x_um(iter_1, iter_2), opt.site_sample_pos_y_um(iter_1, iter_2)], ...
                    false);
                wbim.h_logger.write("MESSAGE", sprintf('%d pulses at (%d, %d, %d) um', tmp_num_pulse, ...
                    round(wbim.sample_xyz_um)));
                pause(1);
                tmp_output_task.start();
                wbim.ablation.si_cleanup_pulse_task(tmp_output_task);
                pause(1);
            end            
        end        
        wbim.ablation.close_ablation_shutter();
        wbim.ablation.h_zaber_controller.set_ablation_do_value(0);
        wbim.move_sample_along_axis_um(1:3, init_pos_um);
        fprintf('Finish ablation task');
        wbim.switch_to_imaging_mode();
    end    
    
    
    
end

end