classdef WBIMWaiter < handle
% When `pause` is execuated, MATLAB yield executation to the event queue.
% The asynchronous processes controlled by the timer can be run during the
% pause period. 
% Notice that this is not always the case for other command. For a
% computationally intensive command, the event might not be execuated until
% the computationally intensive command is finished. However, MATLAB will
% check the event queue after a single iteration of for / while loop is
% done. 
% I am not quite sure whether the extra waiting time is necessary. It was
% necessary previously to ensure the correct order of multiple nested event
% callbacks. 
% Implemented by Xiang Ji on 04/12/2024
    
    properties
        check_period_s = 1.0;
    end
    
    methods(Static)
        function wait_for_imaging_done(wbim, opt)
            arguments
                wbim WBIMControl
                opt.check_period_s = 1.0;
                opt.wait_time_after_finish_s = 1.0;
            end
            running_imaging_Q = true;
            t_tic = tic;
            while running_imaging_Q
                if strcmpi(wbim.hSI.acqState, 'idle')
                    running_imaging_Q = false;
                    pause(opt.wait_time_after_finish_s);
                else
                    pause(opt.check_period_s);
                end
%                 fprintf("Waiting for imaging. Elapsed time is %.2f seconds\n", toc(t_tic));
%                 wbim.h_logger.write("DEBUG", ...
%                     sprintf("Waiting for imaging. Elapsed time is %.2f seconds", ...
%                     toc(t_tic)));
            end            
        end    
        
        function wait_for_ablation_done(wbim, opt)
            arguments
                wbim WBIMControl
                opt.check_period_s = 1.0;
                opt.wait_time_after_finish_s = 1.0;
            end
            running_ablation_Q = true;
            t_tic = tic;
            while running_ablation_Q
                if wbim.ablation.current_state == WBIMMachineState.Idle
                    running_ablation_Q = false;
                    pause(opt.wait_time_after_finish_s);
                else
                    pause(opt.check_period_s);
                end
%                 fprintf("Waiting for ablation. Elapsed time is %.2f seconds\n", toc(t_tic));
%                 wbim.h_logger.write("DEBUG", ...
%                     sprintf("Waiting for ablation. Elapsed time is %.2f seconds", ...
%                     toc(t_tic)));
            end
        end
    end
    
end