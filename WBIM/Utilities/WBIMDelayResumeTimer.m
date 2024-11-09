classdef WBIMDelayResumeTimer < handle
    % Execuate the function handle once with a specific delay (in second)
    properties
        timer (1,1) timer
        im_hdl % WBIMImaging handle
    end

    methods
        function obj = WBIMDelayResumeTimer(delay_s, im_hdl)
            arguments
                delay_s (1, 1) double
                im_hdl 
            end
            obj.im_hdl = im_hdl;
            try
                obj.timer = timer('BusyMode', 'queue', 'ExecutionMode','singleShot');
                obj.timer.Name = 'WBIMDelayResumeTimer';
                obj.timer.TimerFcn = {@obj.TimerFcn, obj};
                obj.timer.StopFcn = {@obj.StopFcn, obj};
            catch ME
                delete(obj.timer);
                rethrow(ME);
            end
            if strcmpi(obj.timer.Running, 'on')
                obj.timer.stop();
            end     
            obj.timer.StartDelay = delay_s;
            obj.timer.start();
        end

        function delete(obj)
            delete(obj.timer);
        end
    end

    methods(Static)
        function TimerFcn(t_obj, evnt, obj)
            arguments
                t_obj timer
                evnt
                obj WBIMDelayResumeTimer
            end
            obj.im_hdl.resume_scanning(false);
        end

        function out = StopFcn(t_obj, evnt, obj)
            
        end
    end
end