classdef WBIMDelayTimer < handle
    % Execuate the function handle once with a specific delay (in second)
    properties
        timer (1,1) timer
        fun_hdl
        fun_args
    end

    methods
        function obj = WBIMDelayTimer(delay_s, fun_hdl, fun_arg)
            arguments
                delay_s (1, :) double = []
                fun_hdl = []
                fun_arg = {}
            end
            try
                obj.timer = timer('BusyMode', 'queue', 'ExecutionMode','singleShot');
                obj.timer.TimerFcn = {@obj.TimerFcn, obj};
                obj.timer.StopFcn = {@obj.StopFcn, obj};
            catch ME
                delete(obj.timer);
                rethrow(ME);
            end
            if ~isempty(delay_s)
                obj.run_with_delay(delay_s, fun_hdl, fun_arg);
            end
        end

        function run_with_delay(obj, delay_s, fun_hdl, fun_arg)
            arguments
                obj (1,1) WBIMDelayTimer
                delay_s (1,1) double {mustBeNonnegative}
                fun_hdl
                fun_arg = {}
            end
            obj.fun_hdl = fun_hdl;
            obj.fun_args = fun_arg;
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
                obj WBIMDelayTimer
            end
            
            if ~isempty(obj.fun_hdl)
                if isempty(obj.fun_args)
                    out = obj.fun_hdl();
                else
                    out = obj.fun_hdl(obj.fun_args{:});
                end
            end
        end

        function out = StopFcn(t_obj, evnt, obj)
            obj.delete();
        end
    end
end