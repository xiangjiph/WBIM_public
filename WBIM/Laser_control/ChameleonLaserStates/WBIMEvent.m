classdef (ConstructOnLoad) WBIMEvent < event.EventData
    properties
        UserData
    end

    methods
        function obj = WBIMEvent(val)
            obj.UserData = val;
        end
    end
end