classdef WBIMAcqPriorityQueue < handle
    
    properties(SetAccess=protected)
        queue
        size (1,1) double        
    end
    
    properties(Access=private)
        
    end
%%
    methods
        function obj = WBIMAcqPriorityQueue()
            obj.clear();
        end
        
        function append(obj, val)
           if iscolun(val)
               val = val .';
           end
           obj.queue = cat(2, obj.queue, val);            
        end
        
        function val = get.size(obj)
           val = numel(obj.queue); 
        end
        
        function val = isEmpty(obj)
           val = (obj.size == 0); 
        end
        
        function inQ = contains(obj, vec)
            inQ = ismember(vec, obj.queue);
        end
        
        function val = getFirst(obj)
           val = obj.queue(1); 
        end
        
        function removeFirst(obj)
            if obj.size >= 1
               obj.queue(1) = []; 
            end
        end
        
        function clear(obj)
            obj.queue = [];
            obj.size = 0;
        end
    end    
end