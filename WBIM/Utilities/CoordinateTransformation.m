classdef CoordinateTransformation < handle
    
    properties
        A (4,4) double
    end
    
    properties(SetAccess=protected)
        R (3,3) double
        t (3, 1) double
        R_xy (2,2) double
        t_xy (2,1) double
        A_xy (3,3) double
        A_xy_local (3,3) double
        A_inv (4,4) double
    end
    
    methods
        function obj = CoordinateTransformation(A)
            obj.A = A;
        end
        
        function v_hstack = forward_transform(obj, v_hstack, is_H_Q)
            if nargin < 2
                is_H_Q = true;
            end
            if is_H_Q
                v_hstack = obj.A * v_hstack;
            else
                v_hstack = obj.R * v_hstack + obj.t;
            end
        end
        
        function v_hstack_H = inverse_transform(obj, v_hstack_H, is_H_Q)
            if nargin < 2
                is_H_Q = true;
            end
            if is_H_Q
                v_hstack_H = obj.A_inv * v_hstack_H;
            else
                v_hstack_H = obj.R' * v_hstack_H - obj.R' *  obj.t;
            end
        end
    end
    
    %%
    methods
        function set.A(obj, val)
            obj.A = val;
            obj.update_by_A();
        end
    end
    
    methods(Hidden)
        function update_by_A(obj)
            obj.R = obj.A(1:3, 1:3);
            obj.t = obj.A(1:3, 4);
            obj.R_xy = obj.A(1:2, 1:2);
            obj.t_xy = obj.A(1:2, 4);
            obj.A_xy = obj.A([1,2,4], [1,2,4]);
            obj.A_xy_local = obj.A_xy;
            obj.A_xy_local([1,2], 3) = 0;
            obj.A_inv = cat(1, cat(2, obj.R', -obj.R' * obj.t), [0,0,0,1]);
        end
    end
    %%
    methods(Static)
        function v_hstack = to_homogeneous_coordiante(v_hstack)
            num_v = size(v_hstack, 2);
            v_hstack = cat(1, v_hstack, ones(1, num_v));
        end
        
        function v_hstack = to_inhomogeneous_coordiante(v_hstack)
            v_hstack = v_hstack(1:(end-1), :) ./ v_hstack(end, :);
        end
        
    end
end