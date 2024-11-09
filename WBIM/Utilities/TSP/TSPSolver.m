classdef TSPSolver < handle    
    
    methods(Static)
        function order_ind = get_node_order_from_distance_matrix_2opt(...
                dist_mat, options)
            % Solve the traveling salesman problem using 2opt algorithm. 
            % Input: 
            %   dist_mat: squared matrix; dist_mat(i, j) is the distance
            %   from node i to node j. This matrix should be squared, but
            %   does not need to be symmetric. 
            arguments
                dist_mat {mustBeNumeric}
                options.max_iteration (1,1) {mustBeNumeric} = 50;
            end
            mat_size = size(dist_mat);
            assert(mat_size(1) == mat_size(2), 'dist_mat should be a square matrix');
            
            num_pts = mat_size(1);
            order_ind = (1 : num_pts).';
            if num_pts <= 2
                return;
            end
            i_count = 1;
            found_improvement = true;
            
            while found_improvement && (i_count < options.max_iteration)
                found_improvement = false;
                for i = 1 : (num_pts - 1)
                    for j = (i+1) : (num_pts - 1)
                        im1 = mod(i - 2, num_pts) + 1;
                        jp1 = mod(j, num_pts) + 1;
                        length_delta = dist_mat(order_ind(im1), order_ind(i)) + ...
                            dist_mat(order_ind(j), order_ind(jp1)) - ...
                            dist_mat(order_ind(im1), order_ind(j)) - ...
                            dist_mat(order_ind(i), order_ind(jp1));
                        if length_delta > eps
                            found_improvement = true;
                            order_ind(i:j) = order_ind(j:-1:i);
                        end
                    end
                end
                i_count = i_count + 1;
            end
        end
        
        function order_idx = get_ordered_position(dist_mat, options)
            arguments
                dist_mat {mustBeNumeric}
                options.NoReturnQ (1,1) logical = true;
                options.FixStartQ (1,1) logical = false;
            end
            num_pos = size(dist_mat, 1);
            if options.NoReturnQ
                dist_mat = cat(1, cat(2, dist_mat, zeros(num_pos, 1)), ...
                    zeros(1, num_pos + 1));
            end
            
            if options.FixStartQ
                dist_mat = cat(1, cat(2, dist_mat, inf(num_pos + 1, 1)), inf(1, num_pos + 2));
                dist_mat(num_pos + 2, 1) = 0;
                dist_mat(num_pos + 1, num_pos + 2) = 0;
            end
            order_idx = TSPSolver.get_node_order_from_distance_matrix_2opt(...
                dist_mat);
            if options.NoReturnQ
                order_idx = order_idx(order_idx <= num_pos);
            end
            if options.FixStartQ
                assert(order_idx(1) == 1, 'The first node index should be 1');
                % Remove the first one
%                 order_idx = order_idx(2:end) - 1;
            end
        end
    end
    
end