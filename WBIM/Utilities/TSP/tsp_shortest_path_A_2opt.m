function order_ind = tsp_shortest_path_A_2opt(dist_mat, verboseQ)
% dist_mat(i,j) is the distance from i to j
% This distance matrix could be asymmetric - when introducing the dummy
% node
% The order of nodes in the distance matrix is the initial order for
% searching the path 
if nargin < 2
    verboseQ = false;
end
mat_size = size(dist_mat);
assert(mat_size(1) == mat_size(2), 'dist_mat should be a square matrix');

num_pts = mat_size(1);
order_ind = (1 : num_pts).';
if num_pts <= 2
    return;
end
i_count = 1;
max_iteration = 50;
found_improvement = true;

while found_improvement && (i_count < max_iteration)
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
if verboseQ
    adj_ind = sub2ind(mat_size, order_ind, circshift(order_ind, -1));
    fprintf('Final iteration %d. Total path length: %.3f.\n', i_count, ...
        sum(dist_mat(adj_ind)));
end
end