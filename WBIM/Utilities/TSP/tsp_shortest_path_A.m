function order_ind = tsp_shortest_path_A(dist_mat)    
    verboseQ = true;
    mat_size = size(dist_mat);
    num_pts = mat_size(1);    
    order_ind = (1 : num_pts).';
    if num_pts <= 2
        return;
    end
    i_count = 1;
    max_iteration = 1e4;
    no_change_count = 0;
    while i_count < max_iteration && no_change_count < num_pts
        flip_length = mod(i_count - 1, num_pts - 3) + 2;
        adj_ind = sub2ind(mat_size, order_ind, circshift(order_ind, -1));
        adj_dist = dist_mat(adj_ind);  
        step_ind = sub2ind(mat_size, order_ind, circshift(order_ind, -flip_length));
        step_dist = dist_mat(step_ind);
        
        adj_dist = cat(1, 0, adj_dist(1:end-1));
        step_dist = cat(1, 0, step_dist(1 : end-flip_length));        
        
        % Cost before flipping: 
        cost_0 = adj_dist(1 : end - flip_length) + adj_dist(flip_length+1 : end);
        % Cost after flipping: 
        cost_1 = step_dist(1:end-1) + step_dist(2:end);
        
%         cost_0 = adj_dist + circshift(adj_dist, -flip_length);
%         cost_1 = step_dist + circshift(step_dist, -1);
        cost_reduction = cost_1 - cost_0;
        % Flip the paths that maximum the cost reduction among all the
        % locally beneficially adjustments 
        
        % The moving minimum window size will effect the number of
        % iterations and the final path length, but its effect on the final
        % quality and performance is unclear. 
        lm_window_size = flip_length * 2 + 1;
%         lm_window_half_size = ceil(flip_length/2);
%         cost_pad = cat(1, cost_reduction(end-lm_window_half_size+1:end), cost_reduction, ...
%             cost_reduction(1 : lm_window_half_size));
        cost_pad = cost_reduction;
        mov_min_reduction = movmin(cost_pad, lm_window_size);
%         mov_min_reduction = mov_min_reduction(lm_window_half_size + 1: end - lm_window_half_size);        
        min_ind = find(mov_min_reduction == cost_reduction & cost_reduction < -eps);
        % Do not adjust the first point
%         min_ind = min_ind(min_ind ~= 1);
        if ~isempty(min_ind)
            for i = 1 : numel(min_ind)
                sequence_ind = min_ind(i) : (min_ind(i) + flip_length - 1);
                sequence_ind = mod(sequence_ind - 1, num_pts) + 1;
                order_ind(sequence_ind) = order_ind(sequence_ind(end:-1:1));
            end
            no_change_count = 0;
        else
            no_change_count = no_change_count + 1;
        end
        i_count = i_count + 1;        
        if verboseQ && mod(i_count, 50) == 0
            fprintf('Iteration %d. Total path length: %.2f\n', i_count, ...
                sum(adj_dist));
        end
    end     
    if verboseQ
       fprintf('Final iteration %d. Total path length: %.2f\n', i_count, ...
           sum(adj_dist));
    end
end