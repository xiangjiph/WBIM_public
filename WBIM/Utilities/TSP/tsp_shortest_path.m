function order_ind = tsp_shortest_path(xy_list)
    dist_name = 'euclidean';
    verboseQ = true;
    visQ = true;

    assert(size(xy_list, 2) == 2, 'xy_list should be a N-by-2 matrix, where each row is the coordinate of a point');
    num_pts = size(xy_list, 1);
    order_ind = (1 : num_pts).';
    if num_pts <= 2
        return;
    end
    i_count = 1;
    max_iteration = 1e4;
    no_change_count = 0;
    while i_count < max_iteration && no_change_count < num_pts
        flip_length = mod(i_count - 1, num_pts - 3) + 2;
        
        adj_dist = step_pair_dist(xy_list(order_ind, :), 1, dist_name);
        step_dist = step_pair_dist(xy_list(order_ind, :), flip_length, dist_name);
        % Cost before flipping: 
        cost_0 = adj_dist(1 : end - flip_length) + ...
            adj_dist(flip_length + 1 : end);
        % Cost after flipping: 
        cost_1 = step_dist(1 : end - 1) + step_dist(2 : end);
        cost_reduction = cost_1 - cost_0;
        % Flip the paths that maximum the cost reduction among all the
        % locally beneficially adjustments 
        
        % The moving minimum window size will effect the number of
        % iterations and the final path length, but its effect on the final
        % quality and performance is unclear. 
        
%         mov_min_reduction = movmin(cost_reduction, flip_length);
        mov_min_reduction = movmin(cost_reduction, 2 * flip_length + 1);
        min_ind = find(mov_min_reduction == cost_reduction & cost_reduction < -eps);
        % Do not adjust the first point
%         min_ind = min_ind(min_ind ~= 1);
        if ~isempty(min_ind)
            for i = 1 : numel(min_ind)
                order_ind(min_ind(i) : min_ind(i) + flip_length - 1) = ...
                    order_ind(min_ind(i) + flip_length-1 : -1 : min_ind(i));
            end
            no_change_count = 0;
        else
            no_change_count = no_change_count + 1;
        end
        i_count = i_count + 1;        
    end     
    if verboseQ
       fprintf('Final iteration %d. Total path length: %.2f\n', i_count, ...
           sum(adj_dist));
    end
    if visQ
       fig_hdl = figure;
       ax_hdl = axes(fig_hdl);
       plot(ax_hdl, xy_list(order_ind, 1), xy_list(order_ind, 2), '-ro');
       hold(ax_hdl, 'on');
       scatter(ax_hdl, xy_list(1, 1), xy_list(1, 2), 'g*');
    end
end

function d = step_pair_dist(xy_list, step, dist_name)
    if nargin < 3
        dist_name = 'euclidean';
    end
    dif = xy_list(1 : end-step, :) - xy_list(step+1 : end, :);
    switch dist_name
        case 'euclidean'
            d = cat(1, 0, sqrt(sum(dif .^2, 2)));
    end
end