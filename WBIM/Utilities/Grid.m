classdef Grid < handle & matlab.mixin.Copyable
    %% To do list:
    
    
    %%
    properties
        % All in pixel value
        space_size
        space_dim
        
        tile_size
        tile_overlap
        tile_spacing
        
        num_tile
        
        grid_size
        
        ul
        ll
        ll_array
        
        mmxx
        mmxx_array
        
        mmll
        mmll_array
        
        center_sub
        array_idx_sub
    end
    
    
    methods
        %% Constructor
        function obj = Grid(space_size, tile_size, tile_overlap)
            obj.update_grid_info(space_size, tile_size, tile_overlap);
        end
        
        function obj = update_grid_info(obj, space_size, tile_size, tile_overlap)
            obj.space_size = space_size;
            obj.space_dim = numel(space_size);
            
            if isscalar(tile_size)
                tile_size = tile_size .* ones(1, obj.space_dim);
            end
            obj.tile_size = tile_size;
            if isscalar(tile_overlap)
                tile_overlap = tile_overlap .* ones(1, obj.space_dim);
            end
            obj.tile_overlap = tile_overlap;
            obj.tile_spacing = tile_size - tile_overlap;
            % Upper left points for each block
            upper_left_pos = cell(obj.space_dim, 1);
            for idim = 1 : obj.space_dim
                upper_left_pos{idim} = 1 : (tile_size(idim) - tile_overlap(idim)) : space_size(idim);
            end
            switch obj.space_dim
                case 3
                    [grid_ul_1, grid_ul_2, grid_ul_3] = ndgrid(upper_left_pos{1}, upper_left_pos{2}, upper_left_pos{3});
                    grid_ul_pos = [grid_ul_1(:) grid_ul_2(:) grid_ul_3(:)];
                case 2
                    [grid_ul_1, grid_ul_2] = ndgrid(upper_left_pos{1}, upper_left_pos{2});
                    grid_ul_pos = [grid_ul_1(:) grid_ul_2(:)];
                case 1
                    [grid_ul_1] = ndgrid(upper_left_pos{1});
                    grid_ul_pos = [grid_ul_1(:)];
            end
            
            obj.ul = grid_ul_pos;
            % [sub1_min, sub2_min, sub1_max, sub2_max]
            obj.mmxx = [grid_ul_pos, min(grid_ul_pos + tile_size - 1, space_size)];
            % [row, col, l1, l2]
            obj.mmll = [grid_ul_pos, obj.mmxx(:,(obj.space_dim+1):(2*obj.space_dim)) - grid_ul_pos + 1 ];
            obj.center_sub = round((obj.mmxx(:,1:obj.space_dim) + obj.mmxx(:, (obj.space_dim + 1): end))/2);
            obj.array_idx_sub =  bsxfun(@rdivide, (grid_ul_pos - 1), (tile_size - tile_overlap)) + 1;
            obj.ll = [obj.mmxx(:,(obj.space_dim+1):(2*obj.space_dim)) - grid_ul_pos + 1];
            obj.grid_size = cellfun(@length, upper_left_pos);
            obj.num_tile = prod(obj.grid_size);
            obj.grid_size = max(obj.array_idx_sub, [], 1);
            
            obj.mmll_array = reshape(permute(obj.mmll, [2,1]), [size(obj.mmll, 2), obj.grid_size]);
            obj.mmxx_array = reshape(permute(obj.mmxx, [2,1]), [size(obj.mmxx, 2), obj.grid_size]);
            obj.ll_array = reshape(permute(obj.ll, [2,1]), [size(obj.ll, 2), obj.grid_size]);
        end
        %% Coordiante transformation
        function grid_sub = pos_to_grid_sub(obj, pos)
            assert(size(pos, 2) == obj.space_dim, 'pos should be a matrix (num_points x dim)');
            grid_sub = ceil(pos ./ obj.tile_spacing);
        end
        
        function grid_ind = pos_to_grid_ind(obj, pos_pxl)
            grid_sub = obj.pos_to_grid_sub(pos_pxl);
            grid_ind = obj.grid_sub2ind(grid_sub);
        end
        
        function grid_sub = grid_ind2sub(obj, grid_ind)
            grid_sub = obj.ind2sub(obj.grid_size, grid_ind);
        end
        
        function grid_ind = grid_sub2ind(obj, grid_sub)
            grid_ind = obj.sub2ind(obj.grid_size, grid_sub);
        end
        %% ROI occupancy grid
        function og = initialize_occupancy_grid(obj, pxl_sub_yx)
            if nargin < 2
                pxl_sub_yx = [];
            end
            og = false(obj.grid_size);
            og = obj.activate_occupancy_grid(og, pxl_sub_yx);
        end
        
        function og = activate_occupancy_grid(obj, og, pxl_sub)
            if ~isempty(pxl_sub)
                grid_ind = obj.pos_to_grid_ind(pxl_sub);
                og(grid_ind) = true;
            end
        end
        
        function tile_info = get_tile_info(obj, grid_ind)
            tile_info = struct;
            tile_info.grid_ind = grid_ind;
            tile_info.tile_size = obj.tile_size;
            tile_info.center_xy = obj.center_sub(grid_ind, [2,1]);
            tile_info.bbox_mmxx = obj.mmxx(grid_ind, :);
        end
        %% ROI acquisition grid
        function ag = initialize_acquisition_grid(obj)
            ag = nan(obj.grid_size);
        end
        
        function [tile_info, cc_info] = analyze_acquisition_grid(obj, ag, axis_order)
            if nargin < 3
                axis_order = [];
            end
            og = obj.convert_acquisition_grid_to_occupancy_grid(ag);
            [tile_info, cc_info] = obj.analyze_occupancy_grid(og, axis_order);
        end
        %% Utilities
        function grid_ind = get_tiles_with_center_in_bbox(obj, bbox_yx_mmxx)
            assert(numel(bbox_yx_mmxx) == 4 &&...
                all(bbox_yx_mmxx(1:2) <= bbox_yx_mmxx(3:4)));
            inQ = obj.center_sub(:, 1) >= bbox_yx_mmxx(1) & ...
                obj.center_sub(:, 2) >= bbox_yx_mmxx(2) & ...
                obj.center_sub(:, 1) <= bbox_yx_mmxx(3) & ...
                obj.center_sub(:, 2) <= bbox_yx_mmxx(4);
            grid_ind = find(inQ);
        end
        
        function grid_ind = get_tiles_completely_in_bbox(obj, bbox_yx_mmxx)
            arguments
                obj (1,1) WBIMImagingGrid2D
                bbox_yx_mmxx (1, 4) {mustBePositive}
            end
            inQ = obj.mmxx(:, 1) >= bbox_yx_mmxx(1) & ...
                obj.mmxx(:, 2) >= bbox_yx_mmxx(2) & ...
                obj.mmxx(:, 3) <= bbox_yx_mmxx(3) & ...
                obj.mmxx(:, 4) <= bbox_yx_mmxx(4);
            grid_ind = find(inQ);            
        end        
        
        function grid_ind = get_tiles_overlap_with_bbox(obj, bbox_yx_mmxx)
            arguments
                obj (1,1) WBIMImagingGrid2D
                bbox_yx_mmxx (1, 4) {mustBePositive}
            end
            inQ = Grid.check_bounding_box_overlap(bbox_yx_mmxx, obj.mmxx);
            grid_ind = find(inQ);
        end        
    end
    
    methods(Static)
        % overload
        function sub = ind2sub(block_size, ind)
            validateattributes(ind, {'numeric'}, {'positive', 'integer'});
            validateattributes(block_size, {'numeric'}, {'positive', 'integer'});
            block_dim = nnz(block_size ~= 1);
            if ind > prod(block_size)
                error('ind is out of range');
            end
            switch block_dim
                case 1
                    sub = ind;
                case 2
                    sub = zeros(numel(ind), 2);
                    [sub(:,1), sub(:,2)] = ind2sub(block_size, ind);
                case 3
                    sub = zeros(numel(ind), 3);
                    [sub(:,1), sub(:,2), sub(:,3)] = ind2sub(block_size, ind);
                otherwise
                    error('Unsupported block size');
            end
        end
        
        function ind = sub2ind(block_size, sub)
            %             space_dim = nnz(block_size ~= 1);
            space_dim = numel(block_size);
            assert(size(sub, 2) == space_dim);
            switch space_dim
                case 1
                    ind = sub;
                case 2
                    ind = sub2ind(block_size, sub(:, 1), sub(:, 2));
                case 3
                    ind = sub2ind(block_size, sub(:, 1), sub(:, 2), sub(:, 3));
                otherwise
                    error('To be implemented');
            end
        end
        
        
        %% Utility
        function og = expand_occupancy_grid(og, margin)
            % Uniformally expand the roi using image dilation
            space_dim = ndims(og);
            switch space_dim
                case 2
                    if isvector(og)
                        og = imdilate(og, strel('line', margin, 0));
                    else
                        og = imdilate(og, strel('disk', margin));
                    end
                case 3
                    og = imdilate(og, strel('sphere', margin));
                otherwise
                    error('To be implemented');
            end
        end
        
        function [bin_cell_array, varargout] = bin_data_to_idx_list(data)
            % fun_bin_data_to_idx_list bin the data according to their values and
            % output the corresponding index list
            % Input:
            %   data: numerical vector
            % Output:
            %   bin_cell_array: cell array, each cell constains a vector, whose
            %   components are the indices of the component of data that have the same
            %   value.
            %   varargout: unique data value
            
            if isempty(data)
                bin_cell_array = {};
                varargout{1} = [];
                return;
            end
            num_data = numel(data);
            if ~issorted(data)
                [data, idx_list] = sort(data(:), 'ascend');
                % To be consistent with the previous behavior
                idx_list = idx_list.';
            else
                idx_list = 1 : num_data;
            end
            est_num_bin = max(2, round(data(end) - data(1)) + 1);
            
            if num_data < 65535
                idx_list = uint16(idx_list);
            elseif num_data < 4294967295
                idx_list = uint32(idx_list);
            end
            bin_value_list = zeros(est_num_bin,1);
            bin_cell_array = cell(est_num_bin,1);
            left_idx = 1;
            bin_data = data(1);
            num_bin = 0;
            for right_idx = 1 : num_data
                tmp_data = data(right_idx);
                if tmp_data ~= bin_data
                    num_bin = num_bin + 1;
                    % copy
                    bin_cell_array{num_bin} = idx_list(left_idx : (right_idx - 1));
                    bin_value_list(num_bin) = bin_data;
                    % Update
                    bin_data = tmp_data;
                    % Re-initialize
                    left_idx = right_idx;
                end
            end
            num_bin = num_bin + 1;
            bin_cell_array{num_bin} = idx_list(left_idx : right_idx);
            bin_cell_array(num_bin + 1 : end) = [];
            bin_value_list(num_bin) = bin_data;
            bin_value_list = bin_value_list(1 : num_bin);
            if nargout > 1
                varargout{1} = bin_value_list;
            end
        end
        
        function binned_data = bin_data_to_cells_by_ind(data, bin_ind)
            % fun_bin_data_to_cells_by_ind reorganize the data into cell arraies
            % accoriding to the bin_ind
            % Input:
            %   data: numerical vector
            %   bin_ind: cell array, indices of the data that need to be put inside
            %   each cell
            % Output:
            %   binned_data: cell array
            %
            num_cell = numel(bin_ind);
            binned_data = cell(num_cell, 1);
            for iter_cell = 1 : num_cell
                binned_data{iter_cell} = data(bin_ind{iter_cell});
            end
        end
        
        function og = convert_acquisition_grid_to_occupancy_grid(ag)
            og = ~isnan(ag);
        end
        
        function [local_im, varargout] = stitch_tiles(bbox_mmxx_pxl, tiles_value, stitch_method)
            if nargin < 3
                stitch_method = 'overwrite';
            end
            [num_tiles, num_dim] = size(bbox_mmxx_pxl);
            num_dim = num_dim / 2;
            assert(num_tiles == numel(tiles_value), 'One bounding box for one entry in tiles_value');
            
            bbox_mm = min(bbox_mmxx_pxl(:, 1:num_dim), [], 1);
            bbox_xx = max(bbox_mmxx_pxl(:, num_dim+1:end), [], 1);
            bbox_ll = bbox_xx - bbox_mm + 1;
            
            local_bbox_mmxx = bbox_mmxx_pxl - [bbox_mm, bbox_mm] + 1;
            
            input_is_cell_Q = iscell(tiles_value);
            if input_is_cell_Q
                im_class = class(tiles_value{1});
            else
                im_class = class(tiles_value);
            end
            % local_im = zeros(bbox_ll, im_class);
            local_im = nan(bbox_ll, 'single');
            for i_tiles = 1 : num_tiles
                tmp_bbox_mmxx = local_bbox_mmxx(i_tiles, :);
                
                if input_is_cell_Q
                    tmp_val = tiles_value{i_tiles};
                else
                    tmp_val = tiles_value(i_tiles);
                end
                if ~issingle(tmp_val)
                    tmp_val = single(tmp_val);
                end
                
                switch stitch_method
                    case 'overwrite'
                        local_im(tmp_bbox_mmxx(1) : tmp_bbox_mmxx(3), tmp_bbox_mmxx(2) : tmp_bbox_mmxx(4)) = ...
                            tmp_val;
                    case 'mip'
                        local_im(tmp_bbox_mmxx(1) : tmp_bbox_mmxx(3), tmp_bbox_mmxx(2) : tmp_bbox_mmxx(4)) = ...
                            max(local_im(tmp_bbox_mmxx(1) : tmp_bbox_mmxx(3), tmp_bbox_mmxx(2) : tmp_bbox_mmxx(4)), ...
                            tmp_val, 'omitnan');
                end
            end
            
            if nargout > 1
                varargout{1} = local_bbox_mmxx;
            end
        end
        
        function neighbor_add_coeff = get_neighbor_add_coeff(array_size, connectivity, remove_center_Q)
            if nargin < 3
                remove_center_Q = true;
            end
            ndim = numel(array_size);
            switch ndim
                case 2
                    switch connectivity
                        case 4
                            cube = zeros(3,3);
                            cube(:, 2) = 1;
                            cube(2, :) = 1;
                            [x,y] = ind2sub([3,3], find(cube(:)));
                        case 8
                            [x,y] = ndgrid(1:3, 1:3);
                        otherwise
                            error('Connectivity for 2D array should be 4 or 8');
                    end
                    neighbor_add_coeff = sub2ind(array_size, x, y);
                case 3
                    switch connectivity
                        case 26
                            [x,y,z] = ndgrid(1:3, 1:3, 1:3);
                        case 6
                            cube = zeros(3,3,3);
                            cube(2,2,:) = 1;
                            cube(:,2,2) = 1;
                            cube(2,:,2) = 1;
                            [x,y,z] = ind2sub([3,3,3], find(cube(:)));
                        case 18
                            cube = zeros(3,3,3);
                            cube(2,:,:) = 1;
                            cube(:,2,:) = 1;
                            cube(:,:,2) = 1;
                            [x,y,z] = ind2sub([3,3,3], find(cube(:)));
                        otherwise
                            error('Connectivity should be 6, 18 or 26');
                    end
                    neighbor_add_coeff = sub2ind(array_size, x, y, z);
            end
            middle_idx = connectivity/2 + 1;
            neighbor_add_coeff = neighbor_add_coeff - neighbor_add_coeff (middle_idx);
            if remove_center_Q
                neighbor_add_coeff(middle_idx) = [];
            end
            neighbor_add_coeff = neighbor_add_coeff(:);
        end
        
        function neighbor_ind = get_neighbor_grid_ind(array_size, connectivity,...
                ctr_ind)
            if isrow(ctr_ind)
                ctr_ind = ctr_ind.';
            end
            ctr_sub = Grid.ind2sub(array_size, ctr_ind);
            ctr_sub_p = ctr_sub + 1;
            add_coeff = Grid.get_neighbor_add_coeff(array_size, connectivity, true).';
            add_coeff_p = Grid.get_neighbor_add_coeff(array_size+2, connectivity, true).';
            ctr_ind_p = Grid.sub2ind(array_size+2, ctr_sub_p);
            neighbor_ind_p = bsxfun(@plus, ctr_ind_p, add_coeff_p);
            neighbor_ind = bsxfun(@plus, ctr_ind, add_coeff);
            neighbor_sub_p = Grid.ind2sub(array_size+2, neighbor_ind_p(:));
            is_valid_Q = all(neighbor_sub_p > 1, 2) & all(neighbor_sub_p < (array_size+2), 2);
            neighbor_ind(~is_valid_Q) = nan;
        end
        
        function overlapQ = check_bounding_box_overlap(bbox_mmxx_1, bboxes_mmxx_2)
            validateattributes(bbox_mmxx_1, {'numeric'}, {'vector'});
            num_dim = numel(bbox_mmxx_1) / 2;
            num_bbox = size(bboxes_mmxx_2, 1);
            assert(size(bboxes_mmxx_2, 2)/num_dim == 2);
            overlapQ = true(num_bbox, 1);
            for i = 1 : num_dim
                overlap_in_this_dim_Q = (bboxes_mmxx_2(:, i + num_dim) >= bbox_mmxx_1(i)) & ...
                    (bboxes_mmxx_2(:, i) <= bbox_mmxx_1(i + num_dim));
                overlapQ = overlapQ & overlap_in_this_dim_Q;
            end            
        end
    end
end