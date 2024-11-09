classdef WBIMImagingGrid2D < Grid
    
   properties
       name char
       grid_axis_order (1,2) double
       acq_grid_axis_order (1,2) double
       stack_size (1,:) double
       stack_size_um (1,:) double
       pixel_size_um (1, :) double
       stack_overlap_size (1, :) double
       stack_overlap_um (1, :) double
       z_list_um (:, 1) double
       num_z_section (1,1) double
       
       connectivity (1,1) double = 8
       
       % The following properties are initialized in WBIMImaging
       piezo_z_list_um (:, 1) double 
       scan_angle_xy (1, 2) double
   end
   
   properties(Hidden)
       priority_score_map double
   end
   
   methods
       function obj = WBIMImagingGrid2D(grid_name, space_xyz_size, tile_xyz_size,...
                pxl_xyz_size_um, overlap_xyz_size, grid_axis_order, acq_grid_axis_order)           
            % If grid_axis_order = [2,1], the first axis of the grid is the
            % y axis and the second axis of the grid is the x axis - to
            % facilitate image stitching 
            % This grid is in the sample coordinate. 
            obj @ Grid(space_xyz_size(grid_axis_order),...
                tile_xyz_size(grid_axis_order), overlap_xyz_size(grid_axis_order));
            
            obj.name = grid_name;
            obj.grid_axis_order = grid_axis_order;
            obj.acq_grid_axis_order = acq_grid_axis_order;
            
            obj.stack_size = tile_xyz_size([grid_axis_order, 3]);
            obj.stack_overlap_size = overlap_xyz_size([grid_axis_order, 3]);
            obj.pixel_size_um = pxl_xyz_size_um([grid_axis_order, 3]);
            obj.init_parameters();
       end
       
       function init_parameters(obj)
           obj.stack_size_um = obj.stack_size .* obj.pixel_size_um;
           obj.z_list_um = (0 : (obj.stack_size(3) - 1)) * obj.pixel_size_um(3);
           obj.num_z_section = numel(obj.z_list_um);
           obj.stack_overlap_um = obj.stack_overlap_size .* obj.pixel_size_um; 
           
           obj.priority_score_map = obj.get_priority_score_map(obj.grid_size, ...
               obj.acq_grid_axis_order);
       end
       
       function info = get_stack_info(obj, grid_ind)
           info = struct;
           info.grid_size = obj.grid_size;
           info.grid_ind = grid_ind;
           info.grid_sub = obj.grid_ind2sub(grid_ind);
           info.stack_size = obj.stack_size;
           info.stack_size_um = obj.stack_size_um;
           info.pixel_size_um = obj.pixel_size_um;
           info.overlap_size = obj.stack_overlap_size;
           info.overlap_size_um = obj.stack_overlap_um;
           info.center_sub = obj.center_sub(grid_ind, :);
           
           info.center_sub_um = info.center_sub .* obj.pixel_size_um(1:2);
           info.center_xy_um = info.center_sub_um(:, obj.grid_axis_order);

           info.tile_mmll_pxl = obj.mmll(grid_ind, :);
           info.tile_mmxx_pxl = obj.mmxx(grid_ind, :);  
           
           info.tile_mmll_um = info.tile_mmll_pxl .* ...
               [obj.pixel_size_um(1:2), obj.pixel_size_um(1:2)];
           info.tile_mmxx_um = info.tile_mmxx_pxl .* ...
               [obj.pixel_size_um(1:2), obj.pixel_size_um(1:2)];
           
           info.neighbor_grid_ind = obj.get_neighbor_grid_ind(grid_ind);
       end
       
       function grid_ind = sort_grid_ind(obj, grid_ind)
           if numel(grid_ind) < 3
               return;
           end
           if obj.name == WBIMMicroscopeMode.Scan
               grid_ind = obj.sort_grid_ind_by_priority_map(grid_ind);
               grid_ind = obj.sort_grid_ind_by_shortest_path(grid_ind);
           elseif obj.name == WBIMMicroscopeMode.Explore
%                candidate_map = false(obj.grid_size);
%                candidate_map(grid_ind) = true;
%                candidate_cc = bwconncomp(candidate_map);
%                cc_stat = regionprops(candidate_cc, '');
%                num_cc = candidate_cc.NumObjects;
                grid_ind = obj.sort_grid_ind_by_shortest_path(grid_ind);               
           end
       end
       
       function grid_ind = sort_grid_ind_by_priority_map(obj, grid_ind)
           [~, sorted_idx] = sort(obj.priority_score_map(grid_ind), 'ascend');
           grid_ind = grid_ind(sorted_idx);
       end
       
       function grid_ind = sort_grid_ind_by_shortest_path(obj, grid_ind, current_grid_ind)
           if nargin < 3
               current_grid_ind = [];
           else
               validateattributes(current_grid_ind, 'numeric', {'positive', 'integer'});
           end
           visQ = false;
           % Use 2opt algorithm to solve the traveling salesman problem
           % without returning to the initial position. 
           % If current_grid_ind is given, force the sequence to start from
           % the current_grid_ind, but return the grid_ind without
           % current_grid_ind
           % Sort first according to the priority map
%            grid_ind = obj.sort_grid_ind_by_priority_map(grid_ind);      
%            grid_ind = grid_ind(randperm(numel(grid_ind)));
           grid_ind = cat(1, current_grid_ind, grid_ind(:));
           grid_sub = obj.grid_ind2sub(grid_ind);
           dist_mat = squareform(pdist(grid_sub, 'seuclidean', [0.75, 1]));
           num_nodes = size(dist_mat, 1);
           % Add one dummy node to decouple the first node and the last
           % node in the chain
           dist_mat = cat(1, cat(2, dist_mat, zeros(num_nodes, 1)), zeros(1, num_nodes+1));
           if ~isempty(current_grid_ind)
               dist_mat = cat(1, cat(2, dist_mat, inf(num_nodes+1, 1)), inf(1, num_nodes+2));
               dist_mat(num_nodes+2, 1) = 0;
               dist_mat(num_nodes+1, num_nodes+2) = 0;
           end
           sorted_grid_idx = tsp_shortest_path_A_2opt(dist_mat, visQ);
           sorted_grid_idx = sorted_grid_idx(sorted_grid_idx <= num_nodes);
           if visQ
               fig_hdl = figure;
               ax_hdl = axes(fig_hdl);
               scatter(ax_hdl, grid_sub(sorted_grid_idx, 2), grid_sub(sorted_grid_idx, 1),...
                   100, 1:num_nodes, 'fill');
               hold(ax_hdl, 'on');
               plot(ax_hdl, grid_sub(sorted_grid_idx, 2), grid_sub(sorted_grid_idx, 1));
               ax_hdl.XLabel.String = 'X';
               ax_hdl.YLabel.String = 'Y';
               cbar_hdl = colorbar(ax_hdl);
               cbar_hdl.Label.String = 'Visiting order';
               grid(ax_hdl, 'on');
           end
           if ~isempty(current_grid_ind)
               assert(sorted_grid_idx(1) == 1);
               sorted_grid_idx = sorted_grid_idx(2:end) - 1;
           end
           grid_ind = grid_ind(sorted_grid_idx);
       end
   end
   %% Coordinate transform
   methods
       function sub = xy_um_to_sub(obj, xy_vstack)
           sub = ceil(xy_vstack(:, obj.grid_axis_order) ./ ...
               obj.pixel_size_um(1:2));
       end
       
       function grid_sub = xy_um_to_grid_sub(obj, xy_vstack)
           pos_pxl = obj.xy_um_to_sub(xy_vstack);
           grid_sub = obj.pos_to_grid_sub(pos_pxl);
       end
       
       function grid_ind = xy_um_to_grid_ind(obj, xy_vstack)
           pos_pxl = obj.xy_um_to_sub(xy_vstack);
           grid_ind = obj.pos_to_grid_ind(pos_pxl);
       end
   end
   %% Utilities    
   methods
       % Stitch multiple tiles
       function tiles_info = get_tiles_local_bbox(obj, grid_ind)
           bbox_mmxx_list = obj.mmxx(grid_ind, :);
           
           tiles_info = struct;           
           tiles_info.pixel_size_um = obj.pixel_size_um(1:2);
           tiles_info.num_tiles = numel(grid_ind);
           tiles_info.g_grid_ind = grid_ind;
           tiles_info.g_grid_sub = obj.grid_ind2sub(grid_ind);
           
           tiles_info.g_bbox_mm = min(bbox_mmxx_list(:, 1:2), [], 1);
           tiles_info.g_bbox_xx = max(bbox_mmxx_list(:, 3:4), [], 1);
           tiles_info.g_bbox_mm_um = tiles_info.g_bbox_mm .* tiles_info.pixel_size_um;
           tiles_info.g_bbox_xx_um = tiles_info.g_bbox_xx .* tiles_info.pixel_size_um;
           
           tiles_info.g_bbox_mmxx = cat(2, tiles_info.g_bbox_mm, ...
               tiles_info.g_bbox_xx);
           tiles_info.g_bbox_mmll = [tiles_info.g_bbox_mmxx(1:2), ...
               tiles_info.g_bbox_mmxx(3:4) - tiles_info.g_bbox_mmxx(1:2) + 1];
           tiles_info.g_grid_bbox_mmxx = [min(tiles_info.g_grid_sub, [], 1), ...
               max(tiles_info.g_grid_sub, [], 1)];
           tiles_info.g_grid_bbox_mmll = [tiles_info.g_grid_bbox_mmxx(1:2), ...
               tiles_info.g_grid_bbox_mmxx(3:4) - tiles_info.g_grid_bbox_mmxx(1:2) + 1];
           
           % Local bbox
           tiles_info.bbox_mmxx_list = bsxfun(@minus, bbox_mmxx_list, ...
               [tiles_info.g_bbox_mm, tiles_info.g_bbox_mm] - 1);
       end
       
       function [local_im, tiles_info] = get_local_tile_image(obj, grid_ind, tiles_value)
           num_tiles = numel(grid_ind);
           assert(num_tiles == numel(tiles_value), 'One tile index for one tile value');
           if num_tiles == 0
               local_im = [];
               return;
           end
           tiles_info = obj.get_tiles_local_bbox(grid_ind);
           [local_im, ~] = obj.stitch_tiles(tiles_info.bbox_mmxx_list, ...
               tiles_value);
       end
       
       function neighbor_grid_ind = get_neighbor_grid_ind(obj, ctr_ind, conn)
           if nargin < 3
               conn = obj.connectivity; 
           end
           neighbor_grid_ind = get_neighbor_grid_ind@Grid(obj.grid_size, ...
               conn, ctr_ind);
       end
       
       function grid_ind = get_tiles_in_polygon_um(obj, vertex_yx_um)
           grid_ctr_yx_um = obj.center_sub .* obj.pixel_size_um(1:2);
           [inQ, onQ] = inpolygon(vertex_yx_um(:, 2), vertex_yx_um(:, 1), ...
               grid_ctr_yx_um(:, 2), grid_ctr_yx_um(:, 1));
           grid_ind = find(inQ | onQ);
       end
       
       function grid_ind = get_tiles_with_center_in_bbox_um(obj, bbox_yx_mmxx_um)
           assert(numel(bbox_yx_mmxx_um) == 4 &&...
               all(bbox_yx_mmxx_um(1:2) <= bbox_yx_mmxx_um(3:4)));
           if iscolumn(bbox_yx_mmxx_um)
               bbox_yx_mmxx_um = bbox_yx_mmxx_um.';
           end
           bbox_yx_mmxx = bbox_yx_mmxx_um ./ [obj.pixel_size_um(1:2), obj.pixel_size_um(1:2)];
           grid_ind = obj.get_tiles_with_center_in_bbox(bbox_yx_mmxx);           
       end
       
       function grid_ind = get_tiles_completely_in_bbox_um(obj, bbox_yx_mmxx_um)
           arguments
               obj (1,1) WBIMImagingGrid2D
               bbox_yx_mmxx_um (1, 4) {mustBeNumeric, mustBePositive}
           end
           bbox_yx_mmxx_pxl = bbox_yx_mmxx_um ./ [obj.pixel_size_um(1:2), obj.pixel_size_um(1:2)];
           grid_ind = obj.get_tiles_completely_in_bbox(bbox_yx_mmxx_pxl);           
       end
       
       function grid_ind = get_tiles_overlap_with_bbox_um(obj, bbox_yx_mmxx_um)
           arguments
               obj (1,1) WBIMImagingGrid2D
               bbox_yx_mmxx_um (1, 4) {mustBeNumeric, mustBePositive}
           end
           bbox_yx_mmxx_pxl = bbox_yx_mmxx_um ./ [obj.pixel_size_um(1:2), obj.pixel_size_um(1:2)];
           grid_ind = obj.get_tiles_overlap_with_bbox(bbox_yx_mmxx_pxl);                      
       end
       
       function tile_info = get_tile_info_at_xy_um(obj, xy_um)
           arguments
               obj (1,1) WBIMImagingGrid2D
               xy_um (:, 2) 
           end
           grid_ind = obj.xy_um_to_grid_ind(xy_um);
           tile_info = obj.get_stack_info(grid_ind);           
       end
   end
   %% 
   methods(Static)
       function psm = get_priority_score_map(grid_size, axis_order)
           if nargin < 2
               axis_order = [1,2];
           end
           num_grid_cell = prod(grid_size);
           psm = reshape(1 : num_grid_cell, grid_size(axis_order));
           for i_col = 1 : grid_size(axis_order(2))
               if mod(i_col, 2) == 0
                   psm(:, i_col) = flip(psm(:, i_col));
               end
           end
           psm = permute(psm, axis_order);
       end
       
       function c_info = add_combined_region_info(c_info)
           tile_dim = size(c_info.tile_mmxx_pxl, 2)/2;
           % Local grid coordinate
           c_info.region_grid_bbox_mm = min(c_info.grid_sub, [], 1);
           c_info.region_grid_bbox_xx = max(c_info.grid_sub, [], 1);
           c_info.region_grid_size = c_info.region_grid_bbox_xx - ...
               c_info.region_grid_bbox_mm + 1;
           c_info.region_grid_sub = c_info.grid_sub - c_info.region_grid_bbox_mm + 1;
           
           c_info.region_grid_ind = Grid.sub2ind(c_info.region_grid_size, ...
               c_info.region_grid_sub);
           
           % Image space coordiante in pixel
           c_info.region_bbox_mm_pxl = min(c_info.tile_mmxx_pxl(:, 1:tile_dim), ...
               [], 1);
           c_info.region_bbox_xx_pxl = max(c_info.tile_mmxx_pxl(:, (tile_dim+1):(2 * tile_dim)), ...
               [], 1);
           c_info.local_tile_bbox_mmxx_pxl = c_info.tile_mmxx_pxl - ...
               [c_info.region_bbox_mm_pxl, c_info.region_bbox_mm_pxl] + 1;
           assert(all(c_info.local_tile_bbox_mmxx_pxl > 0, 'all'));
           c_info.region_bbox_ll_pxl = c_info.region_bbox_xx_pxl - ...
               c_info.region_bbox_mm_pxl + 1;
           
           % Image space coordiante in micron
           if isrowvector(c_info.pixel_size_um)
               c_info.region_bbox_mm_um = c_info.region_bbox_mm_pxl .* ...
                   c_info.pixel_size_um(1:tile_dim);
               c_info.region_bbox_xx_um = c_info.region_bbox_xx_pxl .* ...
                   c_info.pixel_size_um(1:tile_dim);
               c_info.region_bbox_ctr_yx_um = (c_info.region_bbox_xx_um + ...
                   c_info.region_bbox_mm_um) ./2;
               c_info.region_bbox_ll_um = c_info.region_bbox_ll_pxl .* ...
                   c_info.pixel_size_um(1:tile_dim);
               c_info.local_tile_bbox_mmxx_um = c_info.local_tile_bbox_mmxx_pxl .* ...
                   [c_info.pixel_size_um(1:tile_dim), c_info.pixel_size_um(1:tile_dim)];
           else
               warning('Tiles have different pixel size');
           end
       end
       
       function c_info = add_combined_region_info_3D(c_info)
           % Overwrite grid 
           c_info.num_tiles = size(c_info.grid_sub, 1);
           c_info.num_layers = numel(unique(c_info.layer));
           if isscalar(c_info.layer)
               c_info.layer = repelem(c_info.layer, c_info.num_tiles, 1);
           end
           c_info.grid_sub = cat(2, c_info.grid_sub, c_info.layer);
           c_info.grid_size = [c_info.grid_size, c_info.num_layers];
           c_info = rmfield(c_info, {'grid_ind', 'layer', 'num_layers'});
           
           tile_z_mm_um = c_info.sample_xyz_um_done(:, 3);
           tile_z_mm_pxl = round(tile_z_mm_um / c_info.pixel_size_um(3));
           tile_mmll_pxl = cat(2, c_info.tile_mmll_pxl(:, 1:2), tile_z_mm_pxl, ...
               repmat(c_info.stack_size, c_info.num_tiles, 1));
           c_info.tile_mmll_pxl = tile_mmll_pxl;
           c_info.tile_mmxx_pxl = cat(2, c_info.tile_mmll_pxl(:, 1:3), ...
               c_info.tile_mmll_pxl(:, 1:3) + c_info.tile_mmll_pxl(:, 4:6) - 1);           
           c_info.tile_mmll_um = tile_mmll_pxl .* [c_info.pixel_size_um, c_info.pixel_size_um];
           c_info.tile_mmxx_um = c_info.tile_mmxx_pxl .* [c_info.pixel_size_um, c_info.pixel_size_um];
           
           c_info.region_grid_bbox_mm = min(c_info.grid_sub, [], 1);
           region_grid_bbox_xx = max(c_info.grid_sub, [], 1);
           c_info.region_grid_size = region_grid_bbox_xx - ...
               c_info.region_grid_bbox_mm + 1;
           c_info.region_grid_sub = c_info.grid_sub - c_info.region_grid_bbox_mm + 1;
           
           c_info.region_grid_ind = Grid.sub2ind(c_info.region_grid_size, ...
               c_info.region_grid_sub);
           c_info.region_label_array = zeros(c_info.region_grid_size);
           c_info.region_label_array(c_info.region_grid_ind) = 1 : ...
               c_info.num_tiles;
           
           region_bbox_mm_pxl = min(tile_mmll_pxl(:, 1:3), [], 1);
           c_info.region_bbox_mmxx_pxl = c_info.tile_mmxx_pxl - ...
               [region_bbox_mm_pxl, region_bbox_mm_pxl] + 1;           
       end       
   end
   
end