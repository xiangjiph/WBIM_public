classdef WBIMTileManager < matlab.mixin.Copyable
    % To be implemented
    % 1. Display added tiles, imaged tiles in grid / pixel resolution
    % 2. GUI
    %    a. Add tiles by clicking (multi-selections)
    %    b. Deleted un-imaged, added tiles
    % 3. Sort tiles - priority queue implementation
    %
    %% Porperties
    properties
        experiment_group string
        experiment string 
    end
    properties(SetAccess=protected)
       name char
    end
    
    properties        
        layer (1,1) uint16 = 1
        tile_info (:, 1) WBIMTileMetadata = WBIMTileMetadata.empty();
%         h_viewer (1,1) WBIMTileViewer
        num_processed_tile (1,1) uint16 = 0
        tile_properties struct        
        
        grid WBIMImagingGrid2D % Tiles are in sample coordiante
    end
    
    properties(SetObservable)
        % arraies to keep track of the acquisition progress
        candidate_map logical
        % each element is the number of times the tile has been acquired. -
        % unlikely to image the same tile for > 255 times
        imaged_map uint8
    end
    
    properties(Transient, Hidden)
        candidate_map_im logical
        imaged_map_im uint8
        local_mip_yx_um cell
        h_data_manager (1,1) WBIMFileManager = WBIMFileManager
    end
    
    properties(Hidden)                
        % Recording
        lc_tile (1, :) cell
        lc_tile_properties (1, :) cell
        lc_num_tile (1, :) uint16
        lc_num_tile_processed (1, :) uint16
        lc_imaged_map (1, :) cell
        lc_candidate_map (1, :) cell
    end
    
    properties(Dependent)
       num_layer (1,1) double  
    end
    %% Visualization
    properties
%         fig_hdl
%         ax_candidate
%         ax_imaged        
    end
    %% Basic class function
    methods
        function obj = WBIMTileManager(exp_group, exp_name, acq_grid)
            arguments
                exp_group
                exp_name
                acq_grid
            end
            % Do not copy. Pass by reference so that the tile manager grid
            % can get updated as the grid in WBIMImaging
            % But this sound dangerous... 
            obj.experiment_group = exp_group;
            obj.experiment = exp_name;
            obj.grid = acq_grid;
            obj.name = acq_grid.name;
            obj.init_layer();
        end
        
        function delete(obj)
            % To be implemented
        end
        
        function save(obj)
            fp = obj.h_data_manager.fp_tile_manager(obj.experiment_group, ...
                obj.experiment, obj.name);
            [folder, ~] = fileparts(fp);
            if ~isfolder(folder)
                mkdir(folder);
            end
            save(fp, 'obj');
        end
        
        function save_str = saveobj(obj)
           save_str = struct;
           save_str.experiment_group = obj.experiment_group;
           save_str.experiment = obj.experiment;
           save_str.grid = obj.grid;
        end
        
        function compress(obj)
            % Delete buffered data for temporarily storage in the RAM
            obj.local_mip_um = [];
            obj.local_mip_yx_um = [];
            for i = 1 : numel(obj.tile_info)
                obj.tile_info(i).clear_buffer();
            end
        end
    end
    
    methods(Static)
        function obj = loadobj(save_str)
           obj = WBIMTileManager(save_str.experiment_group, ...
               save_str.experiment, save_str.grid);
        end
        
        function obj = load(fp)
            obj = load(fp);
            obj = obj.obj;
        end
    end
    %% Gets and sets
    methods
        function val = get.num_layer(obj)
            val = numel(obj.lc_tile);
        end
    end
    %% Basic operations
    methods
        %% Candidate map
        function add_tile_at_xy(obj, xy_um_vstack, allow_multiple_acqQ)
%             xy_um_vstack is in sample coordinate
            if nargin < 3
                allow_multiple_acqQ = true;
            end
            assert(size(xy_um_vstack, 2) == 2);
            grid_ind = obj.grid.xy_um_to_grid_ind(xy_um_vstack);
            obj.add_tile_by_grid_ind(grid_ind, allow_multiple_acqQ);
        end
        
        function add_tile_by_grid_ind(obj, grid_ind, allow_multiple_acqQ)
            validateattributes(grid_ind, {'numeric'}, {'vector', 'nonnegative'});
            if nargin < 3
                allow_multiple_acqQ = true;
            end
            if ~allow_multiple_acqQ
                grid_ind = grid_ind(~obj.imaged_map(grid_ind));
            end
            % Effectively add a list of unique grid indices, i.e.,
            % duplicated indices are removed
            obj.candidate_map(grid_ind) = true;            
        end
        
        function add_tile_by_grid_sub(obj, grid_sub, allow_multiple_acqQ)
            if nargin < 3
                allow_multiple_acqQ = true;
            end
            grid_ind = obj.grid.grid_sub2ind(grid_sub);
            obj.add_tile_by_grid_ind(grid_ind, allow_multiple_acqQ);            
        end

        function add_tile_in_bbox_um(obj, bbox_mmxx_um)
            % bbox_mmxx_um: [y_min, x_min, y_max, x_max]
            grid_ind = obj.grid.get_tiles_with_center_in_bbox_um(bbox_mmxx_um);
            obj.add_tile_by_grid_ind(grid_ind, false);            
        end
        
        function remove_candidate_tile_at_xy(obj, xy_um_vstack)
            assert(size(xy_um_vstack, 2) == 2);
            grid_ind = obj.grid.xy_um_to_grid_ind(xy_um_vstack);
            obj.candidate_map(grid_ind) = false;
        end
        
        function grid_ind = candidate_map_to_grid_ind(obj, allow_multipe_acqQ, clear_map_Q)
            arguments
                obj (1,1) WBIMTileManager
                allow_multipe_acqQ (1,1) logical = true;
                clear_map_Q (1,1) logical = false;
            end
            grid_ind = find(obj.candidate_map);
            if ~allow_multipe_acqQ
               grid_ind = grid_ind(obj.imaged_map(grid_ind) == 0); 
            end            
            if clear_map_Q
                obj.clear_candidate_map();
            end
        end
        
        function clear_candidate_map(obj)
            obj.candidate_map = false(obj.grid.grid_size);            
        end
        
        function expand_candidate_map(obj, expand_r, expand_dir)
            arguments
                obj WBIMTileManager
                expand_r (1,1) double
                expand_dir = 'all';
            end
            obj.candidate_map = obj.expand_map(obj.candidate_map, ...
                expand_r, expand_dir);
        end
        
        function add_imaged_tile_neighbor_to_candidate_map(obj, expand_r, expand_dir)
           if nargin < 3
               expand_dir = 'all';
           end
           is_imaged_Q = obj.imaged_map > 0;
           obj.candidate_map = ~is_imaged_Q & ...
               obj.expand_map(is_imaged_Q, expand_r, expand_dir);           
        end        
        %% 
        function add_imaged_tile(obj, tile_info)
            if ~isempty(tile_info)
                imaged_grid_ind = [tile_info.grid_ind];
                tile_layer = [tile_info.layer];
                assert(all(tile_layer == obj.layer));
                has_imaged_Q = obj.imaged_map(imaged_grid_ind) > 0;
                obj.tile_info = cat(1, obj.tile_info, tile_info);
                obj.imaged_map(imaged_grid_ind) = obj.imaged_map(imaged_grid_ind) + 1;
            end
            % TODO: Check duplicated tiles by time - only useful when
            % loading tile metadata from disk 
        end
        
        function tile = get_latest_added_tile(obj)
            tile = obj.tile_info(end).copy();
        end
        
        function tile_info = get_unprocessed_tile_info(obj)
            tile_info = obj.tile_info((obj.num_processed_tile + 1) : ...
                numel(obj.tile_info));
        end
        
        function broken_tiles = remove_unprocessed_tile_info(obj, add_to_candidate_Q)
            if nargin < 2
                add_to_candidate_Q = false;
            end
            broken_tiles_Q = ~[obj.tile_info.file_processed_Q] & ...
                ~[obj.tile_info.init_SI_file_exist_Q];
            broken_tiles = obj.tile_info(broken_tiles_Q);
            broken_ind = [broken_tiles.grid_ind];
            obj.imaged_map(broken_ind) = 0;
            obj.tile_info = obj.tile_info(~broken_tiles_Q);            
            
            if add_to_candidate_Q
                obj.add_tile_by_grid_ind(broken_ind);
            end
        end
        %% Get tiles in layer without changing layer setting
        function tile_info = get_latest_tiles_in_layer(obj, layer_idx)
            if layer_idx == obj.layer
                tile_info = obj.select_duplicated_tiles(obj.tile_info);
            else
                if ~isempty(obj.lc_tile{layer_idx})
                    tile_info = obj.select_duplicated_tiles(obj.lc_tile{layer_idx});
                else
                    tile_info = [];
                end
            end            
        end        
        %% Change layers
        function init_layer(obj)
            obj.clear_candidate_map();
            obj.imaged_map = zeros(obj.grid.grid_size, 'uint8');
            obj.tile_info = WBIMTileMetadata.empty();
            obj.num_processed_tile = 0;
            obj.tile_properties = [];
        end
        
        function set.layer(obj, val)
            obj.switch_to_layer(val);            
            obj.layer = val;
        end
        
        %% Tile properties
        function reset_tile_properties(obj)
           obj.num_processed_tile = 0;
           obj.tile_properties = [];
        end
        
        %% Delete tiles
        function clear_tile_in_layer(obj, layer_idx)
            num_del_layer = numel(layer_idx);
            for i = 1 : num_del_layer
                obj.lc_tile{layer_idx(i)} = {};
            end            
        end
        
        function folder_pair = get_tile_sync_folder_pairs(obj, layer_idx)
            if nargin < 2
                layer_idx = obj.layer;
            end
            if layer_idx == obj.layer
                tiles = obj.tile_info;
            else
                tiles = obj.lc_tile{layer_idx};
            end
            folder_pair = arrayfun(@(x) {sprintf('%s\\', x.get_tile_folder); ...
                x.fp_birdstore_folder}, tiles, 'UniformOutput', false);
            folder_pair = cat(2, folder_pair{:});
        end
    end
    
    
    methods(Access=private)
        function switch_to_layer(obj, layer_idx)
            % Called by the set method of layer. Do not call this function
            validateattributes(layer_idx, 'numeric', {'scalar', 'positive', 'integer'});
            if layer_idx ~= obj.layer
                obj.archieve_layer_info(obj.layer);
                if numel(obj.lc_tile) >= layer_idx && ~isempty(obj.lc_tile{layer_idx})
                    obj.load_archieved_layer_info(layer_idx);
                else
                    obj.init_layer();
                end
            end
        end
        
        function archieve_layer_info(obj, layer_idx)
           if nargin < 2
               layer_idx = obj.layer;
           end
           % Copy to clear buffered data
           obj.lc_tile{layer_idx} = obj.tile_info.copy(); 
           obj.lc_imaged_map{layer_idx} = obj.imaged_map;
           if any(obj.candidate_map)
               obj.lc_candidate_map{layer_idx} = obj.candidate_map;
           else
               obj.lc_candidate_map{layer_idx} = [];
           end
           obj.lc_num_tile(layer_idx) = numel(obj.tile_info);
           obj.lc_num_tile_processed(layer_idx) = obj.num_processed_tile;   
           obj.lc_tile_properties{layer_idx} = obj.tile_properties;
        end
        
        function load_archieved_layer_info(obj, layer_idx)
            assert(layer_idx > 0 && layer_idx <= numel(obj.lc_tile), 'Layer index is out of range');
            if isempty(obj.lc_imaged_map{layer_idx})
                assert(isempty(obj.lc_tile{layer_idx}) && ...
                    obj.lc_num_tile_processed(layer_idx) == 0 && ...
                    isempty(obj.lc_tile_properties{layer_idx}));
                fprintf('Layer %d has not been initialized yet. Initializing layer...\n', ...
                    layer_idx);
                obj.init_layer();
            else                
                obj.tile_info = obj.lc_tile{layer_idx};
                obj.num_processed_tile = obj.lc_num_tile_processed(layer_idx);
                obj.tile_properties = obj.lc_tile_properties{layer_idx};
                obj.imaged_map = obj.lc_imaged_map{layer_idx};
                obj.candidate_map = obj.lc_candidate_map{layer_idx};
                if isempty(obj.candidate_map)
                   obj.candidate_map = false(obj.grid.grid_size); 
                end
            end
        end
    end
    
    methods(Static, Hidden)
        function new_map = expand_map(old_map, expand_r, expand_dir)
            if nargin < 3
                expand_dir = 'all';
            end
            switch expand_dir
                case 'all'
                    expand_filter = strel('square', expand_r * 2 + 1).Neighborhood;
                case 'left'
                    expand_filter = cat(2, true(1, expand_r + 1), false(1, expand_r));
                case 'right'
                    expand_filter = cat(2, false(1, expand_r), true(1, expand_r + 1));
                case {'up', 'top', 'above'}
                    expand_filter = cat(1, true(expand_r + 1, 1), false(expand_r, 1));
                case {'down', 'button', 'below'}
                    expand_filter = cat(1, false(expand_r, 1), true(expand_r + 1, 1));
                otherwise
                    error('Unrecognized expand direction')
            end
            % Flip for convolution
            boundary_map = convn(old_map, expand_filter, 'same');
            boundary_map = (boundary_map > 0) & (boundary_map <= expand_r) & ...
                ~old_map;
            new_map = boundary_map | old_map;
        end
    end
    %% Analyze imaged tiles
    methods                    
        function [grid_ind, varargout] = get_tile_grid_ind_in_imaged_mask(obj,...
                local_mask_yx_um, local_tile_info)
            local_mmxx_um = round(local_tile_info.region_bbox_mm_um);
            mask_size = size(local_mask_yx_um);
            assert(all(abs(mask_size - local_tile_info.region_bbox_ll_um) < 1), 'Inconsistent mask size');
            new_bbox_mmxx_um = round(obj.grid.mmxx .* [obj.grid.pixel_size_um(1:2), ...
                obj.grid.pixel_size_um(1:2)]);
            new_bbox_local_mmxx_um = new_bbox_mmxx_um - [local_mmxx_um, local_mmxx_um];
            
            is_in_mask_Q = fun_check_bounding_boxes_overlap([1,1,mask_size], ...
                new_bbox_local_mmxx_um);            
            
            grid_ind = find(is_in_mask_Q);
            new_bbox_local_mmxx_um = new_bbox_local_mmxx_um(grid_ind, :).';
            new_bbox_local_mmxx_um = min(max(new_bbox_local_mmxx_um, 1), [mask_size.'; mask_size.']);
            
            num_tiles = numel(grid_ind);
            in_mask_fraction = nan(num_tiles, 1);
            for i = 1 : num_tiles
                tmp_bbox = new_bbox_local_mmxx_um(:, i);
                tmp_local_mask = local_mask_yx_um(tmp_bbox(1) : tmp_bbox(3), ...
                    tmp_bbox(2):tmp_bbox(4));
                in_mask_fraction(i) = nnz(tmp_local_mask) / numel(tmp_local_mask);
            end
            is_in_mask_Q = in_mask_fraction > 0;
            
            grid_ind = grid_ind(is_in_mask_Q);
            if nargout == 2
                varargout{1} = in_mask_fraction(is_in_mask_Q);
            end
        end
        
        function [local_mask_yx_um, local_tile_info] = get_roi_mask_from_tiles(...
                obj, tile_info, channel_id, visQ)
            % TODO: this function should be moved to somewhere else
            if nargin < 4
                visQ = false;
            end
            if isempty(tile_info)
                local_mask_yx_um = [];
                local_tile_info = [];
            else
                [local_im, local_tile_info] = obj.get_stitched_tile_mip_yx_um(...
                    tile_info, channel_id);
                % To be improve...
                local_im = medfilt2(local_im, [5,5]);
                local_mask_yx_um = local_im > 3500; % Change from 2500. For 30% excitation. 05/02/2024 
                min_cc_fraction = 0.001;
                min_cc_size = min(900, round(numel(local_mask_yx_um) * min_cc_fraction));
                local_mask_yx_um = bwareaopen(local_mask_yx_um, min_cc_size);
                local_mask_yx_um = imdilate(local_mask_yx_um, strel('disk', 50));
%                 local_mask_yx_um = imfill(local_mask_yx_um, 'holes');
                
                if visQ
                    fig_hdl = figure;
                    ax_hdl = axes(fig_hdl);
                    imshowpair(local_im, local_mask_yx_um);
                    ax_hdl.DataAspectRatio = [1,1,1];
                    ax_hdl.XLabel.String = 'Sample X (\mum)';
                    ax_hdl.YLabel.String = 'Sample Y (\mum)';
                    ax_hdl.XAxis.Visible = 'on';
                    ax_hdl.YAxis.Visible = 'on';
                end
            end
        end
    end
    
    methods(Static)
        
    end
    %% Utilities
    methods
        function doneQ = check_tile_processing_state(obj)
            max_waiting_time_s = 300;
            t0 = tic;
            doneQ = false;
            if ~isempty(obj.tile_info)
                while ~doneQ
                    doneQ = all([obj.tile_info.file_processed_Q]);
                    if ~doneQ
                        t_diff_s = toc(t0);
                        if t_diff_s < max_waiting_time_s
                            fprintf('\tWBIMTileManager: Exist unprocessed tiles. Wait for 3 seconds...\n');
                            pause(3.0);
                        else
                            break;
                        end
                    else
                        % Wait for extra 5 seconds
                        pause(5.0);
                    end
                end
            end
        end
        
        function resubmit_processing_task(obj, tiles, options)
            arguments
               obj
               tiles = obj.tile_info;
               options.overwriteQ (1,1) logical = false;
            end
            
            for i = 1 : numel(tiles)
                tmp_fp = fullfile(tiles(i).save_root_folder, ...
                    tiles(i).fprr_tile);
                if ~isfile(tmp_fp) || options.overwriteQ
                    switch tiles(i).acq_mode
                        case WBIMMicroscopeMode.Scan
                            WBIMProcessingPython.processing_scan_data(tiles(i));
                            pause(8);
                        case WBIMMicroscopeMode.Explore
                            WBIMProcessExploreData(tiles(i));
                            pause(0.5);
                    end
                    fprintf('Reprocessing tile %d\n', i);
                end
            end
        end
        
        function [local_im, varargout] = map_to_image(obj, map_array)
            % Expand the grid array. One grid cell is expanded into a tile
            % of size specificed by the bounding box
            grid_ind = find(map_array);
            [local_im, tiles_info] = obj.grid.get_local_tile_image(grid_ind, ...
                map_array(grid_ind));
            if nargout > 1
                varargout{1} = tiles_info;
            end
        end
        
        function [local_im, varargout] = get_latest_imaged_mip_yx_um_in_layer(obj, ...
                channel_id, layer_idx)
            arguments
                obj
                channel_id (1, :) double = [1,2];
                layer_idx (1,1) double = obj.layer
            end
            latest_tile_info = obj.get_latest_tiles_in_layer(layer_idx);
            [local_im, local_tile_info] = obj.get_stitched_tile_mip_yx_um(latest_tile_info, channel_id);
            if nargout > 1
                varargout{1} = local_tile_info;
            end            
        end
        
        function latest_tile_info = get_latest_tiles_info(obj)
            latest_tile_info = obj.select_duplicated_tiles(obj.tile_info, 'latest');
        end
        
        function check_unprocessed_tile(obj)
            for i = 1 : numel(obj.tile_info)
                tmp_info = obj.tile_info(i);
                if ~isfile(fullfile(tmp_info.save_root_folder, tmp_info.fprr_tile))
                    WBIMProcessingPython.processing_scan_data(tmp_info);
                    fprintf('Starting processing tile %d\n', i);
                end
            end                        
        end
        
        function load_tiles_in_layer(obj, layer_idx)
%             if layer_idx ~= obj.layer
%                 obj.switch_to_layer(layer_idx);
%             end
            current_layer = obj.layer;
            obj.layer = layer_idx;
            tile_list = obj.h_data_manager.load_tile_in_layer(obj.experiment_group, ...
                obj.experiment, layer_idx, obj.name);
            obj.add_imaged_tile(tile_list);
            obj.layer = current_layer;
        end
        
        function load_tile_info(obj)
            exp_root_folder = obj.h_data_manager.fp_experiment(obj.experiment_group, ...
                obj.experiment);
            layer_folder_list = dir(exp_root_folder);
            layer_folder_list = layer_folder_list([layer_folder_list.isdir]);
            match_idx = ~cellfun(@isempty, regexpi({layer_folder_list.name}, '[0-9]{5}', 'match'));
            layer_list = arrayfun(@(x) str2double(x.name), layer_folder_list(match_idx));
            current_layer = obj.layer;
            
            num_layers = numel(layer_list);
            wb_hdl = waitbar(0, 'Loading...', 'Name', sprintf('Loading %s tiles', ...
                obj.name));
            for i = 1 : num_layers            
                obj.load_tiles_in_layer(layer_list(i));
                waitbar(i / num_layers, wb_hdl, sprintf('Loading layer %d', ...
                    layer_list(i)));
            end            
            wb_hdl.delete();
            obj.layer = current_layer;
        end
    end
    
    methods(Static)           
        function loaded_data = load_tile_data(tile_info, info_name, channel_id)
            if isempty(tile_info)
                loaded_data = [];
                return;
            end
            if nargin < 3
                channel_id = tile_info(1).channel;
            end
            num_tiles = numel(tile_info);
            num_channel = numel(channel_id);
            loaded_data = cell(num_channel, num_tiles);
            for i = 1 : num_tiles
                loaded_data(:, i) = tile_info(i).load_data(info_name, channel_id);
            end            
        end
        
        function new_tile_info = select_duplicated_tiles(tile_info, priority)
            % Warning: This function also sorts the tiles using their grid
            % indices. This is a side effect that can lead to bug.
           if nargin < 2
               priority = 'latest';
           end
           if isempty(tile_info)
               new_tile_info = [];
           else
               tile_grid_ind = [tile_info.grid_ind];
               tile_acq_count = [tile_info.acq_count];
               % Here the grid ind has been sorted
               [ind_list, bin_val] = Grid.bin_data_to_idx_list(tile_grid_ind);
               
               if numel(tile_grid_ind) == numel(bin_val)
                   % Sort the array - to have consistent output order
%                    ind_list = cat(1, ind_list{:});
%                    new_tile_info = tile_info(ind_list);
                     new_tile_info = tile_info;
               else
                   num_bin = numel(ind_list);
                   new_tile_info = cell(num_bin, 1);
                   for iter_bin = num_bin : -1 : 1
                       tmp_idx = ind_list{iter_bin};
                       if numel(tmp_idx) > 1
                           tmp_act_count_list = tile_acq_count(tmp_idx);
                           switch priority
                               case 'latest'
                                   [~, tmp_selected_idx] = max(tmp_act_count_list);
                               case 'oldest'
                                   [~, tmp_selected_idx] = min(tmp_act_count_list);
                               otherwise
                                   % specified by number
                                   tmp_act_count_list = tmp_act_count_list(tmp_act_count_list <= priority);
                                   [~, tmp_selected_idx] = max(tmp_act_count_list);
                           end
                           tmp_idx = tmp_idx(tmp_selected_idx);
                       end
                       new_tile_info{iter_bin} = tile_info(tmp_idx);
                   end
                   new_tile_info = cat(1, new_tile_info{:});
               end
           end
        end
        
        function [local_im, varargout] = get_stitched_tile_mip(tile_info, channel_id)
            validateattributes(channel_id, {'numeric'}, {'positive', 'scalar', 'integer'});
            loaded_data = WBIMTileManager.load_tile_data(tile_info, 'mip', channel_id);                 
            local_im = Grid.stitch_tiles(cat(1, tile_info.tile_mmxx_pxl), loaded_data, 'mip');          
            if nargout > 1
                local_tile_info = WBIMTileMetadata.combine_tile_info(tile_info);
                varargout{1} = local_tile_info;
            end
        end
        
        function [local_im_yx_um, varargout] = get_stitched_tile_mip_yx_um(tile_info, channel_id)
            [local_im_yx_um, local_tile_info] = WBIMTileManager.get_stitched_tile_mip(tile_info, channel_id);
            if any(size(local_im_yx_um) ~= local_tile_info.region_bbox_ll_um)
                local_im_yx_um = imresize(local_im_yx_um, local_tile_info.region_bbox_ll_um);
            end
            if nargout > 1
                varargout{1} = local_tile_info;
            end
        end
        
        function [stitched_mip, local_tile_info] = get_stitched_step_mip_1_ch(test_tiles, ch_idx, ...
                sec_idx_list, options)
            arguments
                test_tiles
                ch_idx
                sec_idx_list (1, :) double = []
                options.nanTo0Q (1,1) logical = false;
            end
            test_tiles = WBIMTileManager.select_duplicated_tiles(test_tiles);
            combined_stat = WBIMTileManager.load_and_combine_tile_stat(test_tiles, ch_idx);
            assert(isvector(combined_stat.step_mip_size), 'Step MIP stacks have different size');
            num_sections = combined_stat.step_mip_size(3);
            if isempty(sec_idx_list)
                sec_idx_list = 1 : num_sections;
            else
                num_sections = numel(sec_idx_list);
            end
            
            local_tile_info = WBIMTileMetadata.combine_tile_info(test_tiles);
            local_tile_info.step_mip_pixel_yxz_um = combined_stat.step_mip_pixel_size_um;
            
            local_tile_info.stage_z_um_done = unique(local_tile_info.stage_xyz_um_done(:, 3));
            assert(isscalar(local_tile_info.stage_z_um_done), 'Tiles were acquired at different stage z position');
            
            local_tile_info.step_mip_z_r_um = (sec_idx_list - 0.5) * ...
                local_tile_info.step_mip_pixel_yxz_um(3);
            
            local_tile_info.abs_fp_z_um_done = local_tile_info.stage_z_um_done + ...
                local_tile_info.piezo_z0_r_um;
            local_tile_info.step_mip_abs_fp_z_um = local_tile_info.abs_fp_z_um_done + ...
                local_tile_info.step_mip_z_r_um;
            local_tile_info.step_mip_stage_z_um = local_tile_info.stage_z_um_done + ...
                local_tile_info.step_mip_z_r_um;

            mip_mm_pxl_list = round(local_tile_info.tile_mmll_um(:, 1:2) ./ local_tile_info.step_mip_pixel_yxz_um(1:2));
            mip_ll_pxl_list = round(local_tile_info.tile_mmll_um(:, 3:4) ./ local_tile_info.step_mip_pixel_yxz_um(1:2));
            mip_mmxx_um_list = cat(2, mip_mm_pxl_list, mip_mm_pxl_list + mip_ll_pxl_list - 1);
            local_tile_info.step_mip_mmxx_pxl = mip_mmxx_um_list;
            stitched_mip = cell(num_sections, 1);
            for i = 1 : num_sections
                sec_idx = sec_idx_list(i);
                tmp_mip = cellfun(@(x) x(:, :, sec_idx), combined_stat.step_mip, 'UniformOutput', false);
                stitched_mip{i} = Grid.stitch_tiles(mip_mmxx_um_list, tmp_mip, 'mip');
            end
            stitched_mip = cat(3, stitched_mip{:});
            if options.nanTo0Q
                stitched_mip(isnan(stitched_mip)) = 0;
            end
            local_tile_info.step_mip_size = size(stitched_mip);
        end
        
        function [smip_cell, local_tile_info] = get_stitched_step_mip(tiles_in_layer, ch_idx_list, ...
                sec_idx_list, options)
            arguments
                tiles_in_layer (:, 1) WBIMTileMetadata
                ch_idx_list (1, :) double {mustBePositive}
                sec_idx_list (1, :) double = []
                options.nanTo0Q (1,1) logical = false;
            end
            tiles_in_layer = WBIMTileManager.select_duplicated_tiles(tiles_in_layer);
            num_ch = numel(ch_idx_list);
            max_ch_id = max(ch_idx_list);
            smip_cell = cell([max_ch_id, 1]);
            for i = 1 : num_ch
                tmp_ch = ch_idx_list(i);
                if i == 1
                    [smip_cell{tmp_ch}, local_tile_info] = WBIMTileManager.get_stitched_step_mip_1_ch(...
                        tiles_in_layer, tmp_ch, sec_idx_list, 'nanTo0Q', options.nanTo0Q);
                else
                    [smip_cell{tmp_ch}, ~] = WBIMTileManager.get_stitched_step_mip_1_ch(...
                        tiles_in_layer, tmp_ch, sec_idx_list, 'nanTo0Q', options.nanTo0Q);
                end                
            end            
        end
        
        
        function [combined_stat] = load_and_combine_tile_stat(test_tiles, ch_idx)
            ch_name = sprintf('CH%d', ch_idx);
            clearvars tile_stat
            for i = numel(test_tiles) : -1 : 1
                try
                    tile_stat(i) = test_tiles(i).load_stat(ch_idx).(ch_name);
                catch ME
                    fprintf('Failed to load tile stat. Skip. Error message: %s\n', ...
                        getReport(ME, 'extended'));
                end
            end
            combined_stat = struct;
            combined_stat.line_shift = cat(1, tile_stat.line_shift);
            combined_stat.step_mip = {tile_stat.step_mip};
            combined_stat.step_mip_pixel_size_um = cat(1, tile_stat.step_mip_pixel_yxz_um);
            unique_size_um = unique(combined_stat.step_mip_pixel_size_um, 'rows');
            if isvector(unique_size_um)
                combined_stat.step_mip_pixel_size_um = unique_size_um;
            else
                warning('Step MIP stacks have different voxel size');
            end
            if isfield(tile_stat, 'mip')
                combined_stat.mip = cat(1, tile_stat.mip);
            end
            if isfield(tile_stat, 'z_stat')
                combined_stat.z_stat = cat(1, tile_stat.z_stat);
            end
            
            combined_stat.step_mip_size = cellfun(@size, combined_stat.step_mip, ...
                'UniformOutput', false);
            combined_stat.step_mip_size = cat(1, combined_stat.step_mip_size{:});
            mip_size = unique(combined_stat.step_mip_size, 'rows');
            if isvector(mip_size)
                combined_stat.step_mip_size = mip_size;
            end
        end
        
        function tile_info = sort_tile_by_grid_ind(tile_info)
            arguments
                tile_info (:, 1) WBIMTileMetadata
            end
            tile_sub = cat(2, cat(1, tile_info.grid_sub), cat(1, tile_info.layer));
            max_sub = max(tile_sub, [], 1);
            tile_ind_1 = sub2ind(max_sub, tile_sub(:, 1), tile_sub(:, 2), tile_sub(:, 3));
            [~, tile_list_idx] = sort(tile_ind_1, 'ascend');
            tile_info = tile_info(tile_list_idx);
        end
        
    end
    %% Visualization GUI
    methods
        function fig_hdl = visualize_local_tile_grid(obj)
            
        end
        
        function [fig_hdl, varargout] = visualize_imaged_tile_mip_sample_yx_um(...
                obj, options)
            arguments
                obj
                options.channel (1,:) double = [1];
                options.layer (1,1) double = obj.layer;                 
            end
            num_ch = numel(options.channel);
            layer_tile = obj.get_latest_tiles_in_layer(options.layer);   
            local_tile_info = obj.grid.get_tiles_local_bbox(...
                [layer_tile.grid_ind]);
            im_bbox_mm_um = local_tile_info.g_bbox_mmll(1:2) .* ...
                local_tile_info.pixel_size_um;
            im_bbox_xx_um = local_tile_info.g_bbox_mmxx(3:4) .* ...
                local_tile_info.pixel_size_um;
            im_ll_um = im_bbox_xx_um - im_bbox_mm_um + 1;
            num_labels = 5;
            tick_x_val = linspace(1, im_ll_um(2), num_labels);
            tick_x_label = arrayfun(@(x) num2str(round(x), '%d'), ...
                tick_x_val + im_bbox_mm_um(2), 'UniformOutput', false);
            tick_y_val = linspace(1, im_ll_um(1), num_labels);
            tick_y_label = arrayfun(@(x) num2str(round(x), '%d'), ...
                tick_y_val + im_bbox_mm_um(1), 'UniformOutput', false);
            
            fig_hdl = figure;
            fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* 2;
            for i = 1 : num_ch
                if num_ch == 1
                    ax_hdl = axes(fig_hdl);
                else
                    ax_hdl = subplot(1, num_ch, i);
                end
                obj.local_mip_yx_um{options.channel(i)} = obj.get_latest_imaged_mip_yx_um_in_layer(...
                    options.channel(i), options.layer);                
                imagesc(ax_hdl, obj.local_mip_yx_um{options.channel(i)});
                ax_hdl.DataAspectRatio = [1,1,1];
                ax_hdl.XLabel.String = 'CRS / Sample X (\mum)';
                ax_hdl.YLabel.String = 'Galvo / Sample Y (\mum)';
%                 ax_hdl.Title.String = sprintf('Layer Z position: %d \\mum options.channel %d', ...
%                     round(obj.layer_z_um), options.channel(i));
                ax_hdl.XTick = tick_x_val;
                ax_hdl.XTickLabel = tick_x_label;
                ax_hdl.YTick = tick_y_val;
                ax_hdl.YTickLabel = tick_y_label; 
                cbar_hdl = colorbar();
                cbar_hdl.Label.String = 'Intensity';
            end

            if nargout == 2
                varargout{1} = local_tile_info;
            end
        end
        
        function [fig_hdl, varargout] = visualize_imaged_tile_merged_mip_sample_yx_um(obj, ...
                options)
            arguments
                obj
                options.channel (1, :) double = [1,2];
                options.layer (1,1) double = obj.layer
                options.displayQ (1,1) logical = true;
                options.im_filepath char = [];
            end         
            layer_tile = obj.get_latest_tiles_in_layer(options.layer);            
            local_tile_info = obj.grid.get_tiles_local_bbox(...
                [layer_tile.grid_ind]);
            im_bbox_mm_um = local_tile_info.g_bbox_mmll(1:2) .* ...
                local_tile_info.pixel_size_um;
            im_bbox_xx_um = local_tile_info.g_bbox_mmxx(3:4) .* ...
                local_tile_info.pixel_size_um;
            im_ll_um = im_bbox_xx_um - im_bbox_mm_um + 1;
            num_labels = 5;
            tick_x_val = linspace(1, im_ll_um(2), num_labels);
            tick_x_label = arrayfun(@(x) num2str(round(x), '%d'), ...
                tick_x_val + im_bbox_mm_um(2), 'UniformOutput', false);
            tick_y_val = linspace(1, im_ll_um(1), num_labels);
            tick_y_label = arrayfun(@(x) num2str(round(x), '%d'), ...
                tick_y_val + im_bbox_mm_um(1), 'UniformOutput', false);
            
            num_ch = numel(options.channel);   
            mip_cell = cell(1, num_ch);
            for i = 1 : num_ch
                tmp_im = obj.get_latest_imaged_mip_yx_um_in_layer(...
                    options.channel(i), options.layer);   
                tmp_im(isnan(tmp_im)) = 0;
                mip_cell{i} = uint16(tmp_im);
            end
            if options.channel(1) == WBIMChannelName.SHG 
                if num_ch == 2
                    % Set the skull channel to be channel 2...
                    mip_cell{3} = mip_cell{1};
                    mip_cell(1) = [];
                elseif num_ch ==3
                    tmp = mip_cell{2};
                    mip_cell{2} = mip_cell{1};
                    mip_cell{1} = tmp;
                end            
            end
            
            vis_im = fun_merge_image_stacks(mip_cell, 'method', 'GraySkull');
            
            if options.displayQ
                fig_hdl = figure;
            else
                fig_hdl = figure('Visible', 'off');
            end
            fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* 2;
            ax_hdl = axes(fig_hdl);
            image(ax_hdl, vis_im);
            ax_hdl.DataAspectRatio = [1,1,1];
            ax_hdl.XLabel.String = 'CRS / Sample X (\mum)';
            ax_hdl.YLabel.String = 'Galvo / Sample Y (\mum)';
            %                 ax_hdl.Title.String = sprintf('Layer Z position: %d \\mum Channel %d', ...
            %                     round(obj.layer_z_um), channel(i));
            ax_hdl.XTick = tick_x_val;
            ax_hdl.XTickLabel = tick_x_label;
            ax_hdl.YTick = tick_y_val;
            ax_hdl.YTickLabel = tick_y_label;
            ax_hdl.Title.String = sprintf('Layer %d', options.layer);

            if nargout == 2
                varargout{1} = local_tile_info;
            end
            
            if ~isempty(options.im_filepath)
                fun_print_image(fig_hdl, options.im_filepath);
            end            
            if ~options.displayQ
               delete(fig_hdl); 
            end
        end        
        
        function fig_hdl = visualize_candidate_map(obj)
            [fig_hdl, ax_hdl] = obj.visualize_grid(obj.candidate_map, ...
                'Candidate Grid Map');
        end
        
        function fig_hdl = visualize_imaged_tile_map(obj)
            [fig_hdl, ax_hdl] = obj.visualize_grid(obj.imaged_map, ...
                'Imaged Grid Map');
        end
    end
    %% Static visualization method
    methods(Static, Hidden)
        function [fig_hdl, ax_hdl] = visualize_grid(grid_map, grid_name)
            fig_hdl = figure;
            ax_hdl = axes(fig_hdl);
            imagesc(ax_hdl, grid_map);
            ax_hdl.XLabel.String = 'Grid X direction';
            ax_hdl.YLabel.String = 'Grid Y direction';
            ax_hdl.DataAspectRatio = [1,1,1];
            ax_hdl.Title.String = grid_name;
            colorbar(ax_hdl);
        end
    end
    
    methods(Static)
        function [fig_hdl, ax_hdl] = visualize_tile_mip_sample_yx_um(tile_info, ...
                ch_id)
            num_ch = numel(ch_id);
            local_tile_info = WBIMTileMetadata.combine_tile_info(tile_info);
            im_bbox_mm_um = local_tile_info.region_bbox_mm_um(1:2);
            im_ll_um = local_tile_info.region_bbox_ll_um;           
            
            num_labels = 5;
            tick_x_val = linspace(1, im_ll_um(2), num_labels);
            tick_x_label = arrayfun(@(x) num2str(round(x), '%d'), ...
                tick_x_val + im_bbox_mm_um(2), 'UniformOutput', false);
            tick_y_val = linspace(1, im_ll_um(1), num_labels);
            tick_y_label = arrayfun(@(x) num2str(round(x), '%d'), ...
                tick_y_val + im_bbox_mm_um(1), 'UniformOutput', false);
            
            fig_hdl = figure;
            fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* 2;
            for i = 1 : num_ch
                if num_ch == 1
                    ax_hdl = axes(fig_hdl);
                else
                    ax_hdl = subplot(1, num_ch, i);
                end
                local_im = WBIMTileManager.get_stitched_tile_mip_yx_um(...
                    tile_info, ch_id(i));
                imagesc(ax_hdl, local_im);
                ax_hdl.DataAspectRatio = [1,1,1];
                ax_hdl.XLabel.String = 'CRS / Sample X (\mum)';
                ax_hdl.YLabel.String = 'Galvo / Sample Y (\mum)';
%                 ax_hdl.Title.String = sprintf('Layer Z position: %d \\mum options.channel %d', ...
%                     round(obj.layer_z_um), options.channel(i));
                ax_hdl.XTick = tick_x_val;
                ax_hdl.XTickLabel = tick_x_label;
                ax_hdl.YTick = tick_y_val;
                ax_hdl.YTickLabel = tick_y_label; 
                cbar_hdl = colorbar();
                cbar_hdl.Label.String = 'Intensity';
            end
        end        
    end
end