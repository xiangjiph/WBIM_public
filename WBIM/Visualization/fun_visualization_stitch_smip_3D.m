function [smip_im_c, varargout] = fun_visualization_stitch_smip_3D(process_tiles, opt)

arguments
   process_tiles cell
   opt.layer_range (:,1)  = []
   opt.channel (:, 1) = []    
   opt.zero_last_frame_Q (1, 1) logical = false;
   opt.zero_first_frame_Q (1, 1) logical = false;
end

valid_layer_Q = ~cellfun(@isempty, process_tiles);
process_tiles = process_tiles(valid_layer_Q);
process_tiles_s = cellfun(@(x) WBIMTileManager.select_duplicated_tiles(x), process_tiles, 'UniformOutput', false);
num_layer = numel(process_tiles_s);
ch_list = process_tiles_s{1}(1).channel;
num_ch = numel(ch_list);
layer_data = cell(num_layer, 1);
% Somehow the following block cannot be parallelized. Support of OOP is
% limited? 
for i = 1 : num_layer
    tmp_tile = process_tiles_s{i};
    tmp_data = struct;
    tmp_data.smip = cell(1, num_ch);
    for ch = 1 : num_ch
        [tmp_data.smip{ch}, tmp_data.local_tile_info] = WBIMTileManager.get_stitched_step_mip_1_ch(...
            tmp_tile, ch_list(ch));
    end
    arrayfun(@(x) x.clear_buffer(), tmp_tile);
    layer_data{i} = tmp_data;
    fprintf('Finish loading data in layer %d\n', i);
end
layer_data = cat(1, layer_data{:});
%% Compute the bounding box
smip_pixel_yxz_um = layer_data(1).local_tile_info.step_mip_pixel_yxz_um;
stack_size = layer_data(1).local_tile_info.stack_size;
vol_bbox_xy_mm_um = arrayfun(@(x) x.local_tile_info.region_bbox_mm_um, layer_data, 'UniformOutput', false);
vol_bbox_xy_mm_um = cat(1, vol_bbox_xy_mm_um{:});
vol_bbox_z_mm_um = arrayfun(@(x) x.local_tile_info.stage_z_um_done , layer_data);
vol_bbox_xy_xx_um = arrayfun(@(x) x.local_tile_info.region_bbox_xx_um, layer_data, 'UniformOutput', false);
vol_bbox_xy_xx_um = cat(1, vol_bbox_xy_xx_um{:});

vol_bbox_mmxx_um = cat(2, vol_bbox_xy_mm_um, vol_bbox_z_mm_um, vol_bbox_xy_xx_um, vol_bbox_z_mm_um + stack_size(3) - 1);
vol_bbox_mmxx_sp = round(vol_bbox_mmxx_um ./ [smip_pixel_yxz_um, smip_pixel_yxz_um]);
space_mm_sp = min(vol_bbox_mmxx_sp(:, 1:3), [], 1);
space_xx_sp = max(vol_bbox_mmxx_sp(:, 4:6), [], 1);
space_ll_sp = space_xx_sp - space_mm_sp + 1;
vol_bbox_mmxx_sp_l = vol_bbox_mmxx_sp - [space_mm_sp, space_mm_sp] + 1;
%%
stitch_im_type = 'uint16';
num_s_ch = numel(layer_data(1).smip);
smip_im_c = cell(num_s_ch, 1);
for i_c = 1 : num_s_ch
    tmp_smip_im = zeros(space_ll_sp, stitch_im_type);
    for i_l = 1 : num_layer
        tmp_im = layer_data(i_l).smip{i_c};
        if opt.zero_last_frame_Q
            tmp_im(:, :, end) = 0;
        end        
        if opt.zero_first_frame_Q
            tmp_im(:, :, 1) = 0;
        end
        tmp_im_size = size(tmp_im);
        tmp_bbox_mm = vol_bbox_mmxx_sp_l(i_l, 1:3);
        tmp_bbox_xx = tmp_bbox_mm + tmp_im_size - 1;
        tmp_smip_im(tmp_bbox_mm(1) : tmp_bbox_xx(1), tmp_bbox_mm(2) : tmp_bbox_xx(2), tmp_bbox_mm(3) : tmp_bbox_xx(3)) = ...
            max(tmp_smip_im(tmp_bbox_mm(1) : tmp_bbox_xx(1), tmp_bbox_mm(2) : tmp_bbox_xx(2), tmp_bbox_mm(3) : tmp_bbox_xx(3)), uint16(tmp_im));        
        fprintf('Finish adding layer %d of channel %d\n', i_l, i_c);
    end        
    smip_im_c{i_c} = tmp_smip_im;
end

if nargout > 1
    varargout{1} = layer_data;
end
end