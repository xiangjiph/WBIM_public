classdef WBIMROIDetection < handle
    properties
        layer_tile_info (:, 1) WBIMTileMetadata
        num_tile (1,1) double
        mip_cell cell
        im
        im_bbox_mmxx_pxl
        pixel_size_um
        im_bbox_mm_um (1, 2) double
        im_bbox_xx_um (1, 2) double
    end
    
    %%
    methods
        function obj = WBIMROIDetection(varargin)
            if numel(varargin) == 1
               input1 = varargin{1};
               if isa(input1, 'WBIMTileMetadata')
                   tile_info_vec = input1;
               elseif isfolder(input1)
                  layer_acq_folder = input1; 
                  tile_info_vec = obj.load_tile_info_from_acq_folder(layer_acq_folder);
               end
            end
            obj.layer_tile_info = tile_info_vec;
            obj.init();
        end
        
        function obj = init(obj)
            obj.num_tile = numel(obj.layer_tile_info);
            pxl_size = cat(1, obj.layer_tile_info.pixel_size_um);
            pxl_size = unique(pxl_size, 'rows');
            assert(isvector(pxl_size), 'Inconsistent pixel size');
            obj.pixel_size_um = pxl_size;      
        end
        %%        
        function stitch_mip(obj)
            bbox_mmxx_list = cat(1, obj.layer_tile_info.tile_mmxx_pxl);
            overall_mm_pxl = min(bbox_mmxx_list(:, 1:2), [], 1);
            overall_xx_pxl = max(bbox_mmxx_list(:, 3:4), [], 1);
            overall_ll_pxl = overall_xx_pxl - overall_mm_pxl + 1;
            local_mmxx_list = bsxfun(@minus, bbox_mmxx_list, [overall_mm_pxl, overall_mm_pxl] - 1);
            
            obj.im = nan(overall_ll_pxl);
            for i_tile = 1 : obj.num_tile
                tile_mmxx = local_mmxx_list(i_tile, :);
                obj.im(tile_mmxx(1) : tile_mmxx(3), tile_mmxx(2) : tile_mmxx(4)) = ...
                    obj.mip_cell{i_tile};                
            end
            obj.im_bbox_mmxx_pxl = [overall_mm_pxl, overall_xx_pxl];
            obj.im_bbox_mm_um = obj.im_bbox_mmxx_pxl(1:2) .* obj.pixel_size_um(1:2);
            obj.im_bbox_xx_um = obj.im_bbox_mmxx_pxl(3:4) .* obj.pixel_size_um(1:2);
        end

    end
end