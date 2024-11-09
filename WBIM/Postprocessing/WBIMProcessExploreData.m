classdef WBIMProcessExploreData < WBIMAcqPostProcessing
%     Online processing of raw SI tiff files
    %% To do list:
    % 1. *** Add chunk size
    % 2. * Use low-level hdf5
    properties(Hidden, Constant)
        COMPONENT_NAME = 'WBIMProcessExploreData';
    end
    properties
        tile_info WBIMTileMetadata
        mip (1, :) cell
        stat (1, :) cell
    end
    %%
    methods
        function obj = WBIMProcessExploreData(tile_info)
            obj@WBIMAcqPostProcessing();
            persistent DM
            if isempty(DM)
                DM = WBIMFileManager();
            end
            if ~isa(tile_info, 'WBIMTileMetadata')
                if isfile(tile_info)
                    tile_info = DM.load_data(tile_info);
                    assert(isa(tile_info, 'WBIMTileMetadata'));
                else
                    error('The input is neither a WBIMTileMetadata nor a filepath.');
                end
            end
            obj.tile_info = tile_info;
            obj.log_filepath = fullfile(obj.tile_info.save_root_folder, ...
                obj.tile_info.fprr_pplog);
            obj.err_log_filepath = fullfile(obj.tile_info.save_root_folder, ...
                obj.tile_info.fprr_pperr);
            obj.init_logger();
            obj.process();
        end
        
        function exit_code = process(obj)
            obj.h_logger.write("MESSAGE", "Start processing");
            raw_im_stack = obj.parse_SI_tiff();
            if isempty(raw_im_stack)
                obj.h_err_logger.write("ERROR", "Unable to parse the SI Tiff file");
                return
            end
            exit_code_1 = obj.organize_SI_tiff(raw_im_stack);
            exit_code_2 = obj.move_SI_tiff();
            obj.compute_mip_stat();
            exit_code = exit_code_1 || exit_code_2;
        end
        
        function [raw_im_stack] = parse_SI_tiff(obj)
            persistent DM
            if isempty(DM)
                DM = WBIMFileManager();
            end
            raw_im_hdl = [];
            try
                si_source_filepath = obj.tile_info.SI_filepath;
                if ~isfile(si_source_filepath)
                    si_source_filepath = fullfile(obj.tile_info.save_root_folder, ...
                        obj.tile_info.fprr_raw);
                    if ~isfile(si_source_filepath)
                        obj.h_err_logger.write("ERROR", "The SI raw tiff stack does not exist");
                        raw_im_stack = [];
                        return
                    end
                end
                si_metadata_filepath = fullfile(obj.tile_info.save_root_folder, ...
                    obj.tile_info.fprr_si_info);
                obj.h_logger.write("DEBUG", sprintf("Trying to read SI file %s", ...
                    si_source_filepath));
                raw_im_hdl = ScanImageTiffReader.ScanImageTiffReader(si_source_filepath);
                raw_im_stack = raw_im_hdl.data();
                % Flip y
                %                 raw_im_stack = flip(raw_im_stack, 1);
                %                 raw_im_stack = flip(raw_im_stack, 2);
                % Transpose
                % For swapping the X and Y axis - allows tiles to be
                % stitched in the grid coordinate directly
                
                %
                im_info = raw_im_hdl.metadata();
                raw_im_hdl.close();
                obj.write_si_metadata(si_metadata_filepath, im_info);
                
%                 raw_im_stack = WBIMFileManager.load_single_tiff(si_source_filepath);
                raw_im_stack = permute(raw_im_stack, [2,1,3]);
                obj.h_logger.write("MESSAGE", "Finish loading SI tiff");
            catch ME
                if ~isempty(raw_im_hdl) && isvalid(raw_im_hdl)
                    raw_im_hdl.close();
                else
                    obj.h_err_logger.write("ERROR", "raw_im_hdl is not valid!");
                end
                obj.h_err_logger.write("ERROR", getReport(ME, 'extended', 'hyperlinks', 'off'));
                rethrow(ME);
            end
        end
        
        function exit_code = organize_SI_tiff(obj, raw_im_stack)
            exit_code = -1;
            try
                raw_im_type = class(raw_im_stack);
                num_channel = numel(obj.tile_info.channel);
                if mod(size(raw_im_stack, 3), num_channel)
                    obj.h_err_logger.write("ERROR", "Missing data");
                    return
                end
                h5_fp = fullfile(obj.tile_info.save_root_folder, ...
                    obj.tile_info.fprr_tile);
                if isfile(h5_fp)
                    obj.h_logger.write("WARNING", "Target h5 file already exist");
                    h5_info = h5info(h5_fp);
                else
                    h5_info = [];
                end
                obj.mip = cell(1, num_channel);
                for iter_ch = 1 : num_channel
                    % Slicing
                    tmp_im = raw_im_stack(:, :, iter_ch : num_channel : end);
                    if strcmpi(WBIMConfig.SI_CONVERTED_IM_TYPE, raw_im_type)
                        tmp_im = cast(tmp_im, WBIMConfig.SI_CONVERTED_IM_TYPE);
                    end                    
                    tmp_im_size = size(tmp_im);
                    % Write H5
                    tmp_h5_ds = sprintf('%s/%s', obj.tile_info.h5_dataset{iter_ch}, ...
                        'raw');
                    
                    if ~isempty(h5_info) && ismember(tmp_h5_ds(2:end), [h5_info.Datasets.Name], 'rows')
                        warning('%s%s alreadly exist. Overwrite existing dataset',...
                            h5_fp, tmp_h5_ds);
                    else
                        h5create(h5_fp, tmp_h5_ds, tmp_im_size, 'Datatype', raw_im_type);
                    end
                    h5write(h5_fp, tmp_h5_ds, tmp_im);
                    obj.h_logger.write("MESSAGE",  sprintf('Finish writing data in channel %d',...
                        obj.tile_info.channel(iter_ch)));
                    % Write maximum intensity projection
                    tmp_im_mip = max(tmp_im, [], 3);
                    obj.mip{iter_ch} = tmp_im_mip;
                    tmp_mip_fp = fullfile(obj.tile_info.save_root_folder, obj.tile_info.fprr_mip{iter_ch});
                    
                    WBIMFileManager.write_tiff_stack(tmp_im_mip, tmp_mip_fp, 'grayscale');
                    obj.h_logger.write("MESSAGE", sprintf('Finish writing maximum intensity projection of channel %d', ...
                        obj.tile_info.channel(iter_ch)));
                    exit_code = 0;
                end
                
            catch ME
                obj.h_err_logger.write("ERROR", getReport(ME, 'extended', 'hyperlinks', 'off'));
                rethrow(ME);
            end
        end
        
        function exit_code = move_SI_tiff(obj)
            exit_code = WBIMFileManager.move_file_local_w(...
                obj.tile_info.SI_filepath, ...
                fullfile(obj.tile_info.save_root_folder, obj.tile_info.fprr_raw), true);
        end
        
        function compute_mip_stat(obj)
            obj.stat = cell(size(obj.mip));
            for i = 1 : numel(obj.mip)
                tmp_im = obj.mip{i};
                obj.stat{i} = obj.analyze_single_tile_mip(tmp_im);
            end
        end
    end
    %%
    methods(Static)
        function yes = is_background_tile_Q(mip_stat)
            if iscell(mip_stat)
                mip_stat = mip_stat{1};
            end
            % To be refined. 
            yes = (mip_stat.mean_mft < 1000 || mip_stat.std_rf < 0.5);            
        end
        
        function merged_im = merge_channels(im_cell, merge_method)
            if nargin < 2
                merge_method = 'color';
            end
            is_valid_channel_Q = ~cellfun(@isempty, im_cell);
            if ~all(is_valid_channel_Q)
                im_cell = im_cell(is_valid_channel_Q);
            end
            num_channels = numel(im_cell);
            if isnumeric(im_cell) 
                merged_im = im_cell;
            elseif num_channels == 1 && isnumeric(im_cell{1})
                merged_im = im_cell{1};
            else
                switch merge_method
                    case 'mip'
                        for i = 1 : num_channels
                            if i == 1
                                merged_im = rescale(im_cell{1});
                            else
                                merge_im = max(merge_im, rescale(im_cell{i}));
                            end
                        end
                    case 'color'
                        im_size = size(im_cell{1});
                        % Reprocessing 
                        for i = 1 : num_channels
                            im_cell{i} = fun_stretch_contrast(im_cell{i}, 0.01, 0.99, ...
                                'uint8');
                        end
                        switch num_channels
                            case 2
                                merged_im = zeros([im_size, 3], 'uint8');
                                merged_im(:, :, 1) = im_cell{1};
                                merged_im(:, :, 3) = im_cell{2};
                                merged_im(:, :, 2) = (merged_im(:, :, 1) + ...
                                    merged_im(:, :, 3)) ./ 2;
                            case 3
                                merged_im = cat(3, im_cell{:});
                            otherwise
                                error('To be implemented');
                        end
                end                               
            end
        end
    end
end