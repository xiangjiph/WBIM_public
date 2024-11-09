classdef WBIMProcessScanData < WBIMAcqPostProcessing
    % Online processing of SI raw tiff files
    %% To do list:
    % 1. *** Add chunk size
    % 2. * Use low-level hdf5
    properties(Hidden, Constant)
        COMPONENT_NAME = 'WBIMProcessScanData';
    end
    properties
        tile_info WBIMTileMetadata
        mip (1, :) cell
    end
    %%
    methods
        function obj = WBIMProcessScanData(tile_info)
            obj@WBIMAcqPostProcessing();
            DM = WBIMFileManager();
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
            obj.clean_up();
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
            exit_code = exit_code_1 || exit_code_2;
        end
        
        function [raw_im_stack] = parse_SI_tiff(obj)
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
                
                raw_im_hdl = ScanImageTiffReader.ScanImageTiffReader(si_source_filepath);
                raw_im_stack = raw_im_hdl.data();
                % Transpose the image stack
                raw_im_stack = permute(raw_im_stack, [2,1,3]);
                im_info = raw_im_hdl.metadata();
                raw_im_hdl.close();
                obj.write_si_metadata(si_metadata_filepath, im_info);
                obj.h_logger.write("MESSAGE", "Finish loading SI tiff");
            catch ME
                raw_im_hdl.close();
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
                    tmp_im_size = size(tmp_im);
                    % Write H5
                    tmp_h5_ds = obj.tile_info.h5_dataset{iter_ch};
                    
                    if ~isempty(h5_info) && ismember(tmp_h5_ds(2:end), cat(1, h5_info.Datasets.Name), 'rows')
                        obj.h_logger.write("WARNING", sprintf('%s%s alreadly exist. Overwrite existing dataset',...
                            h5_fp, tmp_h5_ds));
                    else
                        h5create(h5_fp, tmp_h5_ds, tmp_im_size, 'Datatype', raw_im_type);
                    end
                    h5write(h5_fp, tmp_h5_ds, tmp_im);
                    obj.h_logger.write("MESSAGE",  sprintf('Finish writing data in channel %d',...
                        obj.tile_info.channel(iter_ch)));
                    % Write maximum intensity projection
                    tmp_im_mip = max(tmp_im, [], 3);
                    tmp_mip_fp = fullfile(obj.tile_info.save_root_folder, obj.tile_info.fprr_mip{iter_ch});
                    
                    WBIMFileManager.write_tiff_stack(tmp_im_mip, tmp_mip_fp, 'grayscale');
                    obj.h_logger.write("MESSAGE", sprintf('Finish writing maximum intensity projection of channel %d', ...
                        obj.tile_info.channel(iter_ch)));
                    obj.mip{iter_ch} = tmp_im_mip;
                end
                exit_code = 1;
            catch ME
                obj.h_err_logger.write("ERROR", getReport(ME, 'extended', 'hyperlinks', 'off'));
                rethrow(ME);
            end
            %% Test
            %             im_written = imread(tmp_mip_fp);
            %             assert(all(im_written == tmp_im_mip, 'all'), 'Loaded data are different from the written data');
            %             im_written_1 = WBIMFileManager.load_single_tiff(tmp_mip_fp);
            %             assert(all(im_written_1 == uint16(tmp_im_mip), 'all'), 'Loaded data are different from the written data');
        end
        
        function exit_code = move_SI_tiff(obj)
            exit_code = WBIMFileManager.move_file_local_w(...
                obj.tile_info.SI_filepath, ...
                fullfile(obj.tile_info.save_root_folder, obj.tile_info.fprr_raw), true);
        end
        
        function clean_up(obj)
            delete(obj.h_logger);
            delete(obj.h_err_logger);
        end
    end
    %%
    methods(Static)
        function exit_code = write_si_metadata(fp, info)
            [folder, ~, ~] = fileparts(fp);
            if ~isfolder(folder)
                mkdir(folder);
            end
            f_hdl = fopen(fp, 'a+');
            fprintf(f_hdl, '%s', info);
            exit_code = fclose(f_hdl);
        end
        
        function yes = is_background_tile_Q(mip_stat)
            if iscell(mip_stat)
                mip_stat = mip_stat{1};
            end
            % To be refined.
            yes = (mip_stat.mean < 500);
        end
        
        function fig_hdl = vis_step_mip_and_mask(im, mask, sec)
            if nargin < 3
                sec = 1;
            end
            fig_hdl = figure;
            tiledlayout(1, 3);
            ax_1 = nexttile;
            imagesc(ax_1, im(:, :, sec));
            ax_1.DataAspectRatio = [1,1,1];
            ax_1.Title.String = 'Image';
            ax_2 = nexttile;
            imagesc(ax_2, mask(:, :, sec));
            ax_2.Title.String = 'Mask';
            ax_2.DataAspectRatio = [1,1,1];
            ax_3 = nexttile;
            imshowpair(rescale(im(:, :, sec)), mask(:, :, sec));
            ax_3.Title.String = 'Overlaid';
        end
    end
end