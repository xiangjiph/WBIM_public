classdef WBIMTCPServer < handle
    % To do 
    % 1. Set hostname and port number when luanching the TCPClient 
    
    
    %%
    properties(Constant)
        HOST = "localhost"
        PORT = 4000
        
        LOG_FILE_THRESHOLD = mlog.Level.DEBUG
        LOG_TERMINAL_THRESHOLD = mlog.Level.DEBUG
        LOG_EVENT_THRESHOLD = mlog.Level.WARNING        
        
        COMPONENT_NAME = "WBIMTCPServer"
        DELIMITERS = "::"
    end
    
    properties
        % Online processing 
        process_tile_running (:, 1) WBIMTileMetadata = WBIMTileMetadata.empty()
        process_tile_done (:, 1) WBIMTileMetadata = WBIMTileMetadata.empty()
        process_tile_failed (:, 1) WBIMTileMetadata = WBIMTileMetadata.empty()
        
        % Sync tiles 
        sync_folder_queue (2, :) cell
        sync_folder_failed (2, :) cell
        running_sync_Q (1, 1) logical = false;
    end
    
    properties(Constant)
        CLIENT_SHUT_DOWN_TOKEN = "WBIMTCPClient OFF"
    end
    
    properties
        server tcpserver.internal.TCPServer
        received_data string
    end
    
    properties(Transient, SetAccess=protected)
       h_logger mlog.Logger        
    end
    %% Events
    events
        ERowShiftDetected % Fires when the python processor detects row shift in the image
    end    
    %% Construction     
    methods
        function obj = WBIMTCPServer(log_filepath)
            arguments
               log_filepath = "" 
            end
            obj.init_logger(log_filepath);
            try
                obj.server = tcpserver(obj.HOST, obj.PORT);
                obj.server.Timeout = 20;
                % 1 will get triggered too many times... There is no a good
                % number for this application. Just use the terminator. 
%                 configureCallback(obj.server, "byte", 1, @obj.BytesAvailableFcn);
                configureCallback(obj.server, "terminator", @obj.BytesAvailableFcn);
                %                 configureTerminator(obj.server, "LF");
                obj.server.ConnectionChangedFcn = @obj.ConnectionChangedFun;
            catch ME
                obj.server.delete();
                rethrow(ME);
            end
        end
        
        function delete(obj)
            obj.shut_down_client();
            pause(2);
            obj.server.delete();
            obj.h_logger.delete();
        end        
    end
    %% Outward communication 
    methods
        function write(obj, data)
            arguments
                obj
                data {char, string}                
            end
            if obj.server.Connected
                obj.server.writeline(data);
            else
                obj.h_logger.write("MESSAGE", "No client is connected to the server");
            end            
        end
        
        function send_data(obj, data)
            if isstring(data) || ischar(data)
                data_header = "String";
            elseif isnumeric(data)
                data_header = "Numeric";
                data = mat2str(data);
            end
            data = strjoin([obj.COMPONENT_NAME, data_header, data], obj.DELIMITERS);
            obj.write(data);
        end
        
        function shut_down_client(obj)
            data = strjoin([obj.COMPONENT_NAME, obj.CLIENT_SHUT_DOWN_TOKEN], obj.DELIMITERS);
            obj.write(data);
        end
        
        function exit_code = launch_client(obj)
           % Launch python TCP client, specify address, ip, delimiters,
           % etc. Run in the background? 
           try
               if ~obj.server.Connected
                   file_manager = WBIMFileManager();
                   python_str = file_manager.PYTHON_EXE_PATH;
                   script_str = fullfile(file_manager.PYTHON_SCRIPT_PATH, 'WBIMTCPClient.py');
                   if ispc
                       cmd_str = sprintf('start /min cmd /c "%s %s&"', ...
                           python_str, script_str);
                   elseif isunix
                       cmd_str = sprintf('cd %s; %s; %s %s&', file_manager.PYTHON_SCRIPT_PATH,...
                           file_manager.CONDA_PATH, python_str, script_str);
                   end
                   exit_code = system(cmd_str);
                   obj.h_logger.write("MESSAGE", "Finish launching WBIMTCPClient");
               else
                   obj.h_logger.write("MESSAGE", "WBIMTCPClient is already running");
               end
           catch ME
               obj.h_logger.write("ERROR", sprintf("Fail to launch the TCP client: %s", ...
                   getReport(ME, 'extended', 'hyperlinks', 'off')));
           end
        end
    end
    %%
    methods
        function set.running_sync_Q(obj, val)
            obj.running_sync_Q = val;
            if val
                obj.h_logger.write("MESSAGE", "Start syncing folders to the server");
            else
                obj.h_logger.write("MESSAGE", "Stop syncing folders to the server");
            end
            obj.sync_next_folder();            
        end
        
    end    
    %% Client functions
    methods
        function submitted_Q = submit_tile_for_processing(obj, tile_info)
            arguments
                obj WBIMTCPServer
                tile_info (1,1) WBIMTileMetadata
            end
            addedQ = obj.add_one_tile(tile_info, WBIMTCPServerTask.ProcessTile, 'check_duplicate_Q', true);
            submitted_Q = false;
            if addedQ
                tile_info_fp = tile_info.fp_info_json;
                cmd = strjoin([obj.COMPONENT_NAME, "ProcessTile", tile_info_fp], ...
                obj.DELIMITERS);
                obj.write(cmd);
                tile_info.processing_state = WBIMProcessingState.Running;
                submitted_Q = true;
            end
        end
        
        function set_client_parameters(obj, para, val)
            assert(isstring(para) || ischar(para));
            if isnumeric(val)
                val = mat2str(val);
            end
            cmd = strjoin([obj.COMPONENT_NAME, "Setting", para, val], ...
                obj.DELIMITERS);
            obj.write(cmd);
        end
        

        function add_sync_folders(obj, folder_pairs, uniqueQ)
            arguments
                obj (1,1) WBIMTCPServer
                % The first row is the soruce folder, the second row is the
                % target folder
                folder_pairs (2, :) cell 
                uniqueQ (1,1) logical = true;
            end
            % Compare the source folder 
            if uniqueQ && ~isempty(obj.sync_folder_queue)
                num_new_pair = size(folder_pairs, 2);
                current_source_queue = obj.sync_folder_queue(1, :);
                is_new_Q = true(1, num_new_pair);
                for i = 1 : num_new_pair
                    in_list_Q = contains(current_source_queue, folder_pairs{1, i});
                    if any(in_list_Q)
                        target_in_list_Q = contains(obj.sync_folder_queue(2, in_list_Q), ...
                            folder_pairs{2, i});
                        if any(target_in_list_Q)
                            is_new_Q(i) = false;
                        end
                    end
                end
                folder_pairs = folder_pairs(:, is_new_Q);
            end            
            obj.sync_folder_queue = cat(2, obj.sync_folder_queue, folder_pairs);
        end

        
        function submittedQ = sync_single_folder(obj, source_folder, target_folder)
            arguments
                obj WBIMTCPServer
                source_folder 
                target_folder
            end
            msg = jsonencode(struct('source_folder', source_folder, ...
                'target_folder', target_folder));
            cmd = strjoin([obj.COMPONENT_NAME, string(WBIMTCPServerTask.SyncTile), ...
                msg], obj.DELIMITERS);
            obj.h_logger.write("MESSAGE", sprintf('Syncing folder %s to %s', ...
                source_folder, target_folder));
            obj.write(cmd);
            submittedQ = true;
        end
        
        function sync_next_folder(obj)
            if obj.running_sync_Q
                if ~isempty(obj.sync_folder_queue)
                    next_folder = obj.sync_folder_queue(:, 1);
                    obj.sync_single_folder(next_folder{1}, next_folder{2});
                else
                    obj.h_logger.write("MESSAGE", "Reach the end of sync folder queue");
                end
            end
        end
    end
    
    %% Job schaduling 
    methods
        function addedQ = add_one_tile(obj, tile_info, task, options)
            arguments
                obj WBIMTCPServer
                tile_info WBIMTileMetadata
                task (1,1) WBIMTCPServerTask = WBIMTCPServerTask.ProcessTile
                options.check_duplicate_Q (1,1) logical = true                
            end
            switch task
                case WBIMTCPServerTask.ProcessTile
                    addedQ = ~options.check_duplicate_Q || isempty(obj.process_tile_running) || ...
                        ~any(obj.tile_is_in_list(tile_info));
                    if addedQ
                        obj.process_tile_running(end+1) = tile_info;
                    else
                        obj.h_logger.write("DEBUG", "The new tile is alreadly in the processing list");
                    end                    
            end
        end
        
        function in_list_Q = tile_is_in_list(obj, new_tile)
            if ~isempty(obj.process_tile_running)
                if isa(new_tile, 'WBIMTileMetadata')
                    in_list_Q = (obj.process_tile_running == new_tile);
                else
                    assert(ischar(new_tile) || isstring(new_tile))
                    tile_fp_list = {obj.process_tile_running.fp_info_json};
                    in_list_Q = contains(tile_fp_list, new_tile);
                end
            else
                in_list_Q = false;
            end
        end     
        
        function event_callback_tile_done(obj, rs_msg)
            % The following properties are hardcoded in python 
            successQ = rs_msg.doneQ;
            tile_fp = rs_msg.tile_info_fp;
            % 
            in_list_Q = obj.tile_is_in_list(tile_fp);
            assert(nnz(in_list_Q) == 1, 'One and only one tile should match the filepath');
            done_tile = obj.process_tile_running(in_list_Q);
            if successQ
                done_tile.processing_state = WBIMProcessingState.Done;
                obj.process_tile_running = obj.process_tile_running(~in_list_Q);
                obj.h_logger.write("INFO", sprintf("Successfully processed %s", tile_fp));
                
                obj.rs_add_new_estimation(rs_msg.num_row_shift);
            else
                % Failed
                obj.process_tile_failed(end+1) = done_tile;
                obj.h_logger.write("INFO", sprintf("Failed to process %s", tile_fp));
            end
            % Do not log the completed tile at the moment
%             obj.process_tile_done(end+1) = done_tile;
        end
        
        function event_callback_tile_failed(obj, tile_fp)
            in_list_Q = obj.tile_is_in_list(tile_fp);
            assert(nnz(in_list_Q) == 1, 'One and only one tile should match the filepath');
            failed_tile = obj.process_tile_running(in_list_Q);
            failed_tile.processing_state = WBIMProcessingState.Failed;
            obj.process_tile_running = obj.process_tile_running(~in_list_Q);
            obj.process_tile_failed(end+1) = failed_tile;
        end
        
        function resubmit_failed_tiles(obj)
            for i = numel(obj.process_tile_failed) : -1 : 1
                obj.submit_tile_for_processing(obj.process_tile_failed(i));
                obj.process_tile_failed(i) = [];
            end
        end
        
        function event_callback_sync_single_folder(obj, msg)
            if msg.return_code == 0
                obj.h_logger.write("MESSAGE", sprintf("Successfully sync folder %s", msg.local_folder));
            else
                obj.h_logger.write("MESSAGE", sprintf("Fail to sync folder %s\nError:\n%s\nCommand:\n%s",...
                    msg.local_folder, msg.error, msg.command));
            end            
            % Sync the next tile
            if strcmpi(obj.sync_folder_queue{1, 1}, msg.local_folder)
                if msg.return_code ~= 0
                    obj.sync_folder_failed(:, end+1) = obj.sync_folder_queue(:, 1);
                end
                obj.sync_folder_queue(:, 1) = [];
                obj.sync_next_folder();
            else
                obj.h_logger.write("MESSAGE", "The synced folder is not in the folder list");
            end
        end
        
    end
    

    %% Inward communication - parsing received data 
    methods
        function parse_received_data(obj, data)
            str_cell = strsplit(data, obj.DELIMITERS);
            if numel(str_cell) == 1
                obj.h_logger.write("DEBUG", sprintf("Receive unrecognized message: %s", ...
                    data));
            else
                assert(str_cell{1} == "WBIMTCPClient", "Unrecognized message head");
                switch char(str_cell{2})
                    case "INFO"
                        obj.h_logger.write("MESSAGE", sprintf("WBIMTCPClient: %s", ...
                            str_cell{3}));
                    case "EVENT"
                        switch char(str_cell{3})
                            case "ProcessTile"
                                obj.event_callback_tile_done(jsondecode(str_cell{4}));
                            case "SyncFolder"
                                obj.event_callback_sync_single_folder(jsondecode(str_cell{4}));
                            otherwise
                                obj.h_logger.write("INFO", sprintf("Unrecognize message: %s", str_cell{4}));
                        end
                    case "DATA"
                        
                    otherwise
                        obj.h_logger.write("DEBUG", sprintf("Receive unrecognized message: %s", ...
                            data));
                end
            end            
        end
    end    
    %% Row shift estimation - to be improve 
    properties        
        max_row_shift (1,1) double = 30
        min_pause_row_shift (1,1) double = 0
    end
    
    methods        
        function rs_add_new_estimation(obj, new_val)
            if new_val ~= 0
                if abs(new_val) > obj.max_row_shift 
                    obj.h_logger.write("DEBUG", sprintf("The new estiamte of row shift is %d - too big", new_val));
                else
                    obj.h_logger.write("INFO", "Row shift artefact detected");
                    if abs(new_val) > obj.min_pause_row_shift
                        notify(obj, "ERowShiftDetected");
                    end
                end
            end
        end        
        
    end    
    %% Callback functions
    methods
        function ConnectionChangedFun(obj, src, connInfo)
            arguments
                obj (1,1) WBIMTCPServer
                src 
                connInfo
            end
            if src.Connected
                obj.h_logger.write("MESSAGE", "Client has connected");
            else
                obj.h_logger.write("MESSAGE", "Client has disconnected");
            end
        end
        
        function data = BytesAvailableFcn(obj, arg1, arg2)
            arguments
                obj (1,1) WBIMTCPServer
                arg1
                arg2
            end
            % If the data contains multiple `\n`, this function will be
            % triggered multiple times. It seems that all the data will be
            % read in the first time? Not sure if need to handle multiple
            % read or not. 
            if obj.server.NumBytesAvailable
                data = obj.server.read(obj.server.NumBytesAvailable, 'char');
                string_cells = strsplit(data, '\n');
                string_cells = string_cells(~cellfun(@isempty, string_cells));
                for i = 1 : numel(string_cells)
                    obj.parse_received_data(string_cells{i});
%                     obj.h_logger.write("DEBUG", sprintf('New data received: %s', string_cells{i}));
                end
            end
            % Parsing data here 
        end
    end   
    %% Loggers    
    methods
        function init_logger(obj, log_filepath)
            arguments
                obj WBIMTCPServer
                log_filepath = "";
            end
            obj.h_logger = mlog.Logger(obj.COMPONENT_NAME, log_filepath);
            obj.h_logger.FileThreshold = obj.LOG_FILE_THRESHOLD;
            obj.h_logger.CommandWindowThreshold = obj.LOG_TERMINAL_THRESHOLD;
            obj.h_logger.MessageReceivedEventThreshold = obj.LOG_EVENT_THRESHOLD;
            obj.h_logger.write("DEBUG", "Finish initializing logger");
        end
        
        function update_logger_filepath(obj, log_filepath)
            arguments
                obj WBIMTCPServer
                log_filepath {char, string}
            end
            obj.h_logger.LogFile = log_filepath;            
        end
    end
end