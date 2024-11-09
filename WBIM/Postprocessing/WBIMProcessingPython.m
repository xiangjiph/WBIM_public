classdef WBIMProcessingPython < handle
   
    properties
%        python_exe_path char
%        tiff_parser_path char
    end
    
    properties(Constant)
        cmd_template = 'start /min cmd /c "%s"';
    end
    
    methods
        function obj = WBIMProcessingPython()
%             obj.file_manager = WBIMFileManager();
                        
        end        
    end
    
    methods(Static)
        function exit_code = processing_scan_data(json_file_str)
            if ~ischar(json_file_str)
               assert(isa(json_file_str, 'WBIMTileMetadata'), 'Unrecognized input')
               json_file_str = json_file_str.fp_info_json;
            end
            file_manager = WBIMFileManager();
            python_str = file_manager.PYTHON_EXE_PATH;
            script_str = fullfile(file_manager.PYTHON_SCRIPT_PATH, 'SITiffParser.py');
            cmd_str = sprintf('start /min cmd /c "ping -n 2 127.0.0.1 > NUL & %s %s -f %s &"', python_str, script_str, ...
                json_file_str);
            exit_code = system(cmd_str);
        end        
    end
end