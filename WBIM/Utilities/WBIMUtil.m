classdef WBIMUtil

    methods(Static)
        function t = get_log_timestamp(log_str_array, ref_t)
            arguments
                log_str_array
                ref_t = []
            end
            log_str_array = arrayfun(@(x) strsplit(x, {' ', '\t'}), log_str_array, 'UniformOutput', false);
            log_str_array = cellfun(@(x) datetime([x{1}, ' ', x{2}]), log_str_array, 'UniformOutput', false);
            log_str_array = cat(1, log_str_array{:});
            if ~isempty(ref_t)
                t = log_str_array - ref_t;
            end
        end
    end

end