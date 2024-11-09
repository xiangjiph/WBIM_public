classdef InlineRefractometer < handle
% Class for controlling and logging data from an inline refractometer

    properties
        port_DEV %(1,1) internal.Serialport 
        LED_on_Q (1,1) logical % This should be read directly. 
        cam webcam
        im
        im_bg
        im_binarized
        im_size
        timerMeasure timer  
        logQ (1, 1) logical = true;
    end
    
    properties(SetAccess=protected, SetObservable)
        BRIX (1,1) double = nan
    end   
    
    properties(Dependent)
        LCD_on_Q (1,1)logical        
        RI (1,1) double
    end
    
    properties (Access=private)
        use_Arduino_Q = false; % True for Arduino, false for vDAQ
        reprocessQ = false
        edge_crop_mask
        dig_len {mustBePositive, mustBeInteger} 
        state {mustBeNonnegative}
    end
    % SI port related
    properties
        h_sidport_LCD_on
        h_sidport_LCD_state
        h_sidport_LED
        h_sidport_measure
        
        h_sidtask_LCD_on
        h_sidtask_read
    end
    
    properties(Hidden)
        hRIC % For using the logger.   
    end
    
    properties(Access=private)
        t_lcd_turned_on (1,1) uint64 = nan;
        t_measurement_sent (1,1) uint64 = nan;
        
        c_try_read (1,1) uint8 = 0;
    end

    properties (Access=private, Constant)
        % keys for 0-9
        % [top_left, bot_left, top, mid, bot, top_right, bot_right]
        dig_keys = [1 1 1 0 1 1 1; %0
                    0 0 0 0 0 1 1; %1
                    0 1 1 1 1 1 0; %2
                    0 0 1 1 1 1 1; %3
                    1 0 0 1 0 1 1; %4
                    1 0 1 1 1 0 1; %5
                    1 1 1 1 1 0 1; %6
                    0 0 1 0 0 1 1; %7
                    1 1 1 1 1 1 1; %8
                    1 0 1 1 1 1 1]; %9
        pulse_length = 0.5; %s
    end

%% Events
    events
        EMeasurementDone % Triggered when RI measurement is done
    end    
    
%% Listeners
    properties(Hidden)
        h_listener = [];        
    end
%%
    methods
        function obj = InlineRefractometer(dev_port)
            arguments
                dev_port = [];
            end
            [~, hostname] = system('hostname');
            hostname = hostname(1:end-1);
            if isempty(dev_port)
                if ~strcmpi(hostname, 'pia')
                    try
                        % Begin serial communications with Arduino
                        obj.port_DEV  = serialport('/dev/ttyACM0',9600);
                        obj.use_Arduino_Q = true;
                    catch ME
                        error("Arduino port not accessible or currently in use.")
                    end
                else
                    try 
                        obj.init_si_do_task();
                    catch ME
                        obj.clean_up_si_related_handles();
                        rethrow(ME);
                    end                                        
                end     
            else
                obj.port_DEV = dev_port;
                % Need to check if dev_port is a serial port object
                obj.use_Arduino_Q = true;
            end
            try
                obj.init_timers();
                obj.init_LED();
                obj.init_camera();
                obj.turn_on_LCD();
                obj.init_compute();
            catch ME
                obj.delete();
                obj.log_info("ERROR", getReport(ME, 'extended', ...
                    'hyperlinks', 'off'));
            end
        end
%% Initialization        
        function init_timers(obj)            
            obj.timerMeasure = timer(...
                'ExecutionMode','fixedSpacing','StartDelay',0, ...
                'Period',1, ...
                'TasksToExecute',1);
            obj.timerMeasure.TimerFcn = @(~,~)obj.stateFunction();
        end
        
        function init_si_do_task(obj)
            rs = dabs.resources.ResourceStore();
            obj.h_sidport_LCD_on = rs.filterByName(WBIMConfig.RI_LCD_ON_PORT_NAME);
            obj.h_sidport_measure = rs.filterByName(WBIMConfig.RI_BRIX_READ_PORT_NAME);
            obj.h_sidport_LED = rs.filterByName(WBIMConfig.RI_RED_LED_PORT_NAME);
            obj.h_sidport_LCD_state = rs.filterByName(WBIMConfig.RI_LCD_STATE_PORT_NAME);
            
            obj.h_sidtask_LCD_on = InlineRefractometer.si_get_do_pulse_task(...
                obj.h_sidport_LCD_on, obj.pulse_length);
            obj.h_sidtask_read = InlineRefractometer.si_get_do_pulse_task(...
                obj.h_sidport_measure, obj.pulse_length);
            obj.log_info("MESSAGE", "Finish initializing SI digital output tasks");
        end     
        
        function init_LED(obj)
            % LED should start turned off
            obj.LED_on_Q = false;
            % Turn on LED and make sure LCD is off
            obj.turn_on_LED();
            obj.turn_off_LCD();
            pause(1);
            obj.log_info("MESSAGE", "Finish initializing LED");
        end
        
        function init_camera(obj)
            obj.cam = webcam(1);
            pause(1);
            preview(obj.cam);
            pause(2)    % Give some time for camera to start
            obj.im_bg = snapshot(obj.cam);
        end
        
        function init_compute(obj)
            try
                pause(7); % Waiting for LCD to turn on.
                % Take images of dashes for length of digits calculation
                im_dashes = obj.im_bg - snapshot(obj.cam);
                obj.im_size = size(im_dashes);
                
                % Use dash image to initialize edge cropping mask
                obj.edge_crop_mask = zeros(obj.im_size(1:2));
                obj.edge_crop_mask(floor(0.1*obj.im_size(1):0.9*obj.im_size(1)),...
                    floor(0.025*obj.im_size(2):0.975*obj.im_size(2))) = 1;
                
                % Digitize dashed image
                im_dashes_gray = obj.binBRIXImage(im_dashes);
                
                % Determine where dashes are using longest connected component
                % after width dilation
                im_dashes_gray = im_dashes_gray.* ...
                    obj.isolateBRIXImage(im_dashes_gray,"longest");
                
                % Extract approximate width of digits
                [im_dash_cropped, dash_col_rising_idx, dash_col_falling_idx] ...
                    = obj.cropBRIXImage(im_dashes_gray);
                dash_horiz_locs = [1, dash_col_rising_idx(2 : end) - dash_col_rising_idx(1) + 1; ...
                    dash_col_falling_idx(1 : (end-1)) - dash_col_rising_idx(1) + 1, size(im_dash_cropped, 2)];
                obj.dig_len = floor(mean(dash_horiz_locs(2, :) - dash_horiz_locs(1, :)));
                obj.log_info("DEBUG", "Finish init_compute of InlineRefractometer.");
            catch ME
                obj.log_info("ERROR", sprintf("Fail in init_comput. Error message: %s", ...
                    getReport(ME, 'extended', 'hyperlinks', 'on')));
            end
        end
%%
        function block_mask = isolateBRIXImage(obj, im_bin, superlative)
            arguments
                obj
                im_bin logical
                superlative {mustBeMember(superlative,{'largest','longest'})} = "largest"
            end

            im_bin = im_bin .* obj.edge_crop_mask;
            block_im = imdilate(im_bin, ...
                strel("rectangle", [1, floor(0.3 * obj.im_size(2))]));
            conn_blocks = bwconncomp(block_im);

            if(superlative == "largest")
                block_sz = cellfun(@numel, conn_blocks.PixelIdxList);
                [~, dig_block_ind] = max(block_sz);
            else
                longest_ind = 0;
                longest_len = 0;
                for i = 1 : conn_blocks.NumObjects
                    [~, block_cols] = ind2sub(obj.im_size, conn_blocks.PixelIdxList{i});
                    temp_len = max(block_cols) - min(block_cols);
                    if longest_len < temp_len
                        longest_len = temp_len;
                        longest_ind = i;
                    end
                end
                dig_block_ind = longest_ind;
            end
            block_mask = zeros(obj.im_size);
            block_mask(conn_blocks.PixelIdxList{dig_block_ind}) = 1;
        end

        function measure_RI(obj, reprocess)
            arguments
                obj
                reprocess logical = false
            end
            
            obj.reprocessQ = reprocess;

            if ~reprocess
                if obj.use_Arduino_Q
                    if ~obj.LCD_on_Q
                        obj.turn_on_LCD();
                    end
                    % send serial command to read LCD
                    obj.read_LCD();
                    % take picture of LCD with camera
                    obj.im = obj.im_bg - snapshot(obj.cam);
                else
                    if obj.LCD_on_Q
                        obj.state = 1;
                        tasks = 2;
                        period = 2.5;
                    else
                        obj.state = 0;
                        tasks = 3;
                        period = 3.75;
                    end
                end
            else
                obj.state = 2;
                tasks = 1;
                period = 0.01;
            end
            
            obj.timerMeasure.Period = period;
            obj.timerMeasure.TasksToExecute = tasks;
            wait(obj.timerMeasure)
            start(obj.timerMeasure)
        end

        function turn_on_LED(obj)
            if ~obj.LED_on_Q
                if obj.use_Arduino_Q
                    obj.port_DEV.writeline("LED")
                    pause(2)  % Give Arduino time to turn on LED
                else
                    obj.h_sidport_LED.setValue(true);
                end
                obj.LED_on_Q = true;
                obj.log_info("DEBUG", "Turning on LED");
            end
        end

        function turn_off_LED(obj)
            if obj.LED_on_Q
                if obj.use_Arduino_Q
                    obj.port_DEV.writeline("LED")
                    pause(2)  % Give Arduino time to turn off LED
                else
                    obj.h_sidport_LED.setValue(false);
                end
                obj.LED_on_Q = false;
                obj.log_info("DEBUG", "Turning on LED");
            end
        end
        
        function turn_on_LCD(obj)
            if ~obj.LCD_on_Q
                if obj.use_Arduino_Q
                    obj.port_DEV.writeline("LCD")
                else    
                    obj.h_sidtask_LCD_on.start();
                    obj.t_lcd_turned_on = tic;
                    pause(0.1);
                    assert(obj.LCD_on_Q, 'LCD is still off');
                end
                obj.log_info("DEBUG", "Turning on LCD");
                % Is this necessary? 7 second seems very long 
                % It took about 4 seconds to turn on. 
%                 pause(5);  
            end
        end

        function turn_off_LCD(obj)
            if obj.LCD_on_Q
                if obj.use_Arduino_Q
                    obj.port_DEV.writeline("LCD")
                    pause(3);   % Give time for LCD to turn off
                else
                    obj.h_sidtask_LCD_on.start();                    
                end
                obj.log_info("DEBUG", "Turning off LCD");
            end
        end

        function read_LCD(obj)
            if obj.LCD_on_Q       
                % Wait for the screen to turn on 
                while toc(obj.t_lcd_turned_on) < 5
                    pause(0.5);
                end
                if obj.use_Arduino_Q
                    obj.port_DEV.writeline("READ")
                    pause(5) % Wait for reading to be made
                else
                    obj.h_sidtask_read.start();
                end
                obj.c_try_read = 0;
                obj.log_info("DEBUG", "Start reading LCD");
            else
                obj.log_info("WARNING", "LCD is off");
                % Turn on here? 
                obj.c_try_read = obj.c_try_read + 1;
                obj.log_info("MESSAGE", sprintf("Try to turn on the LCD. Trial %d", ...
                    obj.c_try_read));
                if obj.c_try_read < 2
                    obj.turn_on_LCD();
                    obj.read_LCD();
                else
                    error('Fail to read LCD');                    
                end
            end
        end

        function delete(obj)
            obj.cam.closePreview();
            obj.turn_off_LED();
            obj.turn_off_LCD();            
            stop([obj.timerMeasure]);
            obj.timerMeasure.delete;
            
            if ~isempty(obj.h_listener)
                delete(obj.h_listener);
            end            
            % Close communications to Serial port
            obj.port_DEV = [];
            obj.clean_up_si_related_handles();
        end
    end
%% SI related
    methods
        function clean_up_si_related_handles(obj)
            [obj.h_sidport_LCD_on, obj.h_sidport_LCD_state, obj.h_sidport_LED, ...
                obj.h_sidport_measure] = deal([]);
            obj.si_cleanup_pulse_task(obj.h_sidtask_LCD_on);
            obj.si_cleanup_pulse_task(obj.h_sidtask_read);
        end
    end
    
    methods(Static)       
        function h_si_task = si_get_do_pulse_task(si_port, pulse_s, clean_up_Q)
            arguments
                si_port
                pulse_s (1,1) double
                clean_up_Q (1,1) logical = false;
            end
            h_si_task = dabs.vidrio.ddi.DoTask(si_port.hDAQ, ...
                'pulse trains');
            try
                h_si_task.addChannel(si_port.name);
                h_si_task.sampleMode = 'finite';
                num_high = round(pulse_s * h_si_task.sampleRate);
                buffer = cat(1, ones(num_high, 1), 0);
                h_si_task.writeOutputBuffer(buffer);
                h_si_task.samplesPerTrigger = numel(buffer);
                if clean_up_Q
                    h_si_task.doneCallback = @InlineRefractometer.si_cleanup_pulse_task;
                end
            catch ME
                InlineRefractometer.si_cleanup_pulse_task(h_si_task);
                rethrow(ME);
            end
        end        
        
        function hDigitalTask = si_cleanup_pulse_task(hDigitalTask)
            if most.idioms.isValidObj(hDigitalTask)
%                 fprintf("Delete SI digital output task\n");
                hDigitalTask.stop;
                hDigitalTask.abort;
                hDigitalTask.delete;
                clear hDigitalTask
                hDigitalTask = [];
            end
        end
    end
%%  Dependent methods
    methods
        function val = get.LCD_on_Q(obj)
            if obj.use_Arduino_Q
                flush(obj.port_DEV,"output");
                obj.port_DEV.writeline("STATE")
                val = str2num(readline(obj.port_DEV));
            else
                val = obj.h_sidport_LCD_state.readValue();
            end
        end
        
        function val = get.RI(obj)
            % convert from BRIX to RI using formula taken from
            % http://www.fruitmanagement.com/brix.html (using inspect
            % element on the calculator)
            val = 1.3330229 + 0.00142117428 * obj.BRIX +...
                0.0000056370904 * (obj.BRIX ^ 2) + ...
                0.0000000154588009 * (obj.BRIX ^ 3);
        end
    end    
    %%
    methods
        function process_new_image(obj)
            if ~obj.reprocessQ
                obj.im = obj.im_bg - snapshot(obj.cam);
            end
            % initial round of binarizing and processing
            im_gray = obj.binBRIXImage(obj.im) .* obj.edge_crop_mask;
            im_gray = im_gray .* obj.isolateBRIXImage(im_gray, "largest");
            [im_cropped,~,~] = obj.cropBRIXImage(im_gray);
            
            % correct for affine rotation with hough transform
            im_cropped = obj.houghRotateImage(im_cropped);
            [im_cropped, col_rising_idx, col_falling_idx] = obj.cropBRIXImage(im_cropped);
            obj.im_binarized = im_cropped;
            
            % locations of BRIX digits
            num_col = size(im_cropped, 2);
            digit_horiz_locs = [1, col_rising_idx(2 : end) - col_rising_idx(1) + 1; ...
                col_falling_idx(1 : (end-1)) - col_rising_idx(1) + 1, num_col];
            dig_horiz_blocks = digit_horiz_locs(2, :) - digit_horiz_locs(1, :);
            
            % initialize array for storing BRIX digits
            dig_val = nan(size(dig_horiz_blocks));
            
            % Artificial Intelligence: rules-based expert system
            for i = 1 : numel(dig_horiz_blocks)
                dig_i = im_cropped(:, digit_horiz_locs(1, i) : digit_horiz_locs(2, i));
                tmp_dig_val = obj.single_digit_image_to_digit(dig_i);
                if ~isempty(tmp_dig_val)
                    dig_val(i) = tmp_dig_val;
                end
            end            
            % remove erronenous small dots
            dig_val = dig_val(isfinite(dig_val));            
            % coalesce results into single BRIX value
            dig_val = dig_val - 1;
            
            obj.BRIX = sum(dig_val .* flip(10.^((1 : numel(dig_val))-2)));
            obj.reprocessQ = false;
        end
        
        function dig_val = single_digit_image_to_digit(obj, dig_i)
            h0 = height(dig_i);
            dig_i = obj.cropBRIXImage(dig_i);
            ih = height(dig_i);
            iw = width(dig_i);
            if (iw < 0.75 * obj.dig_len)
                dig_val = 2; %We subtract by 1 later
                if (ih < 0.5 * h0)
                    dig_val = [];
                end
            else                
                %dig_vec = [top_left, bot_left, top, mid, bot, top_right, bot_right]
                dig_vec = [mode(dig_i(floor(0.15 * ih) : floor(0.5 * ih), 1 : floor(0.35 * iw)), 'all'), ...
                    mode(dig_i(floor(0.5 * ih : end), 1 : floor(0.35 * iw)), 'all'), ...
                    mode(dig_i(1 : floor(0.15 * ih), :),'all'), ...
                    mode(dig_i(floor(0.5 * ih + (0.5 * (-0.15 * ih : 0.15 * ih))), :), 'all'), ...
                    mode(dig_i(floor((1 - 0.15) * ih) : end,:),'all'), ...
                    mode(dig_i(floor(0.15 * ih) : floor(0.5 * (ih)), ...
                    floor((1 - 0.35) * iw) : end), 'all'), ...
                    mode(dig_i(floor(0.5 * ih : floor((1 - 0.15) * ih)), ...
                    floor((1 - 0.35) * iw) : end),'all')];
                [~, dig_val] = min(sum((obj.dig_keys - repmat(dig_vec, [10,1])) .^2, 2));
            end
        end        
        
    end
    
    methods(Static)
        function im_bin = binBRIXImage(im, threshold)
            % Binarize the RGB image, return 2D binary image
            arguments
                im {mustBeNumeric}
                threshold {mustBeMember(threshold,{'global','adaptive'})} = "adaptive"
            end

            im_bin = im2gray(im); %Convert to grayscale

            % Threshold according to either lowest quantiles or adaptive
            if (threshold=="global")
                im_quantiles = quantile(im(:),0:0.05:1); %image intensity quantiles
                im_quantile_idx = find(im_quantiles);
                im_thresh = 0.5*(sum(im_quantiles(im_quantile_idx(1:2))));
                im_bin = im_bin<im_thresh; %Threshold image to binarize
            else
                im_bin = ~imbinarize(255 - im_bin,"adaptive", ...
                    "ForegroundPolarity","bright", ...
                    "Sensitivity",0.7); %Threshold image to binarize
            end
            
            im_bin = bwareaopen(im_bin, 4000); %Remove decimal point
            
            persistent se_disk
            if isempty(se_disk)
                se_disk = strel('disk', 20);
            end
            im_bin = imclose(im_bin, se_disk); %Fill gaps between seven-segement model
        end

        function [im_cropped, col_rising_idx, col_falling_idx] = cropBRIXImage(im_bin)
            arguments
                im_bin logical
            end
            % Determining the bounding box of the digit and crop the image.
            row_idx = diff(sum(im_bin, 2) > 0);
            row_rising_idx = find(row_idx == 1); %rising edge non-zero value
            row_falling_idx = find(row_idx== -1) - 1; %falling edge non-zero value
            if isempty(row_rising_idx)
                row_slicing_idx_1 = 1;
            else
                row_slicing_idx_1 = row_rising_idx(1);
            end
            
            if isempty(row_falling_idx)
                row_slicing_idx_2 = size(im_bin, 1);
            else
                if numel(row_falling_idx) > 1
                    keyboard
                end
                row_slicing_idx_2 = row_falling_idx(end);
            end
            digit_rows = row_slicing_idx_1 : row_slicing_idx_2;

            im_cropped = im_bin(digit_rows, :);

            col_idx = diff(sum(im_cropped, 1) > 0);
            col_rising_idx = find(col_idx == 1); %rising edge non-zero value
            col_falling_idx = find(col_idx == -1)-1; %falling edge non-zero value
            if isempty(col_rising_idx)
                col_slicing_idx_1 = 1;
            else
                col_slicing_idx_1 = col_rising_idx(1);
            end
            
            if isempty(col_falling_idx)
                col_slicing_idx_2 = size(im_cropped, 2);
            else
                col_slicing_idx_2 = col_falling_idx(end);
            end
            digit_cols = col_slicing_idx_1 : col_slicing_idx_2;
            im_cropped = im_bin(digit_rows, digit_cols); % cropped image
        end
        
        function im_warped = houghRotateImage(im_bin)
            % TODO: the projection angle is fixed. Can calculate the right
            % value once with visual inspection and use it later. 
            arguments
                im_bin logical
            end
            im_edges = edge(im_bin);
            [H,T,R] = hough(im_edges);
            P  = houghpeaks(H,20,'threshold',ceil(0.3*max(H(:))));
            shY = -mean(P(logical((P(:,2)<10).*(P(:,2)>-10)),2));
            shX = mean(P(logical((P(:,2)<100).*(P(:,2)>80)),2))-90;
            tform = affine2d([1 tand(shY) 0; tand(shX) 1 0; 0 0 1]);
            im_warped = imwarp(im_bin,tform);
        end
        
    end

    methods(Hidden, Access = protected)
        function stateFunction(obj)
            new_state = 3;
            switch obj.state
                case 0
                    obj.turn_on_LCD();
                    new_state = 1;
                case 1
                    obj.read_LCD();
                    new_state = 2;
                case 2
                    % Get the image of the new measurement on the LCD
                    obj.process_new_image();                    
                    obj.log_info("MESSAGE", sprintf('RI: %.4f\tBRIX: %.1f', obj.RI, obj.BRIX));
                case 3
            end
            obj.state = new_state;
        end
    end
    %% Utility
    methods
        function log_info(obj, msg_type, msg)
            if obj.logQ
                msg_type = upper(msg_type);
                if ~isempty(obj.hRIC) && isvalid(obj.hRIC)...
                        && ~isempty(obj.hRIC.hLogger) && isvalid(obj.hRIC.hLogger)
                    obj.hRIC.log_info(msg_type, msg);
                else
                    fprintf('InlineRefractometer %s: %s\n', msg_type, msg);
                end
            end
        end
    end
end