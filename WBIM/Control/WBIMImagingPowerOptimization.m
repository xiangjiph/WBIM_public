classdef WBIMImagingPowerOptimization < handle
    
    properties
        hSI scanimage.SI
        num_frame_per_z (1,1) double
        
        max_power_fraction (1,1) double = 0.8;
        max_opt_iteration (1,1) double = 100;
        
        current_iteration (1,1) double = 0
        
        target_max_int (1,1) double = double(intmax(WBIMConfig.SI_IM_TYPE)) * 0.95;
        target_max_int_fraction (1,1) double = 0.001;
        opt_channel (1,1) double = 1;
        
        reference_abs_z_um (1,1) double
        exp_length_um (1,1) double
    end
    properties(Dependent)
       current_power_fraction (1,1) double 
    end
    
    
    properties(Hidden)
        buffer_frame cell
        buffer_count (1,1) double = 0
        opt_data struct = struct
        h_OSP_listener = [];
    end
    
    properties(Access=private)
        PIEZO_RANGE_um = [-200, 200];
        init_z_r_um (1,1) double = nan
    end
    %%
    methods
        function obj = WBIMImagingPowerOptimization(hSI)
            obj.hSI = hSI;
        end
        
        function delete(obj)
            obj.delete_listeners(obj.h_OSP_listener);
        end
        
        function val = get.current_power_fraction(obj)
            val = obj.hSI.hBeams.powerFractions;
        end
        
%% SI related
        function set_piezo_park_z_r_um(obj, z_r_um)
            validateattributes(z_r_um, {'numeric'}, {'scalar','finite','real', ...
                '>=', obj.PIEZO_RANGE_um(1), '<=', obj.PIEZO_RANGE_um(2)});
            if obj.hSI.hFastZ.hFastZs{1}.parkPosition ~= z_r_um
                obj.hSI.hFastZ.hFastZs{1}.parkPosition = z_r_um;
                obj.hSI.hFastZ.hFastZs{1}.park();
            end
        end
        
        function move_piezo_to_z_r_um(obj, z_r_um)
            validateattributes(z_r_um, {'numeric'}, {'scalar','finite','real', ...
                '>=', obj.PIEZO_RANGE_um(1), '<=', obj.PIEZO_RANGE_um(2)});
            if obj.hSI.hFastZ.position ~= z_r_um
                obj.hSI.hFastZ.hFastZs{1}.move(z_r_um);
            end
        end
  
%% Single plane optimization
        function optimize_single_plane_power(obj, z_r_um)
            if nargin < 2
                z_r_um = obj.hSI.hFastZ.position;
            end
            obj.init_z_r_um = obj.hSI.hFastZ.position;            
            obj.num_frame_per_z = 10;
            obj.current_iteration = 1;
            obj.osp_init_buffer();
            obj.opt_data = struct;
            obj.opt_data.power_fraction = nan(obj.max_opt_iteration, 1);
            obj.opt_data.avg_saturation_frac = nan(obj.max_opt_iteration, 1);
            
            obj.h_OSP_listener = [addlistener(obj.hSI.hUserFunctions, 'frameAcquired', @obj.osp_single_frame_done), ...
                addlistener(obj.hSI.hUserFunctions, 'acqAbort', @obj.osp_clean_up)];
            % Use focus mode  
            obj.opt_data.power_fraction(obj.current_iteration) = obj.current_power_fraction;
            obj.hSI.hBeams.pzAdjust = scanimage.types.BeamAdjustTypes.None;
            obj.hSI.start('focus');
            obj.move_piezo_to_z_r_um(z_r_um);
        end
        
        function osp_single_frame_done(obj, varargin)
            obj.record_current_SI_image();
            
            if obj.buffer_count == obj.num_frame_per_z
                % Average over frame
                avg_im = mean(cat(3, obj.buffer_frame{:}), 3);
                saturation_fraction = nnz(avg_im >= obj.target_max_int) / ...
                    numel(avg_im);
                obj.opt_data.avg_saturation_frac(obj.current_iteration) = ...
                    saturation_fraction;
                % Processing
%                 saturation_fraction = cellfun(@(x) nnz(x >= obj.target_max_int) / ...
%                     numel(x), obj.buffer_frame);
%                 obj.opt_data.avg_saturation_frac(obj.current_iteration) = ...
%                     mean(saturation_fraction);
                obj.osp_update_power_fraction();
                obj.osp_init_buffer();                
            end            
        end
        
        function osp_init_buffer(obj)
            obj.buffer_frame = cell(obj.num_frame_per_z, 1);
            obj.buffer_count = 0;
        end
        
        function osp_clean_up(obj, varargin)
            obj.delete_listeners(obj.h_OSP_listener);            
            obj.move_piezo_to_z_r_um(obj.init_z_r_um);
        end
        
        function osp_update_power_fraction(obj)
            opt_step_size = 0.03;
            
            if obj.current_iteration <= obj.max_opt_iteration
                if obj.current_iteration == 1 || ...
                        prod(obj.target_max_int_fraction - ...
                        obj.opt_data.avg_saturation_frac([obj.current_iteration-1, obj.current_iteration])) > 0
                    next_power = obj.current_power_fraction + opt_step_size * ...
                        sign(obj.target_max_int_fraction - obj.opt_data.avg_saturation_frac(obj.current_iteration));
                else
                    next_power = (obj.opt_data.power_fraction(obj.current_iteration) + ...
                        obj.opt_data.power_fraction(obj.current_iteration-1)) / 2;
                end
                next_power = round(next_power, 4);                
                if next_power <= obj.max_power_fraction && next_power > 0
                    if all(obj.opt_data.power_fraction ~= next_power) 
                        obj.current_iteration = obj.current_iteration + 1;
                        obj.opt_data.power_fraction(obj.current_iteration) = ...
                            obj.current_power_fraction;
                    else
                        % Get back to an intensity previously searched 
                        [~, min_idx] = min(abs(obj.target_max_int_fraction - ...
                            obj.opt_data.avg_saturation_frac(1:obj.current_iteration)));
                        next_power = obj.opt_data.power_fraction(min_idx);
                        obj.hSI.abort();
                    end
                    obj.hSI.hBeams.powerFractions = next_power;
                    return;                    
                end  
            end     
            obj.hSI.abort();
        end
%% Single plane intensity scan 
        function intensity_scan_single_plane(obj, power_list)
            assert(all(power_list >= 0) && all(power_list <= 1), ...
                'power should be in [0, 1]');
            obj.num_frame_per_z = 10;
            obj.current_iteration = 1;
            obj.opt_data = struct;
            obj.opt_data.power_list = power_list;
            obj.opt_data.num_intensity = numel(power_list);
            obj.opt_data.power_idx = 1;
            obj.opt_data.num_iteration = obj.num_frame_per_z * ...
                obj.opt_data.num_intensity;
            obj.buffer_frame = cell(obj.num_frame_per_z, ...
                obj.opt_data.num_intensity);
            obj.buffer_count = 0;
            
            obj.h_OSP_listener = [addlistener(obj.hSI.hUserFunctions, 'frameAcquired',...
                @obj.issp_single_frame_done), ...
                addlistener(obj.hSI.hUserFunctions, 'acqAbort', ...
                @obj.issp_clean_up)];
            
            obj.hSI.hBeams.pzAdjust = scanimage.types.BeamAdjustTypes.None;
            obj.hSI.hBeams.powerFractions = power_list(1);
            obj.hSI.start('focus');            
        end
        
        function issp_single_frame_done(obj, varargin)
            obj.record_current_SI_image();
            obj.current_iteration = obj.current_iteration + 1;
            if obj.current_iteration < obj.opt_data.num_iteration
                if mod(obj.current_iteration, obj.num_frame_per_z) == 0
                    obj.opt_data.power_idx = obj.opt_data.power_idx + 1;
                    obj.hSI.hBeams.powerFractions = obj.opt_data.power_list(obj.opt_data.power_idx);
                end
            else
                obj.hSI.abort();
            end            
        end
        
        function issp_clean_up(obj, varargin)
            obj.delete_listeners(obj.h_OSP_listener);  
        end
%% Z-stack optimization 

        function optimize_volume_power(obj)
            
        end
        
    end
%% Power functions
    methods
        function set_power_z_function(obj, fun_handle)
            if isa(fun_handle, 'function_handle')
                fun_handle = {fun_handle};
            end
            % This will work if removing the validation step in Line 287 '
            % scanimage.util.validateFunctionHandle(v{i});' of +scanimage/+components/Beams
            obj.hSI.hBeams.pzFunction = fun_handle;
            obj.hSI.hBeams.pzAdjust = scanimage.types.BeamAdjustTypes.Function;
        end
        
        function powers = exp_power(obj, powers, z_abs_um, hBeam)
            % powers is the power at the reference z position 
            % z_abs_um is the vector of absolute z position (stage z +
            % piezo relative z), in micron
            powers = powers .* exp((z_abs_um - obj.reference_abs_z_um)...
                ./ obj.exp_length_um); 
        end
        
        function powers = rectified_exp_power(obj, powers, z_abs_um, hBeam)
            powers = powers .* exp(max(0, (z_abs_um - obj.reference_abs_z_um))...
                ./ obj.exp_length_um);
        end 
    end    
    %% Utilities
    methods
        function record_current_SI_image(obj)
            si_stripe_data = obj.hSI.hDisplay.lastStripeData;
            si_im_cell = obj.si_parse_stripe_img_data(si_stripe_data);
            obj.buffer_count = obj.buffer_count + 1;
            obj.buffer_frame{obj.buffer_count} = si_im_cell{obj.opt_channel};
        end
    end
    
    methods(Static)
        function delete_listeners(h_listener)
            if ~isempty(h_listener)
               num_listeners = numel(h_listener);
               for i = 1 : num_listeners
                   h_listener(i).delete();
               end
            end
        end
        
        function imgData = si_parse_stripe_img_data(stripeData)
            % ScanImage's implementation
            stripeChans = stripeData.channelNumbers;
            if ~isempty(stripeData.roiData)
                imageData = stripeData.roiData{1}.imageData;
            
                chansAvail = 1 : max(stripeData.channelNumbers);

                missingChans = find(ismember(chansAvail, stripeChans)==0);
                % Insert empties for missing channels.
                for i = 1 : numel(missingChans)
                    missingChan = missingChans(i);
                    imageData = {imageData{1:missingChan-1}, {[]}, imageData{missingChan:end}};
                end

                % Remove extra cell layer.
                imgData = cellfun(@(x) x{1}, imageData, 'UniformOutput', false);
            else
                imgData = [];
            end
        end
    end    
end