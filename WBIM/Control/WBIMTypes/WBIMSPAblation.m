classdef WBIMSPAblation <  matlab.mixin.Copyable & dynamicprops
    % WBIM Setting parameter - Albation 
    properties
        material (1, :) WBIMAblationMaterial
        with_diffuser_Q (1,:) logical = false;
        site_scale_factor_fast (1,:) double
        site_scale_factor_slow (1,:) double
        peak_fluence_J_cm2 (1,:) double        
    end
    
    properties(SetAccess=private)
        abl_z_offset_r_um (1,1) double = 0
        abl_range_r_um (1,:) double
        z_axis_step_um (1,:) double
        abl_z_r_um (1, :) double
    end
        
    methods
        function obj = WBIMSPAblation(options)
            arguments
                options.material (1, :) {mustBePositive, mustBeInteger}
                options.with_diffuser_Q (1,:) logical = false;
                options.site_scale_factor_fast (1,:) double
                options.site_scale_factor_slow (1,:) double
                options.peak_fluence_J_cm2 (1,:) double
                
                options.abl_range_r_um (1,:) double = []
                options.z_axis_step_um (1, :) double = []
                options.abl_z_r_um (1,:) double = [];
                % To compensate the surface movement due to flow pressure 
                options.abl_z_offset_r_um (1,1) double = 0 
            end
            obj.set_abl_z_r_um('abl_z_offset_r_um', options.abl_z_offset_r_um, ...
                'abl_z_r_um', options.abl_z_r_um, 'abl_range_r_um', ...
                options.abl_range_r_um, 'z_axis_step_um', options.z_axis_step_um);
            
            cp_fn = {'material', 'with_diffuser_Q', 'site_scale_factor_fast', ...
                'site_scale_factor_slow', 'peak_fluence_J_cm2', 'abl_z_offset_r_um'};
            num_field = numel(cp_fn);
            for i = 1 : num_field
                tmp_fn = cp_fn{i};
                obj.(tmp_fn) = options.(tmp_fn);
            end
        end
        
        function obj = set_abl_z_r_um(obj, options)
           arguments
               obj (1,1) WBIMSPAblation
               options.abl_range_r_um (1,:) double = obj.abl_range_r_um
               options.z_axis_step_um (1,1) double = obj.z_axis_step_um
               options.abl_z_r_um (1, :) double = [];
               options.abl_z_offset_r_um (1,1) double = obj.abl_z_offset_r_um;
           end
           obj.abl_z_offset_r_um = options.abl_z_offset_r_um;
           if ~isempty(options.abl_z_r_um)
               obj.abl_z_r_um = options.abl_z_r_um + options.abl_z_offset_r_um;
               if ~isempty(options.z_axis_step_um) || ~isempty(options.abl_range_r_um)
                   warning('Given abl_z_r_um, z_axis_step_um and abl_range_r_um are ignored');
               end
           else
               obj.z_axis_step_um = options.z_axis_step_um;
               obj.abl_range_r_um = options.abl_range_r_um;
               obj.abl_z_r_um = (obj.abl_range_r_um(1) : obj.z_axis_step_um : obj.abl_range_r_um(2))...
                   + obj.abl_z_offset_r_um;
           end                       
        end
        
        function val = get_abl_z_r_um(obj, options)
            arguments
                obj (1, 1) WBIMSPAblation
                options.abl_range_r_um (1, :) double = [];
            end
            if isempty(options.abl_range_r_um)
                val = obj.abl_z_r_um;
            else
                assert(numel(options.abl_range_r_um) == 2);
                val = options.abl_range_r_um(1) : obj.z_axis_step_um : ...
                    options.abl_range_r_um(2);
                val = val + obj.abl_z_offset_r_um;
                if iscolumnvector(val)
                    val = val.';
                end
            end
            
        end
        
        
        
        function testQ = has_the_same_ablation_parameter_Q(obj, obj2)
            arguments
                obj (1,1) WBIMSPAblation
                obj2 (:, 1) WBIMSPAblation
            end
            test_prop_name = {'with_diffuser_Q', 'site_scale_factor_fast', ...
                'site_scale_factor_slow', 'peak_fluence_J_cm2'};
            num_obj = numel(obj2);
            testQ = true([num_obj, 1]);
            for i = 1 : numel(test_prop_name)
                tmp_n = test_prop_name{i};
                tmp_val_1 = obj.(tmp_n);
                tmp_val_2 = cat(1, obj2.(tmp_n));
                if all(size(tmp_val_1) == size(tmp_val_2))
                    testQ = testQ & (tmp_val_1 == tmp_val_2);
                else
                    testQ = false;
                end
                if ~testQ
                    return;
                end
            end
        end
        
        function testQ = eq(obj, obj2)
            arguments
                obj (1,1) WBIMSPAblation
                obj2 (:,1) WBIMSPAblation
            end            
            testQ = obj.has_the_same_ablation_parameter_Q(obj2);
            for i = 1 : numel(testQ)
               testQ(i) = testQ(i) && (numel(obj.abl_z_r_um) == numel(obj2(i).abl_z_r_um)) ...
                   && all(obj.abl_z_r_um == obj2(i).abl_z_r_um); 
            end
                
        end
    end 
    
    methods(Static)
        function objs_bin = bin_by_ablation_parameter(objs)
            arguments
                objs (:, 1) WBIMSPAblation
            end
            num_objs = numel(objs);
            objs_label = zeros(num_objs, 1);
            label = 0;
            for i = 1 : num_objs
                if ~objs_label(i)
                    label = label + 1;
                    objs_label(i) = label;
                    for j = (i + 1) : num_objs
                        if ~objs_label(j)
                            if objs(i).has_the_same_ablation_parameter_Q(objs(j))
                                objs_label(j) = label;
                            end
                        end
                    end
                end
            end
            objs_idx = fun_bin_data_to_idx_list(objs_label);
            objs_bin = fun_bin_data_to_cells_by_ind(objs, objs_idx);
        end
        
        function [z_r_um, m_list] = merge_ablation_parameter_on_each_plane(objs, options)
            arguments
                objs (:, 1) WBIMSPAblation
                options.z_r_um_range (1, :) double = [];
            end
            if isempty(options.z_r_um_range)
                z_list = arrayfun(@(x) reshape(repmat(x.abl_z_r_um.', 1, numel(x.material)), [], 1), ...
                    objs, 'UniformOutput', false);
            else
                z_list = arrayfun(@(x) reshape(repmat(x.get_abl_z_r_um('abl_range_r_um', options.z_r_um_range).',...
                    1, numel(x.material)), [], 1), objs, 'UniformOutput', false);
            end
            num_z = cellfun(@numel, z_list);
            m_list = arrayfun(@(x) reshape(repmat(x.material, num_z, 1), [], 1), ...
                objs, 'UniformOutput', false);
%             m_list = arrayfun(@(x) reshape(repmat(x.material, numel(x.abl_z_r_um), 1), [], 1), ...
%                 objs, 'UniformOutput', false);
            
            z_list = cat(1, z_list{:});
            m_list = cat(1, m_list{:});
            [z_r_um, ~, tmp_ic] = unique(z_list);
            tmp_z_idx = fun_bin_data_to_idx_list(tmp_ic);
            m_list = fun_bin_data_to_cells_by_ind(m_list, tmp_z_idx);
        end
        
        function selected_parameters = select_by_ablation_material(objs, abl_material)
            arguments
                objs WBIMSPAblation
                abl_material WBIMAblationMaterial
            end
            num_objs = numel(objs);
            selected_Q = false(num_objs, 1);
            for i = 1 : num_objs
                if any(objs(i).material == abl_material, 'all')
                    selected_Q(i) = true;
                end
            end
            selected_parameters = objs(selected_Q);
        end
        
        function is_unique_Q = is_safe_for_referse_refinement(objs)
            arguments
                objs WBIMSPAblation
            end
            material_list = cat(2, objs.material);
            is_unique_Q = (numel(material_list) == numel(unique(material_list)));
        end
    end
    
    %% Computing parameters
    methods(Static)
        function spot_shape = estimate_ablation_area(diffuserQ)
            % Estimation based on hardware parameter. Independent of the
            % machine state except for the diffuser
            spot_shape = struct;
            tf_mag = WBIMConfig.ABLATION_LENS_F_mm / WBIMConfig.OBJECTIVE_F_mm;
            diff_incident_angle_deg = asind(WBIMConfig.ABLATION_WAVELENGTH_nm / 1e3 ...
                / (WBIMConfig.ABLATION_GRATING_SPACING_mm * 1e3));
            d_e2_grating_mm = WBIMConfig.ABLATION_BEAM_DIAMETER_e2_mm;
            if isscalar(d_e2_grating_mm)
                % [nondiffracted direction, diffraction direction]
                d_e2_grating_mm = [d_e2_grating_mm, ...
                    d_e2_grating_mm / cosd(diff_incident_angle_deg)];
            end
            % Further divided by sqrt(2) to account for the effective size
            % of uniform fluence (for gaussian beam)
            im_d_e2_um = d_e2_grating_mm * 1e3 / tf_mag / sqrt(2);
            if WBIMConfig.ABLATION_CYLINDRICAL_LENS_Q
                spot_shape.diff_deff_um = im_d_e2_um(2);
                spot_shape.nondiff_deff_um = 1.22 * WBIMConfig.ABLATION_WAVELENGTH_nm * 1e-3 / ...
                    (max(1, d_e2_grating_mm(1) / WBIMConfig.OBJECTIVE_BACK_APERTURE_SIZE_mm) * ...
                    WBIMConfig.OBJECTIVE_NA);
                if diffuserQ
                    % Add the diameter directly
                    diff_d_addon_um = WBIMConfig.ABLATION_DIFFUSER_FWHM_Deg * pi / 180 * ...
                        WBIMConfig.OBJECTIVE_F_mm * 1e3;
                    spot_shape.diff_deff_um = spot_shape.diff_deff_um + diff_d_addon_um;
                    spot_shape.nondiff_deff_um = spot_shape.nondiff_deff_um + diff_d_addon_um;
                end
            else % Achromatic lens
                spot_shape.diff_deff_um = im_d_e2_um(2);
                spot_shape.nondiff_deff_um = im_d_e2_um(1);
            end
            % Determine the orientation of the ablation pattern w.r.t.
            % stage axes
            if WBIMConfig.ABLATION_DIFFRACTION_STAGE_AXIS == WBIMConfig.STAGE_AXIS_FAST_ID
                spot_shape.deff_um_fast = spot_shape.diff_deff_um;
                spot_shape.deff_um_slow = spot_shape.nondiff_deff_um;
            elseif WBIMConfig.ABLATION_DIFFRACTION_STAGE_AXIS == WBIMConfig.STAGE_AXIS_SLOW_ID
                spot_shape.deff_um_slow = spot_shape.diff_deff_um;
                spot_shape.deff_um_fast = spot_shape.nondiff_deff_um;
            else
                error('ABLATION_DIFFRACTION_STAGE_AXIS should be either 1 or 2');
            end
            
        end
        
        function p_str = compute_ablation_parameters(site_scale_factor_fast, ...
                site_scale_factor_slow, diffuserQ, peak_fluence_J_cm2)
            % Compute ablation parameter for a singel set of parameters
            arguments
                site_scale_factor_fast (1,1) double
                site_scale_factor_slow (1,1) double
                diffuserQ (1,1) logical
                peak_fluence_J_cm2 (1,1) double
            end
            p_str = WBIMSPAblation.estimate_ablation_area(diffuserQ);
            p_str.site_scale_factor_fast = site_scale_factor_fast;
            p_str.site_scale_factor_slow = site_scale_factor_slow;
            p_str.diffuserQ = diffuserQ;
            p_str.peak_fluence_J_cm2 = peak_fluence_J_cm2;
            
            % Speed
            p_str.fast_axis_speed_um_s = p_str.deff_um_fast .* WBIMConfig.ABLATION_REPETITION_RATE_Hz ...
                .* site_scale_factor_fast;
            % Limiting the speed directly does not make sense. The speed is
            % direclty related to the effective fluence 
%             if p_str.fast_axis_speed_um_s > WBIMConfig.FAST_AXES_MAX_SPEED_um_s
%                 p_str.fast_axis_speed_um_s = WBIMConfig.FAST_AXES_MAX_SPEED_um_s;
%             end
            p_str.slow_axis_step_um = p_str.deff_um_slow .* site_scale_factor_slow;
            
            % Acceleration
            stage_max_acc_m_s2 = WBIMConfig.XY_PEAK_THRUST_N ./  WBIMConfig.FAST_AXIS_TOTAL_LOAD_kg;
            res_max_acc_m_s2 = WBIMConfig.ABLATION_RESONANT_FREQUENCY_Hz .* ...
                WBIMConfig.ABLATION_ACC_MAX_FREQUENCY_FRACTION .* (p_str.fast_axis_speed_um_s ./ 1e6);
            max_acc_m_s2 = min([stage_max_acc_m_s2, res_max_acc_m_s2, ...
                WBIMConfig.XY_MAX_ACCELERATION_m_s2]);
            if WBIMConfig.ABLATION_ACC_MAX_M_S2 > max_acc_m_s2
                p_str.fast_acceleration_m_s2 = max_acc_m_s2;
            else
                p_str.fast_acceleration_m_s2 = WBIMConfig.ABLATION_ACC_MAX_M_S2;
            end
            p_str.fast_acceleration_eff_frequency_Hz = p_str.fast_acceleration_m_s2 ./ (p_str.fast_axis_speed_um_s/1e6);
            p_str.fast_acceleration_length_um = (p_str.fast_axis_speed_um_s).^2 ./ ...
                (2 .* p_str.fast_acceleration_m_s2 .* 1e6) ...
                .* WBIMConfig.ABLATION_ACC_LENGTH_EXPANSION_RATIO;
            
            % Ablation laser power after the PBS
            p_str.power_W = WBIMSPAblation.peak_fluence_J_per_cm2_to_output_power_W(...
                peak_fluence_J_cm2, diffuserQ);
        end
        
        function p_str = compute_ablation_parameters_for_multiple_planes(...
                site_scale_factor_fast, site_scale_factor_slow, diffuserQ, peak_fluence_J_cm2)
            
            max_num_elem = max([numel(site_scale_factor_fast), ...
                numel(site_scale_factor_slow), numel(diffuserQ), numel(peak_fluence_J_cm2)]);
            if isscalar(site_scale_factor_fast)
                site_scale_factor_fast = repelem(site_scale_factor_fast, 1, max_num_elem);
            else
                assert(numel(site_scale_factor_fast) == max_num_elem, ...
                    'Mismatch array size');
            end
            if isscalar(site_scale_factor_slow)
                site_scale_factor_slow = repelem(site_scale_factor_slow, 1, max_num_elem);
            else
                assert(numel(site_scale_factor_slow) == max_num_elem, ...
                    'Mismatch array size');
            end
            if isscalar(diffuserQ)
                diffuserQ = repelem(diffuserQ, 1, max_num_elem);
            else
                assert(numel(diffuserQ) == max_num_elem, ...
                    'Mismatch array size');
            end
            if isscalar(peak_fluence_J_cm2)
                peak_fluence_J_cm2 = repelem(peak_fluence_J_cm2, 1, max_num_elem);
            else
                assert(numel(peak_fluence_J_cm2) == max_num_elem, ...
                    'Mismatch array size');
            end
            plane_p = arrayfun(@WBIMSPAblation.compute_ablation_parameters, ...
                site_scale_factor_fast, site_scale_factor_slow, diffuserQ, peak_fluence_J_cm2);
            % Convert structure
            field_names = fieldnames(plane_p(1));
            for i = 1 : numel(field_names)
                tmp_fn = field_names{i};
                p_str.(tmp_fn) = [plane_p.(tmp_fn)];
            end
        end
        
        function power_W = peak_fluence_J_per_cm2_to_output_power_W(fluence, diffuserQ)
            % Estimate the power right after the ablation half wave plate
            num_f = numel(fluence);
            num_ds = numel(diffuserQ);
            max_numel = max([num_f, num_ds]);
            if max_numel > 1
                if num_f == 1
                    fluence = repelem(fluence, 1, max_numel);
                else
                    assert(num_f == max_numel);
                end
                if num_ds == 1
                    diffuserQ = repelem(diffuserQ, 1, max_numel);
                else
                    assert(num_ds == max_numel);
                end
            end
            power_W = zeros(1, max_numel);
            for i = 1 : max_numel
                abl_shape = WBIMSPAblation.estimate_ablation_area(diffuserQ(i));
                power_W(i) = fluence(i) * (abl_shape.diff_deff_um * abl_shape.nondiff_deff_um* ...
                    pi / 4 / 1e8) * WBIMConfig.ABLATION_REPETITION_RATE_Hz;
                if diffuserQ(i)
                    power_W(i) = power_W(i) / WBIMConfig.ABLATION_TRANSMISSION_FRACTION_WITH_DIFFUSER;
                else
                    % No diffuser
                    power_W(i) = power_W(i) / WBIMConfig.ABLATION_TRANSMISSION_FRACTION;
                end
            end
        end
        
    end
end