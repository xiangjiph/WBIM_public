classdef WBIMAVD < handle
    % WBIM Ablation volume detection
    
    properties
        smip_cell (1, :) cell
        local_tile_info % seems that only scan_roi_smip_bbox_pxl is used for removing the ablation traget outside the scan region 
        acquisition_mode (1,1) WBIMMicroscopeMode
        has_skull_Q (1,1) logical
        
        nonempty_cell_idx (1, :) double
        mask_size (1, 3) double {mustBeInteger}
        roi_dist_to_edge_map single
        foreground_mask_2d logical
        
        detection_section_range (1, :) double
        % Copied properties
        region_bbox_ctr_yx_um (1, 2) double
        voxel_size_um (1, 3) double
        abs_fp_z_um_mm (1, :) double
        step_mip_abs_fp_z_um (1, :) double
        
        mask_cell (1, :) cell         
        material_mask_cell (2, :) cell
        labeled_mask 
        itp % mask interpolation object
        section_cc_stat (1, :) cell
    end
    
    % Properties that go into the ablation path generation
    properties
        ablation_parameter
        ablation_mask (:, 1) cell % sample space yx um 
        ablation_z_r_um (:, 1) cell % relative z um for each plane of the ablation mask 
        ablation_path (:, 1) WBIMAblationPath2Dv2        
    end
    
    properties(Access=private)
        detect_ablation_object_Q_ (1,1) logical = false;
    end
    %% Dependent variables
    properties(Dependent)
        mask_size_um (1, 3) double
        need_ablation_Q (1,1) logical
        detect_ablation_object_Q (1, 1) logical
        detection_range_um (1, 2) double
    end
    
    methods
        function val = get.mask_size_um(obj)
            val = obj.mask_size .* obj.voxel_size_um;
        end
        
        function val = get.need_ablation_Q(obj)
            val = ~isempty(obj.ablation_path);
        end
        
        function val = get.detect_ablation_object_Q(obj)
            if ~obj.detect_ablation_object_Q_
                if ~isempty(obj.labeled_mask) && any(obj.labeled_mask(:, :, ...
                        obj.detection_section_range(1):obj.detection_section_range(2)), 'all')
                    obj.detect_ablation_object_Q_ = true;
                end
            end
            val = obj.detect_ablation_object_Q_;
        end
        
        function val = get.detection_range_um(obj)
            if ~isempty(obj.detection_section_range) && ~isempty(obj.voxel_size_um)
                start_z_um = (obj.detection_section_range(1) - 1) * obj.voxel_size_um(3);
                end_z_um = obj.detection_section_range(2) * obj.voxel_size_um(3);
                val = [start_z_um, end_z_um];
            else
                val = nan(1, 2);
            end
        end
        
    end
    
    %%
    methods
        function obj = WBIMAVD(smip_cell, local_tile_info, options)
            arguments
                smip_cell (1, :) cell
                local_tile_info
                options.acquisition_mode (1,1) WBIMMicroscopeMode
                options.sec_range (1,:) double = [];
                options.has_skull_Q (1,1) logical
            end
            % Initialization
            obj.smip_cell = smip_cell;
            obj.local_tile_info = local_tile_info;
            obj.acquisition_mode = options.acquisition_mode;
            obj.has_skull_Q = options.has_skull_Q;
            % Copy fields 
            obj.region_bbox_ctr_yx_um = local_tile_info.region_bbox_ctr_yx_um;
            obj.voxel_size_um = local_tile_info.step_mip_pixel_yxz_um;
            obj.abs_fp_z_um_mm = local_tile_info.abs_fp_z_um_done;
            obj.step_mip_abs_fp_z_um = local_tile_info.step_mip_abs_fp_z_um;
            
            if ~isempty(options.sec_range)
                obj.detection_section_range = options.sec_range;
            end
            % Parse input mask list
            is_not_empty_Q = ~cellfun(@isempty, smip_cell);
            if any(is_not_empty_Q)
                obj.nonempty_cell_idx = find(is_not_empty_Q);
                obj.mask_size = size(smip_cell{obj.nonempty_cell_idx(1)});
                if isempty(obj.detection_section_range)
                    obj.detection_section_range = [1, obj.mask_size(3)];
                end
            end            
            % Construct roi mask and distance map
            obj.construct_roi_dt_map();
        end
        
        function mask = compute_single_channel_masks(obj, options)
            arguments
                obj (1,1) WBIMAVD         
                % There should be no break-through from tdTomato to SHG.
                % Add -0.05 for more conservative segmentation of the SHG
                % mask 
                options.spum_coeff = [    1, -0.03,    0,     0; ...
                                      -0.05,     1,    0,     0; ...
                                          0,     0,    1,     0; ...
                                          0, -0.05,    0,     1];
                options.merge_ch_14_Q (1, 1) logical = true;
                options.sec_range (1, :) double = obj.detection_section_range;
                options.visQ (1,1) logical = false;
                
                options.channel_unmixing_Q (1, 1) logical = true;
            end
            % Apply spectral unmixing to the step MIP images
            obj.smip_cell = WBIMAVD.smip_spectral_unmixing(obj.smip_cell, 'coeff', ...
                options.spum_coeff);
            % Apply channel unmixing before merging the two channels. 
            if options.merge_ch_14_Q && numel(obj.smip_cell) == 4 && ~isempty(obj.smip_cell{1}) && ...
                    ~isempty(obj.smip_cell{4})
                % Treat both vessel and tdTomato-labeled object tissue. Not
                % sure if max will amplifiy the noise too severely or not. 
                obj.smip_cell{1} = max(obj.smip_cell{1}, obj.smip_cell{4});
%                 obj.smip_cell{1} = (obj.smip_cell{1} + obj.smip_cell{4}) ./ 2;
                obj.smip_cell{4} = [];                
            end
            
            switch obj.acquisition_mode
                case WBIMMicroscopeMode.Scan
                    obj.mask_cell = obj.compute_scan_ablation_mask(obj.smip_cell);                    
                case WBIMMicroscopeMode.Explore
                    obj.mask_cell = obj.compute_exploration_ablation_mask(obj.smip_cell, ...
                        'sec_range', options.sec_range);
                    % Exclude regions outside the scan roi
                    assert(~isempty(obj.roi_dist_to_edge_map), 'roi_dist_to_edge_map is empty');
                    obj.mask_cell = obj.exploration_mask_exclude_unscanned_region(...
                        obj.mask_cell, obj.roi_dist_to_edge_map > 0);
                otherwise 
                    error('To be implemented');
            end
            mask = obj.mask_cell;
        end
        
        function labeled_mask = construct_labeled_mask(obj, mask_cell, options)
            arguments
                obj (1,1) WBIMAVD
                mask_cell (1,:) cell = obj.mask_cell
                options.has_skull_Q (1,1) logical = obj.has_skull_Q
                options.detect_vsl_in_skull_Q (1,1) logical = true;
            end
            labeled_cell = obj.construct_material_mask(mask_cell, 'has_skull_Q', ...
                options.has_skull_Q, 'detect_vsl_in_skull_Q', options.detect_vsl_in_skull_Q);
            num_ch = size(labeled_cell, 2);
            labeled_mask = zeros(obj.mask_size, 'uint8');
            for i = 1 : num_ch
                label = labeled_cell{1, i};
                mask = labeled_cell{2, i};
                if ~isempty(mask)
                    labeled_mask = max(labeled_mask, uint8(mask) .* label);
                end
            end
            obj.labeled_mask = labeled_mask;   
        end
        
        function obj = construct_merged_ablation_mask(obj, abl_para_sets, options)
            arguments
                obj (1,1) WBIMAVD
                abl_para_sets (1, :) WBIMSPAblation
                options.visQ (1,1) logical = false
            end
            % Remove the ablation materials absent in the labeled mask 
            exist_material = unique(obj.labeled_mask);         
            abl_para_sets = WBIMSPAblation.select_by_ablation_material(...
                abl_para_sets, exist_material);            
            abl_para_binned = WBIMSPAblation.bin_by_ablation_parameter(abl_para_sets);
            num_abl_sets = numel(abl_para_binned);
            obj.ablation_parameter = cellfun(@(x) x(1), abl_para_binned);
            [obj.ablation_mask, obj.ablation_z_r_um] = deal(cell(num_abl_sets, 1));
            is_valid_Q = false(num_abl_sets, 1);
            for i = 1 : num_abl_sets
                % A little bit wasting time here. Exclude ablation
                % material that are absent in the mask before this step? 
                
                % WARNING: Providing the option of 'z_r_um_range' will
                % overwrite the z_r_um. 
                [z_r_um, m_list] = WBIMSPAblation.merge_ablation_parameter_on_each_plane(...
                    abl_para_binned{i});
%                 [z_r_um, m_list] = WBIMSPAblation.merge_ablation_parameter_on_each_plane(...
%                     abl_para_binned{i}, 'z_r_um_range', obj.detection_range_um);
                mask_sample_yx_um = obj.compute_interpolated_ablation_mask(...
                    z_r_um, m_list);                                
                if any(mask_sample_yx_um, 'all')
                    is_valid_Q(i) = true;
                    % Determine if any layer is emtpy 
%                     is_valid_section_Q = squeeze(any(mask_sample_yx_um, [1,2]));
%                     obj.ablation_mask{i} = mask_sample_yx_um(:, :, is_valid_section_Q);
%                     obj.ablation_z_r_um{i} = z_r_um(is_valid_section_Q);
                    obj.ablation_mask{i} = mask_sample_yx_um;
                    obj.ablation_z_r_um{i} = z_r_um;
                end
            end            
            obj.ablation_mask = obj.ablation_mask(is_valid_Q);
            obj.ablation_z_r_um = obj.ablation_z_r_um(is_valid_Q);
            obj.ablation_parameter = obj.ablation_parameter(is_valid_Q);
        end
        
        function material_cell = construct_material_mask(obj, mask_cell, options)
            arguments
                obj (1,1) WBIMAVD
                mask_cell (1, :) cell = obj.mask_cell
                
                options.detect_vsl_in_skull_Q (1,1) logical = true
                options.max_v2s_dist_um (1, 1) double = 200;
                
                options.has_skull_Q (1,1) logical = obj.has_skull_Q
            end            
            if options.detect_vsl_in_skull_Q && ismember(WBIMChannelName.Vessel, ...
                    obj.nonempty_cell_idx) && ismember(WBIMChannelName.SHG, obj.nonempty_cell_idx)
                max_v2s_dist_pxl = options.max_v2s_dist_um / prod(obj.voxel_size_um)^(1/3);                
                
                [mask_cell{WBIMChannelName.Vessel}, mask_vis, mask_cell{WBIMChannelName.SHG}] ...
                    = WBIMAVD.compute_vessel_in_skull_mask(obj.smip_cell{[WBIMChannelName.Vessel, WBIMChannelName.SHG]}, ...
                    mask_cell{[WBIMChannelName.Vessel, WBIMChannelName.SHG]}, 'max_v2s_dist_pxl', max_v2s_dist_pxl);
            else
                mask_vis = [];
            end
            % Might need to merge different channels that are all tissue
            % later. 
            num_ch = numel(obj.nonempty_cell_idx);
            material_cell = cell(2, 0);            
            for i = 1 : num_ch
                tmp_ch = obj.nonempty_cell_idx(i);
                if ~isempty(mask_cell{tmp_ch})
                    switch tmp_ch
                        case WBIMChannelName.Vessel
                            label = WBIMAblationMaterial.Tissue;
                        case WBIMChannelName.SHG
                            if options.has_skull_Q
                                label = WBIMAblationMaterial.Bone;
                            else
                                label = WBIMAblationMaterial.Membrane;
                            end
                        otherwise
                            label = WBIMAblationMaterial.Tissue;
                            %                         fprintf('To be implemented: construct_material_mask for channel %d\n', ...
                            %                             tmp_ch);
                    end
                    material_cell(:, end+1) = {label, mask_cell{i}};
                end
            end
            
            if ~isempty(mask_vis)
               % Add vessel in skull mask to the label array 
               material_cell(:, end+1) = {WBIMAblationMaterial.TissueInSkull, mask_vis};               
            end            
            % Sort 
            material_list = cat(1, material_cell{1, :});
            [~, sort_idx] = sort(material_list, 'ascend');
            material_cell = material_cell(:, sort_idx);
            
            obj.material_mask_cell = material_cell;
        end
        
        function cc_stat = compute_section_cc_statistics(obj, opts)
            arguments
                obj (1, 1) WBIMAVD
                opts.dist_bin_step (1, 1) double = 5;
                opts.section_range (1, :) double = obj.detection_section_range;
            end
            if isfield(obj.local_tile_info, 'scan_roi_smip_bbox_pxl')
                bbox_pxl = obj.local_tile_info.scan_roi_smip_bbox_pxl;
            else
                bbox_pxl = [1,1, obj.mask_size(1:2)];
            end        
            % Construct scan roi distance transform field 
            
            scan_dt_map = obj.roi_dist_to_edge_map;
            max_scan_dt = max(scan_dt_map(:));
            dist_bin_edge = 0 : opts.dist_bin_step : (ceil(max_scan_dt / opts.dist_bin_step) ...
                * opts.dist_bin_step);            
            
            fg_mask = WBIMAVD.estimate_foreground_mask_2d(obj.smip_cell, 'roi_mask', ...
                scan_dt_map > 0, 'int_th', 1e3, 'imclose_r', 2, 'fill_holes', true);
            fg_dt_map = WBIMAVD.fg_bwdist(fg_mask);
            
            % TODO: Compare the fg map and scan_dt_map to determine if the
            % foreground itself is abnormal e.g. has a big hole inside 
            
            bbox_ll = bbox_pxl(3:4) - bbox_pxl(1:2) + 1;
            bbox_area = prod(bbox_ll);
            pixel_size_um = prod(obj.voxel_size_um(1:2))^(1/2);
            
            num_material = size(obj.material_mask_cell, 2);
            cc_stat = cell(num_material, 1);
            for i = 1 : num_material
                tmp_stat = obj.compute_single_frame_cc_stat(obj.material_mask_cell{2, i}, ...
                    fg_dt_map, 'detection_section_range', opts.section_range, ...
                    'dist_bin_edge', dist_bin_edge);
                tmp_stat.total_area_fraction = tmp_stat.total_area / bbox_area;
                tmp_stat.max_area_fraction = tmp_stat.max_area / bbox_area;
                tmp_stat.max_eff_r_um = tmp_stat.max_eff_r * pixel_size_um;
                
                tmp_stat.avg_dist_um = tmp_stat.avg_dist * pixel_size_um;
                
                tmp_stat.nonedge_eff_r_um = tmp_stat.nonedge_cc_eff_r * pixel_size_um;
                tmp_stat.nonedge_cc_avg_dist_um = tmp_stat.nonedge_cc_avg_dist * pixel_size_um;
                
                tmp_stat.section = opts.section_range(1) : opts.section_range(2);
                % Use -1 (instead of - 0.5) if we don't want to scale the
                % ablation power on the last ablation plane
                tmp_stat.z_rel_um = (tmp_stat.section - 0.5) * obj.voxel_size_um(3);
                tmp_stat.material = obj.material_mask_cell{1, i};                
                cc_stat{i} = tmp_stat;
            end            
            obj.section_cc_stat = cc_stat;
            obj.foreground_mask_2d = fg_dt_map;
        end
                
        function [dist_to_edge_map] = construct_roi_dt_map(obj)
            if isfield(obj.local_tile_info, 'scan_roi_smip_bbox_pxl')
                bbox_pxl = obj.local_tile_info.scan_roi_smip_bbox_pxl;
                roi_mask = false(obj.mask_size(1:2));
                roi_mask(bbox_pxl(1) : bbox_pxl(3), bbox_pxl(2) : bbox_pxl(4)) = true;
            else
                roi_mask = true(obj.mask_size(1:2));
            end
            dist_to_edge_map = WBIMAVD.fg_bwdist(roi_mask);
            obj.roi_dist_to_edge_map = dist_to_edge_map;
        end   
        
        function nonedge_cc_Q = exist_surface_nonedge_cc_Q(obj)
            nonedge_cc_info = obj.analyze_cc_statistics(obj.section_cc_stat);
            nonedge_cc_Q = nonedge_cc_info.exist_Q;
        end
    end
    %% Common
    methods(Static)
        function mask_ht = compute_mask_by_thresholding(smip, options)
            % Assume smip_vsl has the same size as smip_shg
            arguments
                smip
                options.threshold (1,:) double = [2e3, 1e4];
                options.min_cc_vxl (1,1) double = 8; % 8 for skull
                % For tissue
                options.dilation_radius_pxl (1,1) {mustBeInteger} = 0;
                options.cmaxQ (1,1) logical = false;
                options.min_2d_hole_pxl (1,1) {mustBeInteger} = 0;
                % For bones
                options.bwclose_sph_r_pxl (1,1) {mustBeInteger} = 0; % 9 for skull
            end
            % Intensity thresholding 
            if isscalar(options.threshold)
                mask_ht = smip > options.threshold;
            elseif numel(options.threshold) == 2
                mask_ht = hysterisis_thresholding(smip, options.threshold(1), ...
                    options.threshold(2));
            else
                error('threshold shoud be either a scalar or a 2-element vector');
            end
            % Remove small ccs in 3D 
            if options.min_cc_vxl
                mask_ht = bwareaopen(mask_ht, options.min_cc_vxl);
            end
            % Assume the voxel size of both images to be isotropic -
            if options.dilation_radius_pxl
                mask_ht = imdilate(mask_ht, strel('sphere', options.dilation_radius_pxl));
            end
            if options.bwclose_sph_r_pxl
                mask_ht = imclose(mask_ht, strel('sphere', options.bwclose_sph_r_pxl));
            end
            % Remove additional 2D holes on sections at the end
            if options.min_2d_hole_pxl
                mask_ht = WBIMAcqPostProcessing.filter_sections(mask_ht, ...
                    @bw_fill_small_holes, options.min_2d_hole_pxl);
            end
            % For vessels, assume there's always tissue underneath the
            % tissue - not the case for skull 
            if options.cmaxQ
                mask_ht = cummax(mask_ht, 3);
            end
        end
        
        function [abl_mask, varargout] = ablation_mask_top_surface_refinement_frac_z_peak(...
                smip_im, smip_mask, options)
            arguments
                smip_im
                smip_mask
                options.sec_range (1,:) double = [1, size(smip_im, 3)];
                options.max_peak_int (1,1) double = inf;
                options.peak_int_fraction (1,1) double = 0.5;
                % Local maximum detection
                options.local_max_wd_sz (1,1) double = 3
                options.local_max_min_frac_max_int (1,1) double = 0.2;
                options.local_max_min_int (1,1) double = 1e3;
                % Smooth surface map
                options.sf_map_imc_disk_r_pxl (1,1) double = 2;
                % Clean up ablation mask
                options.abl_mask_min_cc_pxl_2d (1,1) double = 0;
                options.visQ (1,1) logical = false;
            end
            if isempty(options.sec_range)
                options.sec_range = [1, size(smip_im, 3)];
            end            
            smip_im(isnan(smip_im)) = 0;
            % Reshape array to get the z intensity profiles
            mask_yx_size = size(smip_im, [1,2]);
            num_sections = size(smip_im, 3);
            sec_idx_list = options.sec_range(1) : options.sec_range(2);
            search_map = any(smip_mask(:, :, sec_idx_list), 3);
            smip_im_mip = max(smip_im(:, :, sec_idx_list), [], 3);
            % Estimate background level. Use median or mean?
            est_bg_med = median(smip_im(~smip_mask), 'omitnan');
            mask_pxl_ind = find(search_map & (smip_im_mip > (est_bg_med + eps(est_bg_med))));
            num_pxl = numel(mask_pxl_ind);
            % Get masked pixel z intensity profile
            pxl_int_profile = reshape(permute(smip_im, [3,1,2]), num_sections, []);
            pxl_int_profile = pxl_int_profile(:, mask_pxl_ind);
            % Search for the local maxima in the intensity profile
            int_lm_Q = movmax(pxl_int_profile, options.local_max_wd_sz, 1);
            % Local maximum should be higher than the minimum masked pixel
            % value, and higher than 20% of the peak intenisty value
            % (along the axial direction)
            z_int_max = max(pxl_int_profile, [], 1);
            int_lm_Q = pxl_int_profile >= int_lm_Q & ...
                pxl_int_profile > max(options.local_max_min_int(1), (z_int_max * options.local_max_min_frac_max_int));
            [~, int_lm_idx] = max(int_lm_Q, [], 1);
            local_max_int = pxl_int_profile(sub2ind(size(pxl_int_profile), int_lm_idx, 1:num_pxl));
            if isfinite(options.max_peak_int)
                local_max_int = min(local_max_int, options.max_peak_int);
            end

            tmp_exp_th = (local_max_int - est_bg_med) * options.peak_int_fraction + est_bg_med;
            % Fractional max intensity point:
            tmp_exp_check_int_ab_th_Q = (pxl_int_profile > tmp_exp_th);            
            hm_Q_diff = diff(cat(1, false(1, num_pxl), tmp_exp_check_int_ab_th_Q, ...
                false(1, num_pxl)));
            % Return the first greater than fractional maximum value position
            [max_Q_value, first_gt_fm_idx] = max(hm_Q_diff, [], 1);
            assert(all(max_Q_value==1));
            % Construct the output
            result = struct;
            result.mask = smip_mask;
            result.pxl_ind = mask_pxl_ind;
            result.peak_int = z_int_max;
            result.local_max_int = local_max_int;
            lm_int_sec_map = zeros(mask_yx_size);
            lm_int_sec_map(mask_pxl_ind) = int_lm_idx;
            result.lm_sec_map = lm_int_sec_map;
            % z-index of the surface in the mask
            surf_sec_map = zeros(mask_yx_size);
            surf_sec_map(mask_pxl_ind) = first_gt_fm_idx;
            % Smooth the map
            surf_sec_map = medfilt2(surf_sec_map);
            if options.sf_map_imc_disk_r_pxl > 0
                surf_sec_map = imclose(surf_sec_map, strel('disk', options.sf_map_imc_disk_r_pxl));
            end
            result.surface_sec_map = surf_sec_map;
            tmp_surf_map = surf_sec_map;
            tmp_surf_map(tmp_surf_map==0) = inf;
            surf_mask = repmat(permute(1:num_sections, [3,1,2]), mask_yx_size);
            % @ge for including the surface; @gt for excluding the surface
            surf_mask = bsxfun(@ge, surf_mask, tmp_surf_map);
            
            result.surface_mask = surf_mask;
            result.abl_mask = result.surface_mask & result.mask;
            if options.abl_mask_min_cc_pxl_2d
                % Remove small component per layer
                result.abl_mask = bw_area_open_in_sections(result.abl_mask, options.abl_mask_min_cc_pxl_2d);
            end
            result.abl_mask_in_range_mip = any(result.abl_mask(:, :, sec_idx_list), 3);
            abl_mask = result.abl_mask;
            if nargout > 1
                varargout{1} = result;
            end
            
            if options.visQ
                fig_hdl = figure;
                num_sp = 4;
                ax_0 = subplot(2,2,1);
                imagesc(ax_0, result.surface_sec_map);
                ax_0.Title.String = 'Surface layer';
                cmap = colormap(ax_0);
                cmap(1, :) = [0,0,0];
                colormap(ax_0, cmap);
                colorbar(ax_0);
                ax_1 = subplot(2,2,2);
                imagesc(ax_1, result.lm_sec_map);
                colorbar(ax_1);
                ax_1.Title.String = 'Peak intensity layer';
                colormap(ax_1, cmap);
                ax_2 = subplot(2,2,3);
                imagesc(ax_2, result.abl_mask_in_range_mip);
                colorbar(ax_2);
                ax_2.Title.String = 'Ablation mask in range';
                colormap(ax_2, 'gray');
                ax_3 = subplot(2,2,4);
                th_map = zeros(size(result.surface_sec_map));
                th_map(mask_pxl_ind) = tmp_exp_th;
                imagesc(ax_3, th_map);
                ax_3.Title.String = 'Threshold map';
                colorbar(ax_3);
                [ax_0.DataAspectRatio, ax_1.DataAspectRatio, ...
                    ax_2.DataAspectRatio, ax_3.DataAspectRatio] = deal([1,1,1]);
                if nargout > 2
                    varargout{2} = fig_hdl;
                end
            end
        end
        
        function [shg_mask_ts, varargout] = ablation_mask_bottom_surface_refinement_skull(...
                s_smip_shg, shg_mask, options)
            % For detecting the bottom surface of the thin skull shell
            % Based on the effective minimum scattering length (intensity
            % drops fastest along the z direciton)
            arguments
                s_smip_shg
                shg_mask logical
                options.min_dI (1,1) double = 250;
                options.max_thin_shell_ls_um (1,1) double = 10;
                options.axial_pixel_size_um (1,1) double = 5;
                options.max_refine_int (1,1) double = 8e3;
                % Only consider normalized gradient with original intensity greater than this value
                options.min_search_int (1,1) double = 1e3;
                options.max_gradient_offset (1,1) double = 1;
                options.mask_sm_r_pxl (1,1) double = 2;
                options.visQ (1,1) logical = false;
            end
            % Derived parameter
            min_neg_dI_n = options.axial_pixel_size_um / options.max_thin_shell_ls_um;
            
            s_smip_size = size(s_smip_shg);
            num_pxl = prod(s_smip_size(1:2));
            num_z = s_smip_size(3);
            shg_mask_v = reshape(permute(shg_mask, [3,1,2]), num_z, []);
            % Compute gradient and normalized gradient along the z
            % direction
            [dI_n, dI] = WBIMAcqPostProcessing.compute_z_normalized_derivative(s_smip_shg);
            % Remove positive gradient and noise
            is_invalid_Q = dI_n > 0 | s_smip_shg < options.min_search_int | abs(dI) < options.min_dI;
            dI_n(is_invalid_Q) = nan;
            [max_neg_dI_n, max_idx] = max(-dI_n, [], 3);
            % Refine masks:
            %   Remove mask pixels below the maximum normalized gradient
            %   point, down to the next mask starting point (if the 1D mask
            %   has more than 1 ablation intervals)
            shg_mask_2d = any(shg_mask, 3);
            s_smip_shg_v_le_max_int_Q = reshape(permute(s_smip_shg, [3,1,2]), num_z, []) ...
                <= options.max_refine_int;
            shg_mask_v_edge = diff(cat(1, false(1, num_pxl), shg_mask_v, false(1, num_pxl)));
            thin_shell_mask_2d = shg_mask_2d & (max_neg_dI_n > min_neg_dI_n);
            thin_shell_ind = find(thin_shell_mask_2d);
            thin_shell_mask_3d = false(size(shg_mask_v));
            for i = 1 : numel(thin_shell_ind)
                tmp_ind = thin_shell_ind(i);
                tmp_max_g_depth_idx = max_idx(tmp_ind);
                for j = (tmp_max_g_depth_idx + options.max_gradient_offset) : num_z
                    if s_smip_shg_v_le_max_int_Q(j, tmp_ind) 
                        thin_shell_mask_3d(j, tmp_ind) = true;
                    end
                    if shg_mask_v_edge(j+1, tmp_ind) == 1
                        break
                    end
                end
            end
            thin_shell_mask_3d = permute(reshape(thin_shell_mask_3d, s_smip_size([3,1,2])), [2,3,1]);
            if options.mask_sm_r_pxl
                thin_shell_mask_3d = imclose(thin_shell_mask_3d, ...
                    strel('disk', options.mask_sm_r_pxl));
            end
            shg_mask_ts = shg_mask & (~thin_shell_mask_3d);
            
            if options.visQ
                max_ls_um = options.axial_pixel_size_um ./max_neg_dI_n;
                max_idx(isnan(max_ls_um)) = nan;
                max_grad_depth_um = max_idx * options.axial_pixel_size_um;
                
                fig_hdl = figure;
                ax_1 = nexttile();
                imagesc(ax_1, max_ls_um);
                ax_1.DataAspectRatio = [1,1,1];
                c_hdl_1 = colorbar(ax_1);
                ax_1.Title.String = 'Minimum scattering legth (\mum)';
                ax_1.CLim = [0, 50];
                ax_2 = nexttile();
                imagesc(ax_2, max_grad_depth_um);
                ax_2.DataAspectRatio = [1,1,1];
                c_hdl_2 = colorbar(ax_2);
                %                 ax_2.CLim = [0, 150];
                ax_2.Title.String = 'Maximum gradient depth (\mum)';
                %                 [ax_1.XLabel.String, ax_2.XLabel.String] = deal(sprintf('X (%d \\mum)',  ...
                %                     s_region_info.step_mip_pixel_yxz_um(2)));
                %                 [ax_1.YLabel.String, ax_2.YLabel.String] = deal(sprintf('Y (%d \\mum)',  ...
                %                     s_region_info.step_mip_pixel_yxz_um(1)));
                ax_3 = nexttile();
                imagesc(ax_3, any(thin_shell_mask_3d, 3));
                ax_3.DataAspectRatio = [1,1,1];
                ax_3.Title.String = sprintf('l_s < %.1f \\mum mask', options.max_thin_shell_ls_um);
                ax_4 = nexttile();
                h_bin = 0 : 1 : 50;
                h_hdl = histogram(ax_4, max_ls_um, 'BinEdges', h_bin);
                ax_4.XLabel.String = 'Minimum scattering length (\mum)';
                ax_4.YLabel.String = 'Counts';
            else
                fig_hdl = [];
            end
            if nargout > 1
                varargout{1} = fig_hdl;
            end
        end
        
        function mask_bottom_map = compute_mask_surface_map(mask_3d, dir, opts)
            arguments
                mask_3d
                dir (1,1) {mustBeMember(dir, [1, -1])} = 1
                opts.visQ (1,1) logical = false;
            end
            % Compute the z position of the lowest point of a 3D mask. If
            % the entire column of voxels are false, return nan
            mask_size = size(mask_3d);
            % Find candidate pixels (2d)
            mask_shg_mip = any(mask_3d, 3);
            mask_shg_ind_2d = find(mask_shg_mip);
            % Get mask profile along z
            mask_z_prof = permute(mask_3d, [3, 1, 2]);
            mask_z_prof = reshape(mask_z_prof, mask_size(3), []);
            mask_z_prof = mask_z_prof(:, mask_shg_ind_2d);
            num_pxl = numel(mask_shg_ind_2d);
            mask_z_prof_d = diff(cat(1, false(1, num_pxl), mask_z_prof, false(1, num_pxl)), 1);
            % Find the bottom
            switch dir
                case 1
                    [tmp_val, surface_idx] = max(mask_z_prof_d, [], 1);
                    assert(all(tmp_val == 1));
                case -1
                    [tmp_val, surface_idx] = min(mask_z_prof_d(end:-1:1, :), [], 1);
                    assert(all(tmp_val == -1));
                    surface_idx = size(mask_z_prof_d, 1) - surface_idx;
            end            
            mask_shg_sub_2d = fun_ind2sub(mask_size(1:2), mask_shg_ind_2d);
            % Interpolate
            skull_bottom_itp = scatteredInterpolant(mask_shg_sub_2d(:, 1), mask_shg_sub_2d(:, 2), surface_idx.', ...
                'linear', 'none');
            [itp_sub_1, itp_sub_2] = ndgrid(1 : mask_size(1), 1 : mask_size(2));
            mask_bottom_map = skull_bottom_itp(itp_sub_1, itp_sub_2);
            
            if opts.visQ
                fig_hdl = figure;
                ax_hdl_1 = nexttile;
                imagesc(ax_hdl_1, mask_shg_mip);
                ax_hdl_1.DataAspectRatio = [1,1,1];
                
                ax_hdl_2 = nexttile;
                imagesc(ax_hdl_2, mask_bottom_map);
                ax_hdl_2.DataAspectRatio = [1,1,1];
                switch dir
                    case 1
                        ax_hdl_2.Title.String = 'Top surface depth';
                    case -1
                        ax_hdl_2.Title.String = 'Bottom surface depth';
                end
            end
        end
        
        function smip_im_d = smip_spectral_unmixing(smip_im, options)
            arguments
                smip_im
                % options.coeff(j, i) is the fraction of the j-th channel
                % intensity appear in the i-th channel 
                options.coeff double = WBIMConfig.IMAGING_CHANNEL_UNMIXING_MATRIX;
            end
            num_ch = numel(smip_im);
            smip_im_d = cell(size(smip_im));
            if ~all(cellfun(@isfloat, smip_im))
                im_type = cellfun(@class, smip_im, 'UniformOutput', false);
                smip_im = cellfun(@single, smip_im, 'UniformOutput', false);
            else
                im_type = {};
            end
                
            for i = 1 : num_ch
                if ~isempty(smip_im{i})
                    tmp_im = smip_im{i};
                    tmp_rf = 1;
                    for j = 1 : num_ch
                        if ~isempty(smip_im{j}) && (i ~= j) && options.coeff(j, i)
                            tmp_im = tmp_im + options.coeff(j, i) .* smip_im{j};
                            tmp_rf = tmp_rf + options.coeff(j, i);
                        end
                    end
                    % Rescale 
                    tmp_im = max(0, tmp_im) ./ tmp_rf;
                    if ~isempty(im_type)
                        tmp_im = cast(tmp_im, im_type{i});
                    end
                    smip_im_d{i} = tmp_im;
                end
            end            
        end
        
        function cc_stat = compute_single_frame_cc_stat(mask_3d, scan_dt_map, options)
            arguments
               mask_3d
               scan_dt_map (:, :) single
               options.nonedge_cc_min_dist_pxl = 5;
               options.detection_section_range = [];
               options.dist_bin_edge (1, :) double = [];
            end
            if isempty(options.detection_section_range)
                num_frame = size(mask_3d, 3);
            else
                assert(numel(options.detection_section_range) == 2);
                mask_3d = mask_3d(:, :, options.detection_section_range(1) : options.detection_section_range(2));
                num_frame = options.detection_section_range(2) - options.detection_section_range(1) + 1;
            end
            
            if isempty(options.dist_bin_edge)
                dt_max = max(scan_dt_map(:));
                bin_edge = 5;
                options.dist_bin_edge = 0: bin_edge : (ceil(dt_max / bin_edge) * bin_edge);
            end
            
            cc_stat = struct;
            cc_stat.dist_bin_val = (options.dist_bin_edge(1:end-1) + ...
                options.dist_bin_edge(2:end)) / 2;
            num_bins = numel(cc_stat.dist_bin_val);
            
            [cc_stat.num_cc, cc_stat.total_area, cc_stat.max_eff_r, cc_stat.max_area, ...
                cc_stat.avg_dist, ...
                cc_stat.nonedge_cc_area, cc_stat.nonedge_cc_avg_dist, ...
                cc_stat.nonedge_cc_eff_r] = ...
                deal(nan(1, num_frame));
            % Reverse accumulated number of pixels at distance from the
            % boundary of the scan_dt_map
            cc_stat.rcumm_dt = nan(num_bins, num_frame);            
            for i = 1 : num_frame
                tmp_im = mask_3d(:, :, i);
                if any(tmp_im, 'all')
                    tmp_cc = bwconncomp(tmp_im);                    
                    cc_size = cellfun(@numel, tmp_cc.PixelIdxList);
                    cc_stat.num_cc(i) = tmp_cc.NumObjects;
                    cc_stat.total_area(i) = sum(cc_size);
                    cc_stat.max_area(i) = max(cc_size);
                    
                    cc_dist2edge = cellfun(@(x) scan_dt_map(x), tmp_cc.PixelIdxList, ...
                        'UniformOutput', false);
                    min_dist_to_edge = cellfun(@min, cc_dist2edge);
                    avg_dist_to_edge = cellfun(@mean, cc_dist2edge);
                    
                    is_nonedge_cc_Q = min_dist_to_edge >= options.nonedge_cc_min_dist_pxl;
                    if any(is_nonedge_cc_Q)
                        nonedge_cc_size = cc_size(is_nonedge_cc_Q);
                        nonedge_cc_avg_dist = avg_dist_to_edge(is_nonedge_cc_Q);
                        cc_stat.nonedge_cc_area(i) = sum(nonedge_cc_size);
                        cc_stat.nonedge_cc_avg_dist(i) = sum(nonedge_cc_avg_dist .* nonedge_cc_size)...
                            ./ cc_stat.nonedge_cc_area(i);
                    end                    
                    pxl_dist_to_edge = scan_dt_map(tmp_im);
                    dist_count = histcounts(pxl_dist_to_edge, options.dist_bin_edge);
                    cc_stat.avg_dist(i) = mean(pxl_dist_to_edge, 'omitnan');
                    cc_stat.rcumm_dt(:, i) = cumsum(dist_count, 'reverse');                    
                end
            end
            cc_stat.max_eff_r = sqrt(cc_stat.max_area / pi);
            cc_stat.nonedge_cc_eff_r = sqrt(cc_stat.nonedge_cc_area / pi);
        end
        
        function fg_mask = estimate_foreground_mask_2d(smip_cell, opts)
            arguments
                smip_cell (1, :) cell
                opts.int_th (1, :) double = 2e3;
                opts.imclose_r (1, 1) double = 0;
                opts.fill_holes (1, 1) logical = true;
                opts.roi_mask = [];
            end
            num_cell = numel(smip_cell);
            if isscalar(opts.int_th)
                opts.int_th = repelem(opts.int_th, 1, num_cell);
            end
            fg_mask = [];
            for i = 1 : num_cell
                tmp_smip = smip_cell{i};
                if ~isempty(tmp_smip)
                    if isempty(fg_mask)
                        fg_mask = false(size(tmp_smip, [1,2]));
                    end
                    fg_mask = fg_mask | any(tmp_smip >= opts.int_th(i), 3);
                end
            end
            if opts.imclose_r
                fg_mask = imclose(fg_mask, strel('disk', opts.imclose_r));
            end
            if opts.fill_holes
                fg_mask = imfill(fg_mask, 'holes');
            end
            if ~isempty(opts.roi_mask)
                fg_mask = fg_mask & opts.roi_mask;
            end
        end
        
        function dist_map = fg_bwdist(roi_mask)
            % Compute the distance to the boundary of the background
            bg_mask = ~roi_mask;
            bg_mask = padarray(bg_mask, [1, 1], 1, 'both');
            dist_map = bwdist(bg_mask);
            dist_map = dist_map(2:end-1, 2:end-1);
        end
        
        % Quality control 
        function nonedge_cc_info = analyze_cc_statistics(cc_stat, opt)
            arguments
                cc_stat (:, 1) cell % output of WBIMADV.compute_section_cc_statistics
                opt.top_range (1, 1) double = 3;
                opt.bottom_range (1, 1) double = 3;
                opt.exclude_last_sec_Q (1, 1) logical = true;
                
                % Mimimal cc size for triggering warning:
                opt.t_s_r_um (1, 1) double = 20; % top, skull, min effective radius 
                opt.t_s_d_um (1, 1) double = 50; % top, skull, min dist to edge
                opt.b_s_r_um (1, 1) double = 40; % bottom, skull, min effective radius
                opt.b_s_d_um (1, 1) double = 50; % bottom, skull, min dist to edge
                
                opt.t_t_r_um = 50; % top, tissue, min effective radius
                opt.t_t_d_um = 100; % top, tissue, min dist to edge
                % opt.b_t_r_um_h = 40;
                % opt.b_t_d_um_h = 50;
            end
            sec_list = cc_stat{1}.section;
            num_sec = numel(sec_list);
            top_idx = 1 : min(opt.top_range, num_sec);
            bottom_idx = num_sec - (opt.bottom_range : -1 : 1) + 1;
            if opt.exclude_last_sec_Q
                bottom_idx = bottom_idx(1 : end - 1);
            end
            nonedge_cc_info = struct('top_skull_Q', false, 'bottom_skull_Q', false, ...
                'top_tissue_Q', false, 'bottom_tissue_Q', false, 'exist_Q', false);
            for i_stat = 1 : numel(cc_stat)
                tmp_stat = cc_stat{i_stat};
                if ~isempty(tmp_stat)
                    % Statistics near the top of the exploration refinement volume
                    top_nonedge_r_um = mean(tmp_stat.nonedge_eff_r_um(top_idx), 'omitnan');
                    top_nonedge_dist_um = mean(tmp_stat.nonedge_cc_avg_dist_um(top_idx), 'omitnan');
                    % Statistics near the bottom of the exploration refinement volume
                    bottom_nonedge_r_um = mean(tmp_stat.nonedge_eff_r_um(bottom_idx), 'omitnan');
                    bottom_nonedge_dist_um = mean(tmp_stat.nonedge_cc_avg_dist_um(bottom_idx), 'omitnan');
                    
                    if tmp_stat.material == WBIMAblationMaterial.Bone
                        if top_nonedge_r_um >= opt.t_s_r_um && ...
                                top_nonedge_dist_um > opt.t_s_d_um
                            nonedge_cc_info.top_skull_Q = true;
                        end
                        if bottom_nonedge_r_um >= opt.b_s_r_um && ...
                                bottom_nonedge_dist_um >= opt.b_s_d_um
                            nonedge_cc_info.bottom_skull_Q = true;
                        end
                    else
                        % Not sure what to do with the tissue at the moment.
                        % Wait for example to occur...
                        if top_nonedge_r_um >= opt.t_t_r_um && ...
                                top_nonedge_dist_um > opt.t_t_d_um
                            nonedge_cc_info.top_tissue_Q = true;
                        end
                    end
                end
            end
            nonedge_cc_info.exist_Q = nonedge_cc_info.top_skull_Q || ...
                nonedge_cc_info.bottom_skull_Q || nonedge_cc_info.top_tissue_Q ...
                || nonedge_cc_info.bottom_tissue_Q;
        end
    end
    
    methods
        function [output] = construct_mask_yx_um_itp(obj, options)
            arguments
                obj (1,1) WBIMAVD
                options.interpolation = 'nearest';
                options.extrapolation = 'nearest';
            end
            sub_1 = (0.5 : obj.mask_size(1) - 0.5) * obj.voxel_size_um(1);
            sub_2 = (0.5 : obj.mask_size(2) - 0.5) * obj.voxel_size_um(2);
            sub_3 = (0.5 : obj.mask_size(3) - 0.5) * obj.voxel_size_um(3);
            [sub_1, sub_2, sub_3] = ndgrid(sub_1, sub_2, sub_3);
            obj.itp = griddedInterpolant(sub_1, sub_2, sub_3, ...
                single(obj.labeled_mask), options.interpolation, options.extrapolation);
            output = obj.itp;
        end
        
        function mask_sample_yx_um = compute_interpolated_ablation_mask(obj, ...
                abl_z_r_um, material_list)
            arguments
                obj (1,1) WBIMAVD
                abl_z_r_um (1, :) double
                material_list (:, 1) cell
            end
            [sub1, sub2, sub3] = ndgrid(1 : obj.mask_size_um(1), ...
                1 : obj.mask_size_um(2), abl_z_r_um);
            num_z = numel(abl_z_r_um);
            assert(num_z == numel(material_list));
            tmp_array = obj.itp(sub1, sub2, sub3);
            mask_sample_yx_um = false(size(tmp_array));
            for i = 1 : num_z
                tmp_m = uint8(material_list{i});
                mask_sample_yx_um(:, :, i) = any(bsxfun(@eq, tmp_array(:, :, i), ...
                    reshape(tmp_m, 1,1,[])), 3);
            end
        end
    end
    %% Scan mode
    methods(Static)
        % This function might be merged with compute_exploration_ablation_mask
        % The differences in options can be merged.
        function abl_vol_cell = compute_scan_ablation_mask(smip_cell, ...
                options)
            arguments
                smip_cell
                options.est_vsl_th = [1e3, 3e3]
                options.est_shg_th = [3e3, 6e3] % Changed from [1.5e3, 3e3] 04/10/2024
                options.smip_mask_min_cc_size_pxl (1,:) double = 27; % in 3D mask
                
                options.vsl_mask_cmax_Q (1,1) logical = true;
                % Segmentation (operations in order)
                options.vsl_bwd_r_pxl (1,1) double = 6;
                options.vsl_bwc_sph_r_pxl (1,1) double = 3;
                options.vsl_min_2d_hole_pxl (1,1) double = 500;
                    % Do not dilate the skull mask 
                options.shg_bwd_r_pxl (1,1) double = 0;
                options.shg_bwc_sph_r_pxl (1,1) double = 0;
                options.shg_min_2d_hole_pxl (1,1) double = 0;
                
                % Surface detection 
                options.refine_surface_Q (1,1) logical = true; % false disables the following parameters
                options.sec_range (1, :) double = [1, size(smip_cell{1}, 3)];
                options.surf2peak_int_ratio (1,1) double = 0.5;
                options.tissue_peak_int = 9e3;
                % Bottom detection 
                options.refine_bottom_Q (1,1)logical = true; % false disables the following parameters
                
                options.vsl_min_cc_pxl_2d (1,1) double = 9;
                options.shg_min_cc_pxl_2d (1,1) double = 4;
                options.visQ (1,1) logical = false;
            end
            num_ch = numel(smip_cell);
            [abl_vol_cell, info_cell] = deal(cell(1, num_ch));
            for i_ch = 1 : num_ch
                if ~isempty(smip_cell{i_ch})
                    switch i_ch
                        % Do not refine top for both - to ensure a clean
                        % surface cut 
                        case WBIMChannelName.Vessel
                            % Compute mip mask array
                            [abl_vol_cell{i_ch}, info_cell{i_ch}] = WBIMAVD.compute_ablation_mask_vessel(...
                                smip_cell{i_ch}, 'threshold', options.est_vsl_th, ...
                                'bwclose_sph_r_pxl', options.vsl_bwc_sph_r_pxl, ...
                                'dilation_radius_pxl', options.vsl_bwd_r_pxl, ...
                                'min_2d_hole_pxl', options.vsl_min_2d_hole_pxl, ...
                                'abl_mask_min_cc_pxl_2d', options.vsl_min_cc_pxl_2d,...
                                'refine_top_Q', options.refine_surface_Q, ...
                                'cmaxQ', options.vsl_mask_cmax_Q, 'visQ', options.visQ);
                        case WBIMChannelName.SHG
                            [abl_vol_cell{i_ch}, info_cell{i_ch}] = WBIMAVD.compute_ablation_mask_SHG(...
                                smip_cell{i_ch}, 'threshold', options.est_shg_th, ...
                                'bwclose_sph_r_pxl', options.shg_bwc_sph_r_pxl,...
                                'dilation_radius_pxl', options.shg_bwd_r_pxl, ...
                                'min_2d_hole_pxl', options.shg_min_2d_hole_pxl, ...
                                'abl_mask_min_cc_pxl_2d', options.shg_min_cc_pxl_2d,...
                                'refine_top_Q', options.refine_surface_Q, ...
                                'refine_bottom_Q', options.refine_bottom_Q);
                        case WBIMChannelName.tdTomato
                            % To be implemented
                        otherwise
                            continue;
                    end
                end
            end
        end
    end
    %% Exploration mode
    methods(Static)
        function abl_vol_cell = compute_exploration_ablation_mask(smip_cell, ...
                options)
            arguments
                smip_cell
                options.est_vsl_th = [1e3, 3e3]
                options.est_shg_th = [7e2, 2e3]
                options.smip_mask_min_cc_size_pxl (1,:) double = 27; % in 3D mask
                options.vsl_mask_cmax_Q (1,1) logical = false;
                % Surface detection 
                options.sec_range (1, :) double = [1, size(smip_cell{1}, 3)];
                options.surf2peak_int_ratio (1,1) double = 0.5;
                options.shg_peak_int = inf;
                options.vsl_peak_int = 9e3;
                
                options.vsl_min_cc_pxl_2d (1,1) double = 25;
                % skull debris in the tissue 
                options.shg_min_cc_pxl_2d (1,1) double = 0; 
                options.visQ (1,1) logical = false;
            end
            if isempty(options.sec_range)
                options.sec_range = [1, size(smip_cell{1}, 3)];
            end
            num_ch = numel(smip_cell);
            [abl_vol_cell, info_cell] = deal(cell(1, num_ch));
            for i_ch = 1 : num_ch
                if ~isempty(smip_cell{i_ch})
                    switch i_ch
                        case WBIMChannelName.Vessel
                            [abl_vol_cell{i_ch}, info_cell{i_ch}] = WBIMAVD.compute_ablation_mask_vessel(...
                                smip_cell{i_ch}, 'threshold', options.est_vsl_th, ...
                                'bwclose_sph_r_pxl', 0, 'dilation_radius_pxl', 0, 'min_2d_hole_pxl', 0, ...
                                'refine_top_Q', true, 'sec_range', options.sec_range, ... 
                                'max_peak_int', options.vsl_peak_int, ...
                                'abl_mask_min_cc_pxl_2d', options.vsl_min_cc_pxl_2d,...
                                'cmaxQ', options.vsl_mask_cmax_Q, 'visQ', options.visQ);
                            % TODO: Add cc analysis to ignore small pieces above
                            % the tissue surface                            
                        case WBIMChannelName.SHG
                            [abl_vol_cell{i_ch}, info_cell{i_ch}] = WBIMAVD.compute_ablation_mask_SHG(...
                                smip_cell{i_ch}, 'threshold', options.est_shg_th, ...
                                'bwclose_sph_r_pxl', 0, 'dilation_radius_pxl', 0, 'min_2d_hole_pxl', 0, ...
                                'refine_top_Q', true, 'sec_range', options.sec_range, ...
                                'surf2peak_int_ratio', options.surf2peak_int_ratio, ...
                                'abl_mask_min_cc_pxl_2d', options.shg_min_cc_pxl_2d,...
                                'refine_bottom_Q', true, 'max_peak_int', options.shg_peak_int, ...
                                'visQ', options.visQ);
                        case WBIMChannelName.tdTomato
                            % To be implemented                            
                        otherwise
                            continue
                    end
                end
            end
            % Combined analysis 
            % Smooth? remove small ccs? 
            
            % Zero regions outside the scan ROI
            
        end
        
        function inQ = check_tiles_in_mask_2D(mask_yx, local_bbox_mmxx_um, mask_pixel_size_um)
            arguments
                mask_yx (:, :) logical
                local_bbox_mmxx_um (:, 4) double
                mask_pixel_size_um (1, 2) double = [1,1] % By default, the mask is at 1 um resolution 
            end
            mask_size = size(mask_yx);
            local_bbox_mmxx_pxl = local_bbox_mmxx_um ./ [mask_pixel_size_um, mask_pixel_size_um];
            local_bbox_mmxx_um = round(min(max(1, local_bbox_mmxx_pxl), [mask_size, mask_size]));
            num_bbox = size(local_bbox_mmxx_um, 1);
            inQ = false(num_bbox, 1);
            for i = 1 : num_bbox
                tmp_bbox_mmxx = local_bbox_mmxx_um(i, :);
                tmp_mask = mask_yx(tmp_bbox_mmxx(1) : tmp_bbox_mmxx(3), ...
                    tmp_bbox_mmxx(2) : tmp_bbox_mmxx(4));
                inQ(i) = any(tmp_mask, 'all');
            end
        end
        
        function abl_vol_cell = exploration_mask_exclude_unscanned_region(...
                abl_vol_cell, scan_roi_mask)
            arguments
               abl_vol_cell (1, :) cell
               scan_roi_mask logical
            end
            for i_c = 1 : numel(abl_vol_cell)
                if ~isempty(abl_vol_cell{i_c})
                    abl_vol_cell{i_c} = logical(abl_vol_cell{i_c} .* scan_roi_mask);
                end
            end
        end
        
        function [exp_abl_p] = adjust_ablation_parameters_based_on_section_cc_stat(exp_abl_p, ...
                cc_stat, options)
            arguments
                exp_abl_p (1, :) WBIMSPAblation
                cc_stat (1, :) cell
                
                options.diff_amp_exp_lz_um (1,1) double = 600;
                options.line_amp_exp_lz_um (1,1) double = 200;
                % Extra amplification at the top
                options.diff_amp_top_extra_ratio (1,1) = 1.25;
                options.line_amp_top_extra_ratio (1,1) = 1.25;
                % Constant shift
                options.line_amp_baseline (1, 1) double {mustBePositive} = 1.25;
                options.abl_z_r_um_list (1, :) double = [];
                
                options.visQ (1, 1) logical = false
            end
            num_stat = numel(cc_stat);
            num_abl_p = numel(exp_abl_p);
            checked_Q = false(num_abl_p, 1);

            for i = 1 : num_stat
                tmp_stat = cc_stat{i};                
                assert(isscalar(tmp_stat.material));
                for j = 1 : numel(exp_abl_p)
                    tmp_p = exp_abl_p(j);
                    tmp_z_range = tmp_p.abl_range_r_um;
                    
                    if ismember(tmp_stat.material, tmp_p.material) && ~checked_Q(j)
                        tmp_dist_to_top_um = tmp_stat.z_rel_um;
                        tmp_z_in_range_Q = tmp_dist_to_top_um >= tmp_z_range(1) & ...
                            tmp_dist_to_top_um <= tmp_z_range(2);
                        if any(tmp_z_in_range_Q)
                            tmp_dist_to_top_um = tmp_dist_to_top_um(tmp_z_in_range_Q);
                            tmp_dist_to_bottom_um = tmp_dist_to_top_um(end) - tmp_dist_to_top_um;
                            
                            if tmp_p.with_diffuser_Q
                                tmp_stat.amp_factor = exp(tmp_dist_to_bottom_um / options.diff_amp_exp_lz_um);
                                tmp_stat.amp_factor = [tmp_stat.amp_factor(1) * options.diff_amp_top_extra_ratio, tmp_stat.amp_factor];
                            else
                                tmp_stat.amp_factor = exp(tmp_dist_to_bottom_um / options.line_amp_exp_lz_um);
                                tmp_stat.amp_factor = [tmp_stat.amp_factor(1) * options.line_amp_top_extra_ratio, tmp_stat.amp_factor];
                                tmp_stat.amp_factor = tmp_stat.amp_factor .* options.line_amp_baseline;
                            end
                            tmp_stat.amp_itp = griddedInterpolant([0, tmp_dist_to_top_um], tmp_stat.amp_factor, 'linear', 'nearest');
                            
                            % Moderest amplification
                            if isempty(options.abl_z_r_um_list)
                                tmp_amp_factor = tmp_stat.amp_itp(tmp_p.abl_z_r_um);
                            else
                                tmp_amp_factor = tmp_stat.amp_itp(options.abl_z_r_um_list);
                            end
                            tmp_p.peak_fluence_J_cm2 = tmp_p.peak_fluence_J_cm2 .* tmp_amp_factor;
                        end
                        checked_Q(j) = true;
                        exp_abl_p(j) = tmp_p;
                    end
                end
            end
        end
        
    end
    
    methods
        function grid_ind = get_exploration_refined_tile_ind(obj)
            arguments
                obj (1,1) WBIMAVD
            end
            overall_mask = any(cat(3, obj.ablation_mask{:}), 3);
            inQ = WBIMAVD.check_tiles_in_mask_2D(overall_mask, ...
                obj.local_tile_info.local_tile_bbox_mmxx_um, [1,1]);
            grid_ind = obj.local_tile_info.grid_ind(inQ);            
        end
    end
    
    %% Vessel
    methods(Static)
        function [vsl_mask, varargout] = compute_ablation_mask_vessel(smip_vsl, options)
            % No need for bottom detection 
            % Assume there's always tissue under tissue -> cmaxQ = true
            arguments
                smip_vsl
                % Segmentation 
                options.threshold (1, :) {mustBeNumeric} = [1e3, 3e3];
                options.smip_mask_min_cc_size_pxl (1,:) double = 27; % in 3D mask
                options.bwclose_sph_r_pxl (1,1) {mustBeInteger, mustBeNonnegative} = 0; % 0 for exploration
                options.dilation_radius_pxl (1,1) {mustBeInteger, mustBeNonnegative} = 0; % 0 for exploration
                options.min_2d_hole_pxl (1,1) {mustBeInteger, mustBeNonnegative} = 0; % 0 for exploration
                options.cmaxQ (1,1) logical = false; % Disable for exploration? 
                % Surface detection for ablation refinement 
                options.refine_top_Q (1,1) logical = false;
                options.sec_range (1, 2) double = [1, size(smip_vsl, 3)];
                options.surf2peak_int_ratio (1,1) double = 0.5;
                options.max_peak_int (1,1) double = 7e3;
                options.abl_mask_min_cc_pxl_2d (1,1) double = 0
                % Visualization
                options.visQ (1,1) logical = false;
            end
            % Do not add voxel to the intensity mask on the lateral
            % direction before top surface detection, as some computations
            % relys on the minimum intensity threshold
            % cmaxQ could be true as it is along z 
            vsl_mask = WBIMAVD.compute_mask_by_thresholding(smip_vsl, ...
                'threshold', options.threshold, ...
                'min_cc_vxl', options.smip_mask_min_cc_size_pxl, ...
                'bwclose_sph_r_pxl', options.bwclose_sph_r_pxl, 'cmaxQ', options.cmaxQ, ...
                'dilation_radius_pxl', options.dilation_radius_pxl,...
                'min_2d_hole_pxl', options.min_2d_hole_pxl);
            
            if options.refine_top_Q
                [vsl_mask, surf_info] = WBIMAVD.ablation_mask_top_surface_refinement_frac_z_peak(...
                    smip_vsl, vsl_mask, 'sec_range', options.sec_range, ...
                    'local_max_min_int', options.threshold(1), ...
                    'peak_int_fraction', options.surf2peak_int_ratio,  ...
                    'max_peak_int', options.max_peak_int, ...
                    'abl_mask_min_cc_pxl_2d', options.abl_mask_min_cc_pxl_2d, ...
                    'visQ', options.visQ);
            else
                surf_info = [];
            end 
            if nargout > 1
                varargout{1} = surf_info;
            end
        end
    end
    %% SHG
    % Scan mode does not require top surface detection?
    % Exploration mode does not require bottom surface detection?
    methods(Static)
        function [shg_mask, varargout] = compute_ablation_mask_SHG(smip_shg, options)
            arguments
                smip_shg
                % Segmentation
                options.threshold (1, :) {mustBeNumeric} = [7e2, 3e3];
                options.smip_mask_min_cc_size_pxl (1,:) double = 27; % in 3D mask
                options.bwclose_sph_r_pxl (1,1) {mustBeInteger, mustBeNonnegative} = 0; % 0 for exploration
                options.dilation_radius_pxl (1,1) {mustBeInteger, mustBeNonnegative} = 0; % 0 for exploration
                options.min_2d_hole_pxl (1,1) {mustBeInteger, mustBeNonnegative} = 0; % 0 for exploration
                
                % Surface detection for ablation refinement
                options.refine_top_Q (1,1) logical = false;
                options.sec_range (1, 2) double = [1, size(smip_shg, 3)];
                options.surf2peak_int_ratio (1,1) double = 0.4;
                options.max_peak_int (1,1) double = inf;
                options.abl_mask_min_cc_pxl_2d (1,1) double = 0
                % Thin shell detection
                options.refine_bottom_Q (1,1) logical = false;
%                 options.min_dI (1,1) double = 250;
%                 options.max_thin_shell_ls_um (1,1) double = 8;
%                 options.axial_pixel_size_um (1,1) double = 5;
%                 options.max_refine_int (1,1) double = 8e3;
%                 options.min_search_int (1,1) double = 1e3;
%                 options.max_gradient_offset (1,1) double = 1;
%                 options.mask_sm_r_pxl (1,1) double = 2;       
                % Visualization 
                options.visQ (1,1) logical = false;
            end
            % Do not add voxel to the intensity mask before top surface
            % detection, as some computations relys on the minimum
            % intensity threshold 
            shg_mask = WBIMAVD.compute_mask_by_thresholding(smip_shg, ...
                'threshold', options.threshold, ...
                'min_cc_vxl', options.smip_mask_min_cc_size_pxl, ...
                'bwclose_sph_r_pxl', options.bwclose_sph_r_pxl, 'cmaxQ', false, ...
                'dilation_radius_pxl', options.dilation_radius_pxl,...
                'min_2d_hole_pxl', options.min_2d_hole_pxl);
            if options.refine_top_Q
                [shg_mask, top_surf_info] = WBIMAVD.ablation_mask_top_surface_refinement_frac_z_peak(...
                    smip_shg, shg_mask, 'sec_range', options.sec_range, ...
                    'local_max_min_int', options.threshold(1), ...
                    'peak_int_fraction', options.surf2peak_int_ratio,  ...
                    'max_peak_int', options.max_peak_int, ...
                    'visQ', options.visQ);
            else
                top_surf_info = [];
            end
            if nargout > 1
                varargout{1} = top_surf_info;
            end            
            if options.refine_bottom_Q
                shg_mask = WBIMAVD.ablation_mask_bottom_surface_refinement_skull(smip_shg, ...
                    shg_mask, 'min_search_int', options.threshold(1), ...
                    'visQ', false);
            end
            
            if options.visQ
                vis_sec = 19;
                vis_im = smip_shg(:, :, vis_sec);
                vis_mask = shg_mask(:, :, vis_sec);
                fig_hdl = figure;
                fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* [4, 2];
                a1 = nexttile;
                imshow(vis_im);
                % a1.DataAspectRatio = [1,1,1];
                a2 = nexttile;
                imshow(vis_mask);
                a3 = nexttile;
                imshowpair(vis_im, vis_mask);
            end            
        end        
    end
    %% Channel unmixing 
    methods(Static)
        function [mask_vsl, mask_vslIs, mask_shg] = compute_vessel_in_skull_mask(...
                im_vsl, im_shg, mask_vsl, mask_shg, options)
            arguments
                im_vsl
                im_shg
                mask_vsl
                mask_shg
                % Unmixing parameters
                options.r_vsl2shg (1,1) double = 0.5; % maximum breakthrough intensity ratio from vsl to shg
                options.r_shg2vsl (1,1) double = 0.5; % maximum breakthrough intensity ratio from shg to vsl
                options.min_int_shg (1,1) double = 5000;
                options.min_int_vsl (1,1) double = 5000;    
                % Masking parameters
                options.max_v2s_dist_pxl (1, 1) double = 25 % by default, 125 um, 5 um/pixel     
                options.mask_sm_r_pxl (1, 1) double = 3
                
                options.rm_shg_from_vsl_Q (1, 1) logical = false;
                options.rm_vsl_from_shg_Q (1, 1) logical = true;
                
                options.vis_sec (1, 1) double = -1
                options.save_fp = [];
            end
            if any(mask_shg, 'all')
                [shg_in_vsl_mask, vsl_in_shg_mask] = WBIMAcqPostProcessing.compute_step_mip_spectrum_unmixing_mask(...
                    im_vsl, im_shg, 'r_vsl2shg', options.r_vsl2shg, 'r_shg2vsl', options.r_shg2vsl, ...
                    'min_int_shg', options.min_int_shg, 'min_int_vsl', options.min_int_vsl);
                mask_shg_rf = mask_shg & ~vsl_in_shg_mask;
                mask_shg_rf_dt = bwdist(mask_shg_rf);
                shg_proximity_mask = (mask_shg_rf_dt < options.max_v2s_dist_pxl) & ~mask_shg_rf;
                if options.mask_sm_r_pxl
                    %                 shg_proximity_mask = WBIMAcqPostProcessing.filter_sections(shg_proximity_mask, ...
                    %                     @(x) imclose(x, strel('disk', options.mask_sm_r_pxl)));
                    shg_proximity_mask = imclose(shg_proximity_mask, strel('sphere', options.mask_sm_r_pxl));
                end
                mask_vslIs = (shg_proximity_mask) & mask_vsl;
                mask_vsl = mask_vsl & ~mask_vslIs;
                
                if options.rm_vsl_from_shg_Q
                    mask_shg = mask_shg_rf;
                end
                if options.rm_shg_from_vsl_Q
                    mask_vsl = mask_vsl & ~shg_in_vsl_mask;
                end
            else
                [vsl_in_shg_mask, mask_vslIs] = deal(false(size(mask_shg)));
            end
            
            if options.vis_sec > 0
                fig_hdl = figure;
                fig_hdl.Position(3) = fig_hdl.Position(3) * 3;
                ax1 = nexttile;
                imshowpair(im_vsl(:, :, options.vis_sec), im_shg(:, :, options.vis_sec));
                ax1.Title.String = "Merged image";
                ax2 = nexttile;
                imshowpair(im_vsl(:, :, options.vis_sec), mask_shg(:, :, options.vis_sec));
                ax2.Title.String = "Vessel image with SHG mask";
%                 ax3 = nexttile;
%                 imshowpair(im_vsl(:, :, options.vis_sec), mask_vsl(:, :, options.vis_sec));
%                 ax3.Title.String = "Vessel image with vessel mask";
                ax4 = nexttile;
                imshowpair(im_vsl(:, :, options.vis_sec), vsl_in_shg_mask(:, :, options.vis_sec));
                ax4.Title.String = "Vessel image with SHG leakage mask";
                ax5 = nexttile;
                imshowpair(im_vsl(:, :, options.vis_sec), mask_vslIs(:, :, options.vis_sec));
                ax5.Title.String = "Vessel image with fine ablation mask";
                
                if ~isempty(options.save_fp)
                   fun_print_image_in_several_formats(fig_hdl, options.save_fp); 
                end
            end            
        end
        
    end
    %% Visualization
    methods
        function fig_hdl = vis_single_section_image_and_mask(obj, ch_id, sec, options)
            % TODO: add coordiante 
            arguments
                obj WBIMAVD
                ch_id (1,1) double
                sec (1,1) double
                options.displayQ (1,1) logical = true
                options.save_dir = []
            end
            if isempty(obj.mask_cell{ch_id})
                return
            end
            im = obj.smip_cell{ch_id}(:, :, sec);
            mask = obj.mask_cell{ch_id}(:, :, sec);
            mask = repmat(mask, 1, 1, 3);
%             mask(:, :, 3) = 0;
%             mask(:, :, 2) = 0;
            if options.displayQ
                fig_hdl = figure('Visible', 'on');
            else
                fig_hdl = figure('Visible', 'off');
            end
            t = tiledlayout(1,1, 'Padding', 'compact');
            ax_hdl = axes(fig_hdl);
            im_hdl_1 = imagesc(ax_hdl, im);
            colormap(ax_hdl, 'jet');
            hold(ax_hdl, 'on');
            im_hdl_2 = image(ax_hdl, mask);
            im_hdl_2.AlphaData = 0.4;
            ax_hdl.DataAspectRatio = [1,1,1];
            ax_hdl.Title.String = sprintf('%s Section %d', WBIMChannelName(ch_id), sec);
            ax_hdl.XAxis.Visible = 'off';
            ax_hdl.YAxis.Visible = 'off';
            
            if ~isempty(options.save_dir)
                filename = sprintf('CH%d_Sec_%d.jpg', ch_id, sec);
                if ~isfolder(options.save_dir)
                    mkdir(options.save_dir);
                end
                filepath = fullfile(options.save_dir, filename);
                exportgraphics(fig_hdl, filepath, 'Resolution', 300);
            end
            if ~options.displayQ
                delete(fig_hdl);
            end
        end
        
        function vis_single_channel_image_and_mask(obj, ch_id, options)
            arguments
                obj WBIMAVD
                ch_id (1,1) double
                options.section_list (1, :) double = 1 : obj.mask_size(3)
                options.displayQ (1,1) logical = true
                options.save_dir = []
            end
            for i = options.section_list
                obj.vis_single_section_image_and_mask(ch_id, i, 'displayQ', options.displayQ, ...
                    'save_dir', options.save_dir);
            end            
        end
    end
end