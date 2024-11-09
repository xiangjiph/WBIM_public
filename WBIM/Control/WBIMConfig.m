classdef WBIMConfig < handle
    % Configuration parameters
    properties(Constant=true, Hidden)
        %% Objective & Optics
        OBJECTIVE_F_mm = 7.2;
        OBJECTIVE_RADIUS_um = 18e3; % 18 mm
        OBJECTIVE_POS_WRT_HOLE_mm = [1.9, 1.8] * 25.4; % Upper left corner of the hole (sample image coordinate)
        OBJECTIVE_POS_WRT_STAGE_um = [];
        OBJECTIVE_DIST_TO_BREADBOARD_BUTTON_mm = 28.6;
        OBJECTIVE_WORKING_DISTANCE_mm = 8;
        OBJECTIVE_BACK_APERTURE_SIZE_mm = 13.68;
        OBJECTIVE_NA = 0.95;
        
        POST_SCANNER_MAG = 3.25;
        SCANNER_RANGE_DEG = [26, 20]; % bi-directional.         
        %% Discovery Laser
        IMAGING_LASER_PORT_ID = "COM7";
        IMAGING_LASER_BAUD_RATE = 19200;        
        ILSM_PERIOD_S = 2.5;
        ILSM_BUFFER_SIZE = 60;
        
        IMAGING_LASER_GDD = 15000;
        
        IMAGING_LASER_PULSE_SPACING_s = 1 / 80e6; % 1/repetition rate, s
        
        IMAGING_LASER_MIN_NORMAL_POWER_mW = 2000;
        IMAGING_LASER_NORMAL_POWER_STD_mW = 50;
        IMAGING_LASER_WAIT_DURATION_s = 10;
        IMAGING_LASER_WAIT_COUNT = 30;
        
        %% Astrella Laser
        ABLATION_WAVELENGTH_nm = 800;
        ABLATION_LASER_PULSE_ENERGY_mJ = 5;
        ABLATION_REPETITION_RATE_Hz = 1000;   
        ABLATION_BEAM_DIAMETER_e2_mm = 11 * 1.25; % After 1.25x beam expander
        % The z-position difference between the imaging focal plane and the
        % ablation focal plane
        FOCAL_PLANE_OFFSET_Z_um = [];
        ASTRELLA_GATE_PORT_NAME = '/vDAQ0/D1.0'
        ASTRELLA_GATE_PORT_CHANNEL_NAME = 'D1.0'
        %% Piezo
        PIEZO_RANGE_um = [-200, +200];
        PIEZO_REFERENCE_POS_um = -50;
        %% 3-axes stage parameters
        STAGE_SERIAL_PORT_ID = "COM4";
        STAGE_LIMIT_um = [75000, 75000, 75000];
        XY_STAGE_LIMIT_um = [75000, 75000];
        % Axis 1: 38.2 mm; Axis 2: 37.35 mm. Add extra distance for safty 
        XY_STAGE_HOME_MIN_POSITION_um = [38200, 37350] + 3000; 
        FAST_AXES_MAX_SPEED_um_s = 50000; % 5.0 cm/s - manually set
        FAST_AXIS_MAX_ACCELERATION_m_s2 = 0.1; 
        SLOW_AXIS_MAX_ACCELERATION_m_s2 = 0.1; 
        SLOW_AXES_MAX_SPEED_um_s = 5000; % 0.5 cm/s - manually set
        XY_PEAK_THRUST_N = 13;
        XY_MAX_CONTINUOUS_THRUST_N = 6;
        XY_STAGE_LOAD_kg = 0.43;
        XY_MAX_ACCELERATION_m_s2 = 1.0;        
        
        STAGE_EST_SAMPLE_LOAD_kg = 0.7;
        FAST_AXIS_TOTAL_LOAD_kg = 1.13; % STAGE_EST_SAMPLE_LOAD_kg + XY_STAGE_LOAD_kg
        
        % Affine transformation matrix (homogeneous coordinate) 
        % Configuration 1: 
            % Slow stage (X) positive direction toward Astrella, paralle to
            % the shorter edge of the table
            % Fast stage (Y) positive direction toward the lasers
%         A_STAGE_TO_SAMPLE = [-1, 0, 0, 75000;...
%             0, 1, 0, 0;...
%             0, 0, 1, 0;...
%             0, 0, 0, 1];
        % Configuration 2:
            % Slow stage (X) positive direction toward Astrella, paralle to
            % the shorter edge of the table 
            % Fast stage (Y) positive direction toward the lasers
        A_STAGE_TO_SAMPLE = [0, 1, 0, 0; ...
            1, 0, 0, 0;...
            0, 0, 1, 0;...
            0, 0, 0, 1];
        
        STAGE_AXIS_SLOW_ID = 1;
        STAGE_AXIS_FAST_ID = 2;
        STAGE_AXIS_Z_ID = 3;
        STAGE_AXIS_ACUATOR_ID = 4;
        STAGE_NUM_AXIS = 4;
        
        %
        %% Configuration space of the actuator
        ACTUATOR_DICHROIC_POS_um = 43000;
        ACTUATOR_POS_MAX_um = 44000;
        ACTUATOR_POS_TOLERENCE_um = 50;    
        ACTUATOR_MIRROR_POS_um = 0;        
        %% Configuration space of the XY stages  
        % For sample holder
        STAGE_X_MAX_um = 37500 + 15000;
        STAGE_X_MIN_um = 37500 - 15000;
        STAGE_Y_MAX_um = 37500 + 20700;
        STAGE_Y_MIN_um = 37500 - 20700;
        % For imaging without ablation in disk 
%         STAGE_X_MAX_um = 37500 + 18000;
%         STAGE_X_MIN_um = 37500 - 18000;
%         STAGE_Y_MAX_um = 37500 + 20700;
%         STAGE_Y_MIN_um = 37500 - 20700;        
%         STAGE_Z_MAX_um = 40000;
        STAGE_Z_MAX_um = 62000;
        %% Configuration space of the Z stages
        DIST_BREADBOARD_BUTTON_TO_TABLE_mm = 285.5;
        STAGE_ADAPTER_PLATE_THICKNESS_mm = 12.7;
        STAGE_MINIMUM_HEIGHT_mm = 187.46;
        SAMPLE_ADAPTER_PLATE_THICKNESS_mm = 3/8 * 25.4;
        SAMPLE_HOLDER_THICKNESS_mm = 10;
        %% Microscope Logger
        LOGGER_NAME = "WBIM";
        LOGGER_FILE_THRESHOLD = mlog.Level.DEBUG
        LOGGER_COMMAND_WINDOW_THRESHOLD = mlog.Level.MESSAGE
        LOGGER_EVENT_THRESHOLD = mlog.Level.WARNING;
        LOGGER_TIME_FORMAT = 'yyyy-mm-dd hh:MM:ss.FFF';
        %% Imaging
        IMAGING_MAX_NUM_CHANNEL = 4;
        % Direction of the stage in the image coordiante
        % Image X: parallel to the shorter edge of the table. This grid
        % axis order allows stithcing the acquired tile in the sample space
        % directly 
        IMAGING_2DGRID_AXIS_ORDER = [2,1];
        % Acquisition grid axis order: this order is in the grid
        % coordinate. Therefore, [2,1], combined with the
        % IMAGING_2DGRID_AXIS_ORDER bing [2,1] means move along the X
        % direction, then along the Y axis in the sample coordinate.
        IMAGING_ACQ_AXIS_ORDER = [2,1];
        % (Pockel cell) x (Mirrors and lens) x (Objective)
        IMAGING_LASER_POWER_TRANSMISSION_RATE = 0.51 * 0.7 * 0.8;
        
        IMAGING_SCAN_VOLUME_FLYBACK_TIME_s = 0;
        IMAGING_EXPLORE_VOLUME_FLYBACK_TIME_s = 0;
        % IMAGING_CHANNEL_UNMIXING_MATRIX(j, i) is the fraction of the j-th
        % channel intensity appear in the i-th channel
        % 2023/06/08
        IMAGING_CHANNEL_UNMIXING_MATRIX = [    1, -0.12,    0, -0.12; ...
                                           -0.02,     1,    0, -0.02; ...
                                               0,     0,    1,     0; ...
                                               0,     0,    0,     1];
                                       
        IMAGING_MPPC_STATE_CH1_PORT_NAME = '/vDAQ0/D2.0';
        
        % Resonant scanner fill factor
        IMAGING_RS_NORMINAL_FREQUENCY_HZ = 7910;
        IMAGING_PULSE_PER_PIXEL = 4;
        IMAGING_RS_FILL_FACTRACTION_TIME = 1024 * ...
            WBIMConfig.IMAGING_PULSE_PER_PIXEL * WBIMConfig.IMAGING_LASER_PULSE_SPACING_s ...
            * 2 * WBIMConfig.IMAGING_RS_NORMINAL_FREQUENCY_HZ;
        %% Ablation
        ABLATION_RESONANT_FREQUENCY_Hz = 15.5493;
        ABLATION_ACC_MAX_FREQUENCY_FRACTION = 0.75;
        ABLATION_ACC_MAX_M_S2 = 1.0;
        ABLATION_ACC_LENGTH_EXPANSION_RATIO = 1.05;
        ABLATION_GRATING_SPACING_mm = 1/830;
        ABLATION_DIFFRACTION_STAGE_AXIS = 1;
        ABLATION_LENS_F_mm = 700;
        ABLATION_CYLINDRICAL_LENS_Q = true;
        
        ABLATION_DIFFUSER_FWHM_Deg = 0.25;        
        
        ABLATION_SHUTTER_IDX = 2;
        % 0.0075 second appear to work for length 400, 600 um, but not 1000
        % mm
        ABLATION_ZABER_PRETRIGGER_TIME_s = 0.0075;
        ABLATION_ZABER_PRETRIGGER_TIME_UNCERTAINTY_ms = 5;
        ABLATION_PUMP_PORT_CHANNEL_NAME = '/vDAQ0/D1.1';
        % Diffuser: 
        ABLATION_DIFFUSER_CONTROL_PORT_NAME = '/vDAQ0/D0.6';
        ABLATION_DIFFUSER_STATE_PORT_NAME = '/vDAQ0/D2.7';
        
        ABLATION_MAX_NUM_EXPLORATION_REFINEMENT_TIMES = 3;
        % Power meters
        ABLATION_THERMAL_PM_NAME = 'PM102'
        ABLATION_PHOTODIODE_PM_NAME = 'PM100D'
        %% Motorized rotation stage
        ABLATION_HWP_STAGE_SERIAL_NUMBER = '27255786';
        % 2022/08/29
%         ABLATION_OTH_DATA_FP = './Calibration/AlbWPPwrVsAgl/AlbOthPowVsAng_20220829.mat';
%         ABLATION_HWP_DELTA_THETA_deg = -31.16;        
%         ABLATION_HWP_PMax_W = 5.08;
%         ABLATION_OTH_HWP_PMax_mW = 10.49;        
%         ABLATION_OTH_HWP_DELTA_THETA_deg = 10.07;  
%         ABLATION_TRANSMISSION_FRACTION = 0.66 * 0.88 * 0.8 * 0.95; % 830/mm grating, w/o diffuser, 09/01/2022
        % 2022/09/22
%         ABLATION_OTH_DATA_FP = './Calibration/AlbWPPwrVsAgl/AlbOthPowVsAng_20220922.mat';
%         ABLATION_HWP_DELTA_THETA_deg = -31.1099;
%         ABLATION_HWP_PMax_W = 4.8884;
%         ABLATION_OTH_HWP_PMax_mW = 12.30;            
%         ABLATION_OTH_HWP_DELTA_THETA_deg = 10.08;  
        ABLATION_OTH_DATA_FP = './Calibration/AlbWPPwrVsAgl/AlbOthPowVsAng_20240306.mat';
        ABLATION_HWP_DELTA_THETA_deg = -31.28;
        ABLATION_HWP_PMax_W = 4.8046;
        ABLATION_OTH_HWP_PMax_mW = 27.52;            
        ABLATION_OTH_HWP_DELTA_THETA_deg = 12.87;  
        % Transmission efficiency 
        % (after grating) x objective x (before the grating) 2022/09/22
%         ABLATION_TRANSMISSION_FRACTION = 0.5378 * 0.8 * 0.97;
        % (after grating) x objective x (before the grating) 2023/03/11
%         ABLATION_TRANSMISSION_FRACTION = 0.5867 * 0.8 * 0.977;
        % After replacing Astrella Pockel cells (05/26/2023)
        ABLATION_TRANSMISSION_FRACTION = 0.5867 * 0.8 * 0.969; 
        ABLATION_TRANSMISSION_FRACTION_WITH_DIFFUSER = 0.4683 * 0.8 * 0.969;
        
        ABLATION_HWP_DEFAULT_ENERGY_FRACTION = 0.1;
        ABLATION_DIFFUSER_TRANSMISSION_RATE = 0.9;
        % (Grating + diffuser efficiency) x (After grating cutoff) x (Objective) x (others) 
%         ABLATION_TRANSMISSION_FRACTION = 0.6909 * 0.89 * 0.8 * 0.9; 
%         ABLATION_TRANSMISSION_FRACTION = 0.53 * 0.89 * 0.8 * 0.9; % Grating degraded. Measured 08/2022
        
        IMAGING_HWP_STAGE_SERIAL_NUMBER = '27257631';
        IMAGING_HWP_DELTA_THETA_deg = -16.1932;
        IMAGING_HWP_DEFAULT_ENERGY_FRACTION = 0.5;
        %% Acquisition
        SI_FILESTEM = 'WBIM'; % Has to be char, not string
        SI_IM_TYPE = 'int16';
        SI_CONVERTED_IM_TYPE = 'uint16';
        SI_FRAME_FLYBACK_TIME_s = 2e-3; % 2 ms
        SI_MAX_IDLE_TIME_s = 120;
        SI_PIEZO_FLYBACK_TIME_s = 100e-3; % 50 ms
        FILENAME_TIME_FORMAT = 'YYYYmmDDHHMMss';
        % Estsimated time delay between acqAbort and resuming the
        % acquisition. This number should be smaller than
        % WBIMConfig.SI_MAX_IDLE_TIME_s
        SI_ACQ_RESTART_WAITING_TIME_s = 10;
        SI_LOST_FRAME_RESTART_WAITING_TIME_s = 45;
        
        SI_MAX_REIMAGE_TIME = 1;
        SI_MAX_REIMAGE_TIME_LASER_ABNORMAL = 3;
        SI_MAX_REIMAGE_TIME_ROW_SHIFT = 1;
        %% Post processing
        ACQ_POST_PROCESSING_CLUSTER = 'local';
        ACQ_POST_PROCESSING_TIMEOUT_s = 120;
        %% Pipeline control 
        MACHINE_MAX_IDLE_TIME_s = 600;
        %% Index Correction
        RI_LCD_ON_PORT_NAME = '/vDAQ0/D1.7';
        RI_BRIX_READ_PORT_NAME = '/vDAQ0/D3.6';
        RI_RED_LED_PORT_NAME = '/vDAQ0/D1.6';
        RI_LCD_STATE_PORT_NAME = '/vDAQ0/D0.7';
        RI_VALVE_PORT_NAME = '/vDAQ0/D3.7';
        
        RI_LOGGER_NANME = 'WBIMRI';
        RI_SET_POINT = 1.4286;
        RI_WAIT_TIME_s = 600;
        RI_CORRECTION_PERIOD_s = 60;
        RI_VALVE_V_5s_mL = 8.25; % Sincheng's magic number
        RI_TOTAL_VOLUME_mL = 1500;
    end
    %% Handles
    properties(Hidden)
        T_stage_to_sample CoordinateTransformation
    end    
    %% Dependent properties
    properties(Hidden, Dependent)
       Z_STAGE_MAX_HEIGHT_mm         
    end    
    
    methods
        function val = get.Z_STAGE_MAX_HEIGHT_mm(obj)
            val = obj.DIST_BREADBOARD_BUTTON_TO_TABLE_mm - ...
                obj.OBJECTIVE_DIST_TO_BREADBOARD_BUTTON_mm - ...
                obj.OBJECTIVE_WORKING_DISTANCE_mm - ...
                obj.STAGE_MINIMUM_HEIGHT_mm - ...
                obj.SAMPLE_ADAPTER_PLATE_THICKNESS_mm - ...
                obj.SAMPLE_HOLDER_THICKNESS_mm;
            
            val = min(val, obj.STAGE_LIMIT_um(3)/1e3);
        end        
    end
    %%
    methods
        function obj = WBIMConfig()
            obj.T_stage_to_sample = CoordinateTransformation(obj.A_STAGE_TO_SAMPLE);
        end
    end
%% Coordinate transform
    methods
        function xyz_hstack = sample_to_stage_xyz_um(obj, xyz_hstack)
            if isvector(xyz_hstack)
                if isrow(xyz_hstack)
                    xyz_hstack = xyz_hstack';
                end
            else
                assert(size(xyz_hstack, 1) == 3);
            end
            xyz_hstack = obj.T_stage_to_sample.forward_transform(xyz_hstack, false);
        end
        
        function xyz_hstack = stage_to_sample_xyz_um(obj, xyz_hstack)
            if isvector(xyz_hstack)
                if isrow(xyz_hstack)
                    xyz_hstack = xyz_hstack';
                end
            else
                assert(size(xyz_hstack, 1) == 3);
            end
           xyz_hstack = obj.T_stage_to_sample.inverse_transform(xyz_hstack, false); 
        end 
        
        function stage_bbox_xyz = sample_to_stage_bbox_xyz_um(obj, sample_bbox_xyz_mmxx_um)
            sample_bbox_xyz_mm = sample_bbox_xyz_mmxx_um(1:3);
            sample_bbox_xyz_xx = sample_bbox_xyz_mmxx_um(4:6);
            stage_bbox_xyz_mm = obj.sample_to_stage_xyz_um(sample_bbox_xyz_mm).';
            stage_bbox_xyz_xx = obj.sample_to_stage_xyz_um(sample_bbox_xyz_xx).';
            for i = 1 : 3
                if obj.A_STAGE_TO_SAMPLE(i, i) < 0
                    tmp = stage_bbox_xyz_xx(i);
                    stage_bbox_xyz_xx(i) = stage_bbox_xyz_mm(i);
                    stage_bbox_xyz_mm(i) = tmp;
                end
            end
            assert(all(stage_bbox_xyz_mm <= stage_bbox_xyz_xx));
            stage_bbox_xyz = cat(2, stage_bbox_xyz_mm, stage_bbox_xyz_xx);
        end
        
        function im2D_stage = sample_to_stage_im_yx_um(obj, im_sample_um)
            % Transform the image in sample coordinate [y, x] to stage
            % coordiante [y,x]
            % Might need to 0 the last column for local transform
            % Only rotation is needed. One transpose for inverting the
            % rotation matrix, another transpose for converting to MATLAB's
            % def of T. Need to double check if the direction of the axis
            % is flipped instead of rotated. 
            num_sec = size(im_sample_um, 3);
            im2D_stage = cell(num_sec, 1);
            for i = 1 : num_sec
                to_stage_2D = affine2d(obj.T_stage_to_sample.A_xy_local);
                im2D_stage{i} = imwarp(im_sample_um(:, :, i), to_stage_2D);
            end
            if num_sec == 1
                im2D_stage = im2D_stage{1};
            else
               im2D_stage = cat(3, im2D_stage{:});
            end
        end
    end    
    %% Configuration space of the stage
    methods
        function feasibleQ = stage_xyz_um_is_in_c_space_Q(obj, xyz_um)
            % Check if a position in the stage coordiante is in the
            % configuration space
            % Input: 
            %   xyz_um: 3-element row vecotr or N-by-3 matrix
            % 
            % TODO: replace the conditions by a 3D volumetric
            % reprensentation of the configuration space. 
            if iscolumn(xyz_um)
                xyz_um = xyz_um.';
            else
                assert(size(xyz_um, 2) == 3, 'xyz_um should have 3 columns')
            end
            feasibleQ = xyz_um(:, 1) > obj.STAGE_X_MIN_um & ...
                xyz_um(:, 1) < obj.STAGE_X_MAX_um & ...
                xyz_um(:, 2) > obj.STAGE_Y_MIN_um & ...
                xyz_um(:, 2) < obj.STAGE_Y_MAX_um & ...
                xyz_um(:, 3) < obj.STAGE_Z_MAX_um;
        end
        
        function feasibleQ = sample_xyz_um_is_in_c_space_Q(obj, xyz_um)
            % Check if a position in the sample coordiante is in the
            % configuration space
            % Input: 
            %   xyz_um: 3-element row vecotr or N-by-3 matrix
            if iscolumnvector(xyz_um)
                xyz_um = xyz_um.';
            else
                assert(size(xyz_um, 2) == 3, 'xyz_um should have 3 columns')
            end
            xyz_um_stage = obj.sample_to_stage_xyz_um(xyz_um.').';
            feasibleQ = obj.stage_xyz_um_is_in_c_space_Q(xyz_um_stage);            
        end
        
    end
    %%
    methods(Access=protected)
        function result = is_valid_z_pos_Q(obj, pos_um)
            validateattributes(pos_um, {'numeric'}, {'nonnegative', 'scalar'});
            result = (obj.Z_STAGE_MAX_HEIGHT_mm >= pos_um);
        end        
        
        function result = is_valid_actuator_pos_Q(obj, pos_um)
            validateattributes(pos_um, {'numeric'}, {'scalar', 'nonnegative'});
            result = (pos_um <= obj.ACTUATOR_POS_MAX_um);
        end
    end

end