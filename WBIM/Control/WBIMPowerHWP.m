classdef WBIMPowerHWP < ThorlabsPRM1Z8
    properties(Hidden, Constant)
        COMPONENT_NAME = 'WBIMPowerHWP'
    end
    
    properties
        theta_0_deg (1,1) double
        
    end
    
    properties(Dependent)
        theta_deg (1,1) double
        fractional_power (1,1) double
    end
    %%
    methods
        function obj = WBIMPowerHWP(serial_num, coeff_delta_theta, default_fraction)            
            obj@ThorlabsPRM1Z8(serial_num);
            obj.theta_0_deg = coeff_delta_theta;
            if nargin == 3
                try
                   obj.set_fraction_power(default_fraction); 
                catch ME
                   fprintf("Fail to set the default fractional power\n");
                   fprintf("%s", getReport(ME, 'extended'));
                end
            end           
        end
        
        function delete(obj)
            try
                delete@ThorlabsPRM1Z8(obj);
            catch ME
                delete@ThorlabsPRM1Z8(obj);
                rethrow(ME);
            end            
        end
    end
    %% Convert betweem fractional energy and angle
    methods
        function val = get_fractional_power_at_angle(obj, theta_deg)
            val = cosd(2 * (theta_deg + obj.theta_0_deg)) .^ 2;
        end
        
        function theta = get_angle_at_fractional_power(obj, pwr_r)
            validateattributes(pwr_r, {'float', 'double', 'single'}, ...
                {'nonnegative', '<=', 1});
            % Force the range of theta to be around 0 degree.
            theta = - 1/2 * (acosd(sqrt(pwr_r))) - obj.theta_0_deg;
            theta = mod(theta, 360);
        end
        
        function obj = set_fraction_power(obj, pwrf)
            validateattributes(pwrf, {'float', 'double', 'single'}, ...
                {'nonnegative', '<=', 1});
            new_theta_deg = obj.get_angle_at_fractional_power(pwrf);
            delta_theta_deg = new_theta_deg - obj.theta_deg;
            if delta_theta_deg > 180
                delta_theta_deg = delta_theta_deg - 360;
            end
            obj.moveRelative(delta_theta_deg);
%             obj.moveTo(new_theta_deg);  
        end
        
        function val = get_othogonal_path_angle(obj, abl_agl)
            if nargin < 2
                abl_agl = obj.theta_deg;
            end
            theta_0_oth = obj.theta_0_deg + 45;
            agl = mod(abl_agl + theta_0_oth, 90);
            if agl > 45
                agl = agl - 45;
            end
            val = agl - theta_0_oth;
        end
    end
    %% 
    methods
        function val = get.theta_deg(obj)
            if ~isempty(obj.device) && isvalid(obj.device)
               val = obj.getPosition(); 
            else
                val = nan;
            end
        end
        
        function val = get.fractional_power(obj)
           val = obj.get_fractional_power_at_angle(obj.theta_deg); 
        end
    end
end