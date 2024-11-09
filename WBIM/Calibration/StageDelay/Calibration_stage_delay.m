exp_root_folder = 'D:\My Drive\Microscope\Calibration\Zaber_controller_delay\Y';
exp_folder_list = dir(sprintf('%s\\Exp*', exp_root_folder));
num_exp = numel(exp_folder_list);
exp_result = cell(num_exp, 1);

for exp_idx = 1 : num_exp
    exp_name = exp_folder_list(exp_idx).name;
    exp_folder = fullfile(exp_root_folder, exp_folder_list(exp_idx).name);
    vis_folder = fullfile(exp_folder, 'vis');
    exp_files = dir(sprintf('%s\\*.csv', exp_folder));
    num_files = numel(exp_files);
    exp_data = cell(num_files, 1);
    for i = 1 : num_files
        tmp_data = readtable(fullfile(exp_files(i).folder, exp_files(i).name));
        exp_data{i} = table2struct(tmp_data, 'ToScalar', true);
    end
    para_file = fullfile(exp_folder, 'Parameters.txt');
    fid = fopen(para_file);
    try
        tmp_c = textscan(fid, '%s %f %s', 'Delimiter', ' ');
        fclose(fid);
    catch
        fclose(fid);
    end
    para_str = struct;
    for i = 1 : numel(tmp_c{1})
        para_str.(tmp_c{1}{i}) = tmp_c{2}(i);
    end
    
    %%
    vis_Q = false;
    result = cell(num_files, 1);
    for file_idx = 1 : num_files        
        tmp_str = struct;
        tmp_str.exp_name = exp_name;
        tmp_str.file_idx = file_idx;
        tmp_str.target_v_m_s = para_str.Maximum_speed * 1e-3;
        tmp_str.target_acc_m_s2 = para_str.Acceleration;
        tmp_str.target_acc_time_s = tmp_str.target_v_m_s / tmp_str.target_acc_m_s2 ;
        tmp_str.target_acc_l_m = tmp_str.target_acc_m_s2 * tmp_str.target_acc_time_s^2/2;
        
        tmp_exp_data = exp_data{file_idx};
        
        tmp_str.l_um = diff(tmp_exp_data.encoder_pos__m_([1, end]));
        tmp_str.dt_s = unique(round(diff(tmp_exp_data.Time_ms_ * 1e-3), 4));
        assert(isscalar(tmp_str.dt_s));
        tmp_str.T_s = tmp_exp_data.Time_ms_(end) * 1e-3;
        
        tmp_exp_data.delay_dist_um = tmp_exp_data.pos__m_ - tmp_exp_data.encoder_pos__m_;
        tmp_exp_data.v_m_s = gradient(tmp_exp_data.pos__m_ * 1e-6) ./ tmp_str.dt_s;
        tmp_exp_data.ev_m_s = gradient(tmp_exp_data.encoder_pos__m_ * 1e-6) ./ tmp_str.dt_s;
        tmp_exp_data.a_m_s2 = gradient(tmp_exp_data.v_m_s) ./ tmp_str.dt_s;
        tmp_exp_data.ea_m_s2 = gradient(tmp_exp_data.ev_m_s) ./ tmp_str.dt_s;
        
        tmp_exp_data.delay_dist_time_s = tmp_exp_data.delay_dist_um * 1e-6 ...
            / tmp_str.target_v_m_s;
        % Encoder acceleration duration
        tmp_exp_data.a_gt_0_idx = find(tmp_exp_data.a_m_s2 > 0.1 * tmp_str.target_acc_m_s2);
        tmp_str.e_pa_idx_range = tmp_exp_data.a_gt_0_idx([1, end]);
        tmp_str.e_pa_t_range_s = tmp_str.e_pa_idx_range * tmp_str.dt_s;
        % Deceleration is less important 
        tmp_exp_data.a_lt_0_idx = find(tmp_exp_data.a_m_s2 < -0.1 * tmp_str.target_acc_m_s2);
        tmp_str.e_na_idx_range = tmp_exp_data.a_lt_0_idx([1,end]);
        tmp_str.e_na_t_range_s = tmp_str.e_na_idx_range * tmp_str.dt_s;
        % Displacement difference at the end of the acceleration phase
        tmp_str.eoa_delay_dist_um = tmp_exp_data.delay_dist_um(tmp_str.e_pa_idx_range(2));
        tmp_str.eoa_delay_dist_um_2_at_s = tmp_str.eoa_delay_dist_um * 1e-6 / ...
            tmp_str.target_acc_m_s2 / tmp_str.e_pa_t_range_s(2);
                
        tmp_str.const_v_t_idx = find(abs(tmp_exp_data.v_m_s - tmp_str.target_v_m_s) < 1e-6);
        if ~isempty(tmp_str.const_v_t_idx)
            tmp_str.const_v_t_s_0 = tmp_exp_data.Time_ms_(tmp_str.const_v_t_idx(1));
        else
            tmp_str.const_v_t_s_0 = [];
        end
        tmp_v_m_s = diff(tmp_exp_data.pos__m_) * 1e-6 ./ ...
            diff(tmp_exp_data.Time_ms_ * 1e-3);
        tmp_str.const_v_disp_diff_um = tmp_exp_data.delay_dist_um(tmp_str.const_v_t_idx);
        tmp_str.disp_diff_um_stat = fun_analysis_get_basic_statistics(tmp_str.const_v_disp_diff_um);
        
        tmp_str.const_v_disp_delay_s = tmp_exp_data.delay_dist_time_s(tmp_str.const_v_t_idx);
        tmp_str.dis_diff_delay_s_stat = fun_analysis_get_basic_statistics(tmp_str.const_v_disp_delay_s);
        result{file_idx} = tmp_str;
        if vis_Q
            fig_hdl = figure;
            fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* [2, 1.5];
            ax_1 = nexttile();
            yyaxis(ax_1, 'left');
            plt_ideal_pos = plot(ax_1, tmp_exp_data.Time_ms_, tmp_exp_data.pos__m_);
            hold(ax_1, 'on');
            plt_encoder_pos = plot(ax_1, tmp_exp_data.Time_ms_, tmp_exp_data.encoder_pos__m_);
            ax_1.YLabel.String = 'Position (\mum)';
            yyaxis(ax_1, 'right');
            plt_diff = plot(ax_1, tmp_exp_data.Time_ms_, tmp_exp_data.delay_dist_um);
            ax_1.YLabel.String = 'Distance (\mum)';
            if tmp_str.e_pa_idx_range(2) < tmp_str.e_na_idx_range(1)-1
                l1 = line(ax_1, tmp_str.e_pa_t_range_s([2,2])*1e3, ...
                    ax_1.YLim, 'Color', 'g');
                l2 = line(ax_1, tmp_str.e_na_t_range_s([1,1])*1e3, ...
                    ax_1.YLim, 'Color', 'g');
            end
            legend(ax_1, [plt_ideal_pos, plt_encoder_pos, plt_diff],...
                'Ideal', 'Encoder', 'Difference delay', 'Location', 'best');
            ax_1.Title.String = sprintf('Experiment %d', file_idx);
            ax_1.XLabel.String = 'Time (ms)';
            
            ax_2 = nexttile();
            yyaxis(ax_2, 'left');
            plot(ax_2, tmp_exp_data.Time_ms_, tmp_exp_data.ev_m_s);
            hold(ax_2, 'on');
            plot(ax_2, tmp_exp_data.Time_ms_, tmp_exp_data.v_m_s);
            ax_2.XLabel.String = 'Time (ms)';
            ax_2.YLabel.String = 'Speed (m/s)';
            yyaxis(ax_2, 'right');
            plot(ax_2, tmp_exp_data.Time_ms_, tmp_exp_data.a_m_s2);
            ax_2.YLabel.String = 'Acceleration (m/s^2)';
            legend(ax_2, 'Encoder', 'Ideal', 'Ideal Acc', 'Location', 'best');
            
            if ~isempty(tmp_str.const_v_disp_diff_um)
                ax_3 = nexttile();
                histogram(ax_3, tmp_str.const_v_disp_diff_um);
                ax_3.XLabel.String = 'Displacement (\mum)';
                ax_3.YLabel.String = 'Count';
                tmp_disp_str = fun_analysis_basic_stat_str_to_string(tmp_str.disp_diff_um_stat);
                legend(ax_3, tmp_disp_str, 'Location', 'best');
                ax_4 = nexttile();
                histogram(ax_4, tmp_str.const_v_disp_delay_s);
                ax_4.XLabel.String = 'Time delay (s)';
                ax_4.YLabel.String = 'Count';
                tmp_disp_str = fun_analysis_basic_stat_str_to_string(tmp_str.dis_diff_delay_s_stat);
                legend(ax_4, tmp_disp_str, 'Location', 'best');
            end
            fig_fp = fullfile(vis_folder, sprintf('Stage_delay_file_%d.png', file_idx));
            fun_print_image_in_several_formats(fig_hdl, fig_fp);
        end
    end
    %%
    result = cat(1, result{:});
    exp_result{exp_idx} = result;
end
%%
% fig_hdl = figure;
% ax_1 = nexttile();
% plot(ax_1, 
%% Questions
% 1. Is the delay time stable w.r.t. different acceleration & target speed
% 2. How does the displacement distance depends on the acceleration
% duration? 
target_v_list = cellfun(@(x) x(end).target_v_m_s, exp_result(3:5));
target_a_list = cellfun(@(x) x(end).target_acc_m_s2, exp_result(3:5));
const_v_delay_t_list = cellfun(@(x) x(end).dis_diff_delay_s_stat.median, exp_result(3:5));
eoa_delay_dist_um_list = cellfun(@(x) x(end).eoa_delay_dist_um_2_at_s, exp_result(3:5));

