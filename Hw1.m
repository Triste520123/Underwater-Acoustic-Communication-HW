clear;close all;clc;
communication_distance_calculator()
function communication_distance_calculator()
    % 水下声通信理论距离计算（参数化版本）
    % 作者: 基于范嘉晨和朱江的文档

    % 默认参数设置（可修改）
    SL = 170;        % 声源级 (dB)
    SNR_req = 20;    % 所需接收信噪比 (dB)
    f_min = 100;     % 最低频率 (Hz)
    f_max = 200;     % 最高频率 (Hz)
    f_step = 10;     % 频率步长 (Hz)
    s_values = 0;   % 航运活动程度
    w_values = 1; % 风速/波高值
    BW = 100;        % 带宽 (Hz)
    gamma_values = [1, 1.5, 2]; % 扩展损失系数
    
    % 用户输入参数
    prompt = {'声源级 SL (dB):', '所需信噪比 SNR (dB):', ...
              '最低频率 f_min (Hz):', '最高频率 f_max (Hz):', ...
              '频率步长 f_step (Hz):', '带宽 BW (Hz):', ...
              '航运活动程度 (0-1，支持输入区间如0：0.5：1):', ...
              '风速/波高值 (0-10，支持输入区间如0：5：10):'};
    dlgtitle = '输入参数';
    dims = [1 50];
    definput = {num2str(SL), num2str(SNR_req), num2str(f_min), ...
                num2str(f_max), num2str(f_step), num2str(BW), ...
                '0', '1'};
    answer = inputdlg(prompt, dlgtitle, dims, definput);
    
    if isempty(answer)
        disp('用户取消输入，使用默认参数。');
    else
        SL = str2double(answer{1});
        SNR_req = str2double(answer{2});
        f_min = str2double(answer{3});
        f_max = str2double(answer{4});
        f_step = str2double(answer{5});
        BW = str2double(answer{6});
        s_values = str2num(answer{7});
        w_values = str2num(answer{8});
    end
    
    freqs = f_min:f_step:f_max; % 频率范围 (Hz)
    
    % 初始化结果存储
    results = cell(length(s_values), length(w_values), length(gamma_values));
    
    % 计算每个情况和频率下的通信距离
    for s_idx = 1:length(s_values)
        s = s_values(s_idx);
        
        for w_idx = 1:length(w_values)
            w = w_values(w_idx);
            
            for g_idx = 1:length(gamma_values)
                gamma = gamma_values(g_idx);
                
                % 为每个频率计算距离
                distances = zeros(size(freqs));
                for f_idx = 1:length(freqs)
                    f = freqs(f_idx);
                    distances(f_idx) = calculate_distance(SL, SNR_req, f, s, w,gamma);
                end
                
                % 存储结果
                results{s_idx, w_idx, g_idx} = distances;
            end
        end
    end
    
    % 绘制结果
    plot_results(freqs, s_values, w_values, gamma_values, results, SL, SNR_req, BW);
end

function d = calculate_distance(SL, SNR_req, f, s, w, gamma)
    % 计算特定条件下的通信距离
    
    % 转换为kHz用于公式计算
    f_kHz = f / 1000;
    
    % 计算吸收系数 (dB/m)
    term1 = (0.11 * f_kHz^2) / (1 + f_kHz^2);
    term2 = (44 * f_kHz^2) / (4100 + f_kHz^2);
    term3 = 2.75e-4 * f_kHz^2;
    alpha = (term1 + term2 + term3 + 0.003) / 1000;
    
    % 计算噪声级
    f_integrate = 100:1:200; % 1Hz步长
    %f_integrate = 20e3; % 1Hz步长
    total_noise_power = 0;
    
    for f_i = f_integrate
        f_kHz_i = f_i / 1000;
        
        % 计算每个频率点的噪声谱密度（线性值）
        N_turb_i = 10^(1.7 - 3*log10(f_kHz_i));
        N_ship_i = 10^(4 + 2*(s-0.5) + 2.6*log10(f_kHz_i) - 6*log10(f_kHz_i + 0.03));
        N_wave_i = 10^(5 + 0.75*sqrt(w) + 2*log10(f_kHz_i) - 4*log10(f_kHz_i + 0.4));
        N_ther_i = 10^(-1.5 + 2*log10(f_kHz_i));
        
        % 总噪声谱密度（线性值）
        N_total_psd = N_turb_i + N_ship_i + N_wave_i + N_ther_i;
        
        % 积分：噪声谱密度乘以频率步长
        total_noise_power = total_noise_power + N_total_psd *1; % 1Hz步长
    end
    
    % 将积分后的总噪声功率转换为分贝值
    NL = 10*log10(total_noise_power); 
    % 使用数值方法求解距离d
    % 方程: SL - [gamma*10*log10(d) + d*alpha] - NL = SNR_req
    % 重新排列: gamma*10*log10(d) + d*alpha = SL - NL - SNR_req
    
    target = SL - NL - SNR_req;
    
    % 使用fzero求解方程
    options = optimset('Display', 'off');
try
    % 先尝试小范围
    d = fzero(@(d) propagation_loss(d, gamma, alpha) - target, [1, 10e4], options);
catch
    try
        % 失败后尝试更大范围
        d = fzero(@(d) propagation_loss(d, gamma, alpha) - target, [1, 1e7], options);
    catch
        d = NaN;
    end
end
end

function loss = propagation_loss(d, gamma, alpha)
    % 计算传播损失
    loss = gamma * 10 * log10(d) + d * alpha;
end

function plot_results(freqs, s_values, w_values, gamma_values, results, SL, SNR_req, BW)
    % 绘制结果图表 - 优化版
    
    % 使用更丰富的颜色和标记
    colors = lines(7); % 使用MATLAB内置的lines颜色图
    markers = {'o', 's', 'd', '^', 'v', '<', '>', 'p', 'h', '*'};
    gamma_names = {'柱面扩展', '浅海声传播', '球面扩展'};
    
    % 创建图表 - 调整大小和位置
    fig = figure('Position', [100, 100, 1000, 600], 'Name', '水下声通信理论距离分析', 'Color', 'w');
    
    % 使用tiledlayout创建更紧凑的布局
    t = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % 为每种传播模型创建子图
    for g_idx = 1:length(gamma_values)
        ax = nexttile;
        hold(ax, 'on');
        grid(ax, 'on');
        box(ax, 'on');
        
        % 为每个航运活动程度和风速组合创建曲线
        line_count = 1;
        for s_idx = 1:length(s_values)
            for w_idx = 1:length(w_values)
                distances = results{s_idx, w_idx, g_idx};
                if all(~isnan(distances))
                    % 使用更清晰的绘图样式
                    color_idx = mod(line_count-1, size(colors, 1)) + 1;
                    marker_idx = mod(line_count-1, length(markers)) + 1;
                    
                    plot(ax, freqs, distances/1000, ...
                        'Color', colors(color_idx, :), ...
                        'Marker', markers{marker_idx}, ...
                        'LineStyle', '-', ...
                        'LineWidth', 1.5, ...
                        'MarkerSize', 6, ...
                        'MarkerFaceColor', colors(color_idx, :), ...
                        'DisplayName', sprintf('s=%.1f, w=%d', s_values(s_idx), w_values(w_idx)));
                    
                    line_count = line_count + 1;
                end
            end
        end
        
        title(ax, sprintf('传播模型: %s', gamma_names{g_idx}), 'FontSize', 12);
        xlabel(ax, '频率 (Hz)', 'FontSize', 10);
        ylabel(ax, '通信距离 (km)', 'FontSize', 10);
        set(ax, 'YScale', 'log');
        
        % 添加图例 - 放在外部以避免遮挡
        if g_idx == 1
            lg = legend(ax, 'Location', 'eastoutside');
            lg.FontSize = 8;
            lg.Box = 'off';
        end
    end
    
    % 添加参数表格 - 使用更紧凑的显示方式
    ax = nexttile;
    axis(ax, 'off');
    
    % 创建参数文本 - 使用更紧凑的格式
    param_text = {...
        ['声源级 SL: ', num2str(SL), ' dB'], ...
        ['所需信噪比 SNR: ', num2str(SNR_req), ' dB'], ...
        ['带宽 BW: ', num2str(BW), ' Hz'], ...
        ['频率范围: ', num2str(min(freqs)), '-', num2str(max(freqs)), ' Hz'], ...
        '', ...
        '传播模型说明:', ...
        '  γ=1: 柱面扩展', ...
        '  γ=1.5: 浅海声传播', ...
        '  γ=2: 球面扩展', ...
        '', ...
        '航运活动程度 s:', ...
        '  0-低, 0.5-中等, 1-高', ...
        '', ...
        '风速/波高 w:', ...
        '  0-平静, 5-中等, 10-大风浪'};
    
    text(ax, 0.05, 0.95, param_text, ...
        'VerticalAlignment', 'top', ...
        'FontSize', 9, ...
        'Interpreter', 'none');
    
    title(ax, '参数设置', 'FontSize', 10);
    
    % 添加总标题
    title(t, '水下声通信理论距离计算', 'FontSize', 14, 'FontWeight', 'bold');
    
    % 调整布局
    t.TileSpacing = 'tight';
    t.Padding = 'compact';
    
    % 保存图表
    saveas(fig, 'communication_distance_analysis.png');
    
    % 保存数据到Excel
    save_to_excel(freqs, s_values, w_values, gamma_values, results, SL, SNR_req, BW);
end

function save_to_excel(freqs, s_values, w_values, gamma_values, results, SL, SNR_req, BW)
    % 保存结果到Excel文件
    
    % 创建数据表格
    data = {};
    headers = {'频率 (Hz)', '航运活动 s', '风速/波高 w', '传播模型 γ', '通信距离 (km)'};
    data{1, 1} = headers{1};
    data{1, 2} = headers{2};
    data{1, 3} = headers{3};
    data{1, 4} = headers{4};
    data{1, 5} = headers{5};
    
    row = 2;
    for f_idx = 1:length(freqs)
        f = freqs(f_idx);
        for s_idx = 1:length(s_values)
            s = s_values(s_idx);
            for w_idx = 1:length(w_values)
                w = w_values(w_idx);
                for g_idx = 1:length(gamma_values)
                    gamma = gamma_values(g_idx);
                    distances = results{s_idx, w_idx, g_idx};
                    distance_km = distances(f_idx) / 1000;
                    
                    data{row, 1} = f;
                    data{row, 2} = s;
                    data{row, 3} = w;
                    data{row, 4} = gamma;
                    data{row, 5} = distance_km;
                    
                    row = row + 1;
                end
            end
        end
    end
    
    % 写入Excel - 使用兼容性更好的方法
    filename = 'communication_distance_results.xlsx';
    
    % 检查是否有写入Excel的功能
    if ispc && exist('xlswrite', 'file')
        % 写入数据到第一个工作表
        xlswrite(filename, data, '数据结果');
        
        % 准备参数信息
        params = {...
            '参数名称', '值', '单位/说明';...
            '声源级 SL', SL, 'dB';...
            '所需信噪比 SNR', SNR_req, 'dB';...
            '带宽 BW', BW, 'Hz';...
            '频率范围', sprintf('%d-%d', min(freqs), max(freqs)), 'Hz';...
            '航运活动程度 s', '0-低, 0.5-中等, 1-高', '';...
            '风速/波高 w', '0-平静, 5-中等, 10-大风浪', '';...
            '传播模型 γ', '1-柱面, 1.5-浅海, 2-球面', ''};
        
        % 写入参数到第二个工作表
        xlswrite(filename, params, '参数设置');
    else
        % 如果没有Excel支持，保存为CSV文件
        filename = 'communication_distance_results.csv';
        fid = fopen(filename, 'w');
        
        % 写入数据
        for i = 1:size(data, 1)
            for j = 1:size(data, 2)
                if j == size(data, 2)
                    fprintf(fid, '%s', num2str(data{i, j}));
                else
                    fprintf(fid, '%s,', num2str(data{i, j}));
                end
            end
            fprintf(fid, '\n');
        end
        fclose(fid);
        
        % 保存参数到另一个文件
        param_filename = 'communication_distance_parameters.txt';
        fid = fopen(param_filename, 'w');
        fprintf(fid, '参数设置\n');
        fprintf(fid, '声源级 SL: %d dB\n', SL);
        fprintf(fid, '所需信噪比 SNR: %d dB\n', SNR_req);
        fprintf(fid, '带宽 BW: %d Hz\n', BW);
        fprintf(fid, '频率范围: %d-%d Hz\n', min(freqs), max(freqs));
        fclose(fid);
    end
    
    disp(['结果已保存到文件: ' filename]);
end