% verify_RZ_AMI_PSD.m
% 验证：单极性 RZ 与 双极性 AMI 的 PSD（仿真 vs 理论）
clear; close all; clc;

%% ----- 参数 -----
A   = 1;           % 脉冲幅度
Tb  = 1e-3;        % 比特周期 (s)
Fs  = 500e3;       % 采样率 (Hz) - 要比 1/Tb 高很多
Nsym = 2e5;        % 比特数（尽量大以减小统计误差）

duty = 0.5;        % RZ 占空比 50%
Nps = round(Tb * Fs);    % 每符号采样点数
if Nps < 6
    error('Fs 太小，导致每符号采样点数过少，请增大 Fs 或增大 Tb。');
end
pwidth = floor(duty * Nps);  % 脉冲宽度（采样点）

%% ----- 生成比特流 -----
bits = randi([0,1], 1, Nsym);  % 0/1 等概率独立

%% ----- 生成单极性 RZ 波形 -----
sig_rz = zeros(1, Nsym * Nps);
for k = 1:Nsym
    if bits(k) == 1
        idx = (k-1)*Nps + (1:pwidth);
        sig_rz(idx) = A;
    end
end

%% ----- 生成双极性 AMI 波形 -----
sig_ami = zeros(1, Nsym * Nps);
next_polar = A;  % 下一个“1”脉冲的极性
for k = 1:Nsym
    if bits(k) == 1
        idx = (k-1)*Nps + (1:Nps);
        sig_ami(idx) = next_polar;
        next_polar = -next_polar;  % 翻转极性
    end
end

t = (0:length(sig_rz)-1)/Fs;

%% ----- 理论 PSD（连续部分） -----
% 1) RZ：P_rz(f) = (A^2 * Tb / 16) * sinc^2( f*Tb/2 )
%    (使用 MATLAB 的 sinc，其中 sinc(x) = sin(pi*x)/(pi*x))
Nfft = 2^18;
fvec = linspace(-Fs/2, Fs/2, Nfft);
P_rz_th = (A^2 * Tb / 16) .* ( sinc( fvec * Tb / 2 ).^2 );

% 2) AMI：P_ami(f) = A^2 * Tb * sin^2(pi f Tb) * sinc^2(f Tb)
%    注意 matlab sinc 的定义，所以表达为：
P_ami_th = (A^2 * Tb) .* ( sin(pi * fvec * Tb).^2 ) .* ( sinc( fvec * Tb ).^2 );

%% ----- 仿真 PSD：使用 pwelch -----
% 为了比较连续部分，去掉直流均值（RZ 有直流项，AMI 理论上均值接近 0）
sig_rz_zm  = sig_rz  - mean(sig_rz);   % 去直流
sig_ami_zm = sig_ami - mean(sig_ami);  % AMI 均值应约为 0

window = hamming(16 * Nps);
noverlap = round(length(window)/2);
nfft_pwelch = 8192;

[pxx_rz, f_rz] = pwelch(sig_rz_zm, window, noverlap, nfft_pwelch, Fs, 'centered');
[pxx_ami, f_ami] = pwelch(sig_ami_zm, window, noverlap, nfft_pwelch, Fs, 'centered');

% 将理论谱插值到 pwelch 的频点以便比较
P_rz_th_interp = interp1(fvec, P_rz_th, f_rz, 'linear', 0);
P_ami_th_interp = interp1(fvec, P_ami_th, f_ami, 'linear', 0);

%% ----- 绘图比较 -----
fplot_lim = 1.5e4; % 绘图频带范围，可根据 Tb 调整

figure('Units','normalized','Position',[0.05 0.05 0.9 0.85]);

% RZ 连续部分比较（线性）
subplot(2,2,1);
plot(f_rz, pxx_rz, 'b'); hold on;
plot(f_rz, P_rz_th_interp, 'r--','LineWidth',1.4);
xlim([-fplot_lim fplot_lim]);
xlabel('Frequency (Hz)'); ylabel('PSD (linear)');
title('单极性 RZ：仿真 vs 理论（线性）');
legend('仿真 (pwelch, 去直流)', '理论连续部分');
grid on;

% RZ 比例（dB）
subplot(2,2,2);
plot(f_rz, 10*log10(pxx_rz), 'b'); hold on;
plot(f_rz, 10*log10(P_rz_th_interp), 'r--','LineWidth',1.4);
xlim([-fplot_lim fplot_lim]);
xlabel('Frequency (Hz)'); ylabel('PSD (dB/Hz)');
title('单极性 RZ：仿真 vs 理论（dB）');
legend('仿真', '理论连续部分');
grid on;

% AMI 连续部分比较（线性）
subplot(2,2,3);
plot(f_ami, pxx_ami, 'b'); hold on;
plot(f_ami, P_ami_th_interp, 'r--','LineWidth',1.4);
xlim([-fplot_lim fplot_lim]);
xlabel('Frequency (Hz)'); ylabel('PSD (linear)');
title('双极性 AMI：仿真 vs 理论（线性）');
legend('仿真 (pwelch, 去直流)', '理论');
grid on;

% AMI 比例（dB）
subplot(2,2,4);
plot(f_ami, 10*log10(pxx_ami), 'b'); hold on;
plot(f_ami, 10*log10(P_ami_th_interp), 'r--','LineWidth',1.4);
xlim([-fplot_lim fplot_lim]);
xlabel('Frequency (Hz)'); ylabel('PSD (dB/Hz)');
title('双极性 AMI：仿真 vs 理论（dB）');
legend('仿真', '理论');
grid on;

%% ----- 直流/均值检验 -----
mean_rz = mean(sig_rz);
mean_ami = mean(sig_ami);
fprintf('RZ: mean(sig) = %.6g, mean^2 = %.6g, 理论直流系数 A^2/16 = %.6g\n', ...
    mean_rz, mean_rz^2, A^2/16);
fprintf('AMI: mean(sig) = %.6g (应接近 0)\n', mean_ami);

% 显示 pwelch 含直流时 f=0 附近的能量（可见 RZ 的直流谱峰）
[pxx_rz_all, f_rz_all] = pwelch(sig_rz, window, noverlap, nfft_pwelch, Fs, 'centered');
[~, idx0_rz] = min(abs(f_rz_all));
fprintf('pwelch(含直流) 在 f≈0 估计 PSD = %g (unit^2/Hz)\n', pxx_rz_all(idx0_rz));

%% ----- 说明提示 -----
disp('注意：由于窗口、离散频点与有限样本长度，仿真曲线与理论曲线会有小偏差；增大样本数 Nsym 与 nfft 能改善匹配。');