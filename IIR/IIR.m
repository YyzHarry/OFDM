%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  IIR滤波器设计，IIR数字滤波器从源信号中提取特定频率的信号                       %%
%  手动录制：以 8kHz 采样频率录制自己的声音，叠加高斯白噪声，信噪比小于20dB        %%
%                                                                             %%
%  完成两种低通滤波器设计：脉冲响应不变法 + 切比雪夫I型； 双线性变换法 + 巴特沃斯   %%
%                                                                             %%
%                                   Yuzhe Yang                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%------------------------------------- 录音部分 ----------------------------------%
% 
% 运行平台：Windows 10 64bit，MATLAB R2014a
% 录音5秒钟，采样频率默认为 8kHz
%{
recObj = audiorecorder;
disp('Start speaking.')
recordblocking(recObj, 5);
disp('End of Recording.');
play(recObj);       % 回放录音数据
myRecording = getaudiodata(recObj);         % 获取录音数据

plot(myRecording);          % 绘制录音数据时域波形
plot(abs(fftshift(fft(myRecording))));      % 绘制录音数据频谱
%}


%--------------------------------- 数据处理部分 ----------------------------------%
load myrecord.mat myRecording           % 导出已录制好的音频

%------------------------------- 做FFT，画出原信号幅度谱
N1 = 40000;     %  采样点数，时间长度 5s
fs = 8000;      %  采样频率 8kHz
n1 = 0:N1-1;
t1 = n1/fs;     %  时间序列

myrecord_F = fftshift(fft(myRecording,N1));    %  对信号进行FFT,并转移到中央
mag1 = abs(myrecord_F);             %  幅度谱
F1 = n1*fs/N1 - fs/2;       %  频谱序列长度，且将零频移到中心

figure
plot(t1,myRecording);          %  绘制录音数据时域波形
title('源信号时域波形');
figure
plot(F1, mag1, '-r');          %  绘制录音数据频谱，并局限至 2kHz内（人声音频谱范围85~1100Hz）
axis([-2000 2000 0 1500]), xlabel('f/Hz'), title('源信号频谱');

%------------------------------- 对原信号加高斯白噪声，信噪比5dB（<20dB）
SNR = 5;        %  信噪比 5dB
myrecord_noise = awgn(myRecording,SNR,'measured');
mag2 = abs(fftshift(fft(myrecord_noise,N1)));

figure
plot(t1,myrecord_noise);          %  绘制加噪声后时域波形
title('加噪声后信号时域波形');
figure
plot(F1, mag2, '-r');             %  绘制加噪声后数据频谱，并局限至 2kHz内
axis([-2000 2000 0 1500]), xlabel('f/Hz'), title('加噪声后信号频谱');


%------------------------------- 方法1：脉冲响应不变法 设计 低通切比雪夫I型 -------------------------------%
%
%   滤波器参数：通带上限300Hz, 阻带临界500Hz，抽样频率8kHz
%              通带最大衰减1dB，阻带最小衰减50dB
wp = 2*pi*300;
ws = 2*pi*500;      %  数字滤波器特征转化到模拟滤波器频率特征
Rp = 1;
Rs = 50;
Fs = 8000;
[N, Wn] = cheb1ord(wp,ws,Rp,Rs,'s');        %  选择滤波器的最小阶数
[Z, P, K] = cheb1ap(N,Rp);                  %  创建低通切比雪夫I滤波器
[A, B, C, D] = zp2ss(Z,P,K);                %  归一化低通滤波器
[At, Bt, Ct, Dt] = lp2lp(A,B,C,D,Wn);       %  去归一化，得到截止频率 Wn 的
[num_ana,den_ana] = ss2tf(At,Bt,Ct,Dt);     %  得到模拟滤波器传递函数（分子、分母多项式系数）

[num_dig,den_dig] = impinvar(num_ana,den_ana,Fs);      %  脉冲响应不变法 转换为数字滤波器传递函数（分子、分母多项式系数）
disp(num_dig)
disp(den_dig)                %  IIR系统函数零极点
figure
[H, W] = freqz(num_dig,den_dig,N1,'whole');            %  滤波器幅频响应
plot(W*Fs/2/pi,20*log10(abs(H)));
axis([0 600 -80 0])
grid; xlabel('f/Hz'); title('方法1：脉冲响应不变法，切比雪夫I型低通滤波器');
figure
subplot(121)
plot(W,angle(H)); xlabel('\omega'); ylabel('\phi'); title('相位响应'); axis([0,pi,-pi,pi]);      %  相频响应和系统零极图
subplot(122)
zplane(num_dig,den_dig); axis([-1.1,1.1,-1.1,1.1]); title('零级图');

%---------------------------------- 加噪信号经过此滤波器
out_mag = abs(mag2.*fftshift(H));
figure
plot(F1, out_mag, '-r');             %  绘制滤波后数据频谱，并局限至 2kHz内
axis([-2000 2000 0 1500]), xlabel('f/Hz'), title('滤波后信号频谱');


%------------------------------- 方法2：双线性变换法 设计 低通巴特沃斯 --------------------------------%
%
%   滤波器参数：通带上限300Hz, 阻带临界500Hz，抽样频率8kHz
%              通带最大衰减1dB，阻带最小衰减50dB
Wp = 0.075*pi;            %  通带边界频率(归一化)：wp=fp*2*pi/fs
Ws = 0.125*pi;            %  阻带边界频率(归一化)：ws=fr*2*pi/fs
Rp = 1;                   %  通带波纹
Rs = 50;                  %  阻带衰减
Ts = 0.000125;            %  Fs = 8kHz
%   转换为模拟滤波器指标（预畸变）
OmegaP = (2/Ts)*tan(Wp/2);        %  模拟低通原型滤波器通带频率
OmegaS = (2/Ts)*tan(Ws/2);        %  模拟低通原型滤波器阻带频率

[N, OmegaC] = buttord(OmegaP,OmegaS,Rp,Rs,'s');     %  模拟巴特沃斯滤波器的阶数和-3dB截止频率计算

[Z, P, K] = buttap(N);                              %  设计归一化巴特沃兹滤波器低通原型
num_ana = K * real(poly(Z));
den_ana = real(poly(P));                            %  原型（归一化）滤波器系数
[num_ana, den_ana] = lp2lp(num_ana,den_ana,OmegaC); %  去归一化处理，得到截止频率 OmegaC 的

[num_dig, den_dig] = bilinear(num_ana,den_ana,Fs);  %  双线性变换法 转换为数字滤波器传递函数
disp(num_dig)
disp(den_dig)                %  IIR系统函数零极点
figure
[H, W] = freqz(num_dig,den_dig,N1,'whole');         %  滤波器幅频响应
plot(W*Fs/2/pi,20*log10(abs(H)));
axis([0 600 -80 0])
grid; xlabel('f/Hz'); title('方法2：双线性变换法，巴特沃斯低通滤波器');

%---------------------------------- 加噪信号经过此滤波器
out_mag = abs(mag2.*fftshift(H));
figure
plot(F1, out_mag, '-r');             %  绘制滤波后数据频谱，并局限至 2kHz内
axis([-2000 2000 0 1500]), xlabel('f/Hz'), title('滤波后信号频谱');
figure
subplot(121)
plot(W,angle(H)); xlabel('\omega'); ylabel('\phi'); title('相位响应'); axis([0,pi,-pi,pi]);      %  相频响应和系统零极图
subplot(122)
zplane(num_dig,den_dig); axis([-1.1,1.1,-1.1,1.1]); title('零级图');
