%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  对OFDM系统进行仿真，得到误比特率曲线                                                          %%
%  采用BPSK进行调制，系统对128个BPSK调制符号进行IFFT运算                                         %%
%  IFFT运算结果构成一个OFDM符号系统对128个BPSK调制符号进行IFFT运算，IFFT运算结果构成一个OFDM符号   %%
%  收发端：FIR滤波器，等波纹滤波器                                                              %%
%  信道：AWGN信道                                                                              %%
%                                                                                             %%
%                                           Yuzhe Yang                                        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

%---------------------------------------- 系统参数定义 ----------------------------------------%
nFFT = 128;             % n点fft的长度
nBpsk_Ofdm = 128;       % 每个OFDM符号带有的BPSK符号
nBit_Sym = 128;         % 每个OFDM符号所带的bit数目（对于BPSK而言和 nBpsk_Ofdm 一样）
nSym = 10^4;            % OFDM 符号数目（bit number > 10^5，即symbol > 10^3，此处取10^4保证数目精确度）

EbN0dB = 0:0.5:12;      % Bit to Noise Ratio
EsN0dB = EbN0dB + 10*log10(nBpsk_Ofdm/nFFT) - 10*log10(8);         % 转换为每个symbol的SNR, Nsample = 8


%---------------------------------------- 滤波器设计 及 其幅相频响应 ----------------------------------------%
[n,f,a,w] = firpmord([3.5*10^6 4.5*10^6], [1, 0], [0.005, 0.001], 24*10^6);
Hd = firpm(n,f,a,w);                        % 设计滤波器的系数
disp(Hd)                                    % FIR系统函数零点的显示

freqz(Hd, 1, 'whole');                      % 滤波器幅频\相频响应
figure
zplane(Hd, 1); axis([-1.1,1.1,-1.1,1.1]); title('零极点图');        % 滤波器的零极点图
grpdelay(Hd);  title('群延迟图');                                   % 滤波器的群延迟
stem(Hd); title('滤波器时域响应/冲激响应');                          % 滤波器时域响应

%---------------------------------------- OFDM 系统误比特率仿真 ----------------------------------------%
for ii = 1:length(EbN0dB)

    % 发射端
    s_bit = round(rand(1,nBit_Sym*nSym));               % 随机0, 1比特序列
    s_bit_bpsk = 2*s_bit - 1;                           % BPSK调制： 0 --> -1, 1 --> +1
    s_bit_bpsk = reshape(s_bit_bpsk,nBit_Sym,nSym).';   % 将调制后信号重新分组为 nSym*nBit_Sym 矩阵，便于操作
    % disp(s_bit_bpsk(1:10));
  
    % IFFT变换 实现OFDM每个符号
    % (nFFT/sqrt(nBpsk_Ofdm))是将每个发送符号的能量归一化到1
    s_ofdm = (nFFT/sqrt(nBpsk_Ofdm))*ifft(fftshift(s_bit_bpsk.')).';    % ifft前需要fftshift到中心，且由于为矩阵的shift，需要作转置

    % 插入循环前缀CP，设定其长度为 L=16
    s_cp = [s_ofdm(:,113:128) s_ofdm];

    % 将多个OFDM符号串联到一个一维向量
    s_cp = reshape(s_cp.',1,nSym*(16+128));
    
    % 将实部与虚部信号分别处理
    s_I = real(s_cp);
    s_Q = imag(s_cp);
    
    % 8倍上采样
    s_Iup = upsample(s_I,8);
    s_Qup = upsample(s_Q,8);
    
    s_Iup = reshape(s_Iup, 8*(nBit_Sym+16), nSym).';      % reshape
    s_Qup = reshape(s_Qup, 8*(nBit_Sym+16), nSym).';
    
    %{
    % 经过发射端低通滤波器
    s_I_filter = zeros(nSym,144*8);
    s_Q_filter = zeros(nSym,144*8);
    for j = 1:nSym
        s_I_filter(j,:) = filter(Hd,1,s_Iup(j,:));      % 卷积, conv()函数有问题？？暂时用filter
        s_Q_filter(j,:) = filter(Hd,1,s_Qup(j,:));
    end
    
    s_I_filter = reshape(s_I_filter.', 1, 8*nSym*(16+128));     % reshape
    s_Q_filter = reshape(s_Q_filter.', 1, 8*nSym*(16+128));
    %}

    % 经过发射端低通滤波器
    % Hd = my_filter;
    s_I_filter = zeros(nSym,144*8 + n);
    s_Q_filter = zeros(nSym,144*8 + n);
    for j = 1:nSym
        s_I_filter(j,:) = conv(s_Iup(j,:),Hd);                  % conv卷积
        s_Q_filter(j,:) = conv(s_Qup(j,:),Hd);
    end
    
    % 由于卷积引入了多余项，需要去掉（滤波器阶数 = 点数 = n）
    s_I_filter = s_I_filter(:,(n/2+1):end-n/2);
    s_Q_filter = s_Q_filter(:,(n/2+1):end-n/2);
    
    s_I_filter = reshape(s_I_filter.', 1, 8*nSym*(16+128));     % reshape
    s_Q_filter = reshape(s_Q_filter.', 1, 8*nSym*(16+128));
    
    % 两路正交信号经过LPF后合并
    s_filter = s_I_filter + 1j*s_Q_filter;
    
    % 通过AWGN信道
    s_awgn = awgn(s_filter, EsN0dB(ii), 'measured');
    
    % 接收端将实部与虚部信号分别处理
    r_I = real(s_awgn);
    r_Q = imag(s_awgn);
    
    r_I = reshape(r_I, 8*(nBit_Sym+16), nSym).';        % reshape
    r_Q = reshape(r_Q, 8*(nBit_Sym+16), nSym).';
    %{
    % 经过接收端低通滤波器
    r_I_filter = zeros(nSym,144*8);
    r_Q_filter = zeros(nSym,144*8);
    for j = 1:nSym
        r_I_filter(j,:) = filter(Hd,1,r_I(j,:));      % 卷积
        r_Q_filter(j,:) = filter(Hd,1,r_Q(j,:));
    end
    
    r_I_filter = reshape(r_I_filter.', 1, 8*nSym*(16+128));     % reshape
    r_Q_filter = reshape(r_Q_filter.', 1, 8*nSym*(16+128));
    %}
    
    % 经过接收端低通滤波器
    r_I_filter = zeros(nSym,144*8 + n);
    r_Q_filter = zeros(nSym,144*8 + n);
    for j = 1:nSym
        r_I_filter(j,:) = conv(r_I(j,:),Hd);            % conv卷积
        r_Q_filter(j,:) = conv(r_Q(j,:),Hd);
    end
    
    % 由于卷积引入了多余项，需要去掉（滤波器阶数 = 点数 = n）
    r_I_filter = r_I_filter(:,(n/2+1):end-n/2);
    r_Q_filter = r_Q_filter(:,(n/2+1):end-n/2);
    
    r_I_filter = reshape(r_I_filter.', 1, 8*nSym*(16+128));     % reshape
    r_Q_filter = reshape(r_Q_filter.', 1, 8*nSym*(16+128));
    
    % 8倍下采样
    r_Idown = downsample(r_I_filter,8);
    r_Qdown = downsample(r_Q_filter,8);
    
    % 两路正交信号合并
    r_filter = r_Idown + 1j*r_Qdown;
    
    r_cp = reshape(r_filter,144,nSym).';        % 将调制后信号重新分组为 nSym*(nBit_Sym+L) 矩阵，便于操作
    
    % 丢弃循环前缀CP，L=16
    r_fft = r_cp(:,17:end);
    
    % 转换到频域，利用FFT将OFDM符号转换为待解调的信息符号
    r_F = (sqrt(nBpsk_Ofdm)/nFFT)*fftshift(fft(r_fft.')).';

    % BPSK 解调
    % 将 +value --> 1, 将 -value --> -1
    r_demod = real(r_F);
    r_demod(find(r_demod > 0)) = +1;
    r_demod(find(r_demod < 0)) = -1;

    % 将解调后的信息符号转换为比特形式
    r_bit = (r_demod + 1)/2;
    r_bit = reshape(r_bit.',nBit_Sym*nSym,1).';

    % 计算错误bits数目
    nErr(ii) = size(find(r_bit - s_bit),2);
        
end

simBer = nErr/(nSym*nBit_Sym);                          % 仿真误比特率
theoryBer = (1/2)*erfc(sqrt(10.^(EbN0dB/10)));          % 理论误比特率 BPSK

figure
semilogy(EbN0dB,simBer,'--*r','LineWidth',2);           % 仿真误比特率
hold on
semilogy(EbN0dB,theoryBer,'--*b','LineWidth',2);        % 理论误比特率 BPSK
% axis([0 10 10^-5 1])
grid on
legend('simulation', 'theory');
xlabel('Eb/No, dB')
ylabel('Bit Error Rate')
title('Bit error probability curve for BPSK using OFDM')
