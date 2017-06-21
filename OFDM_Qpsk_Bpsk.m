%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  对思考题2作出仿真来验证结论的正确性                                                           %%
%  系统带宽3M，128点IFFT，CP长度为32，试分析当采用QPSK调制时，低通滤波器的参数是否需要改变？为什么？%%
%                                                                                             %%
%  收发端：FIR滤波器，等波纹滤波器（与之前的滤波器相同）                                          %%
%  信道：AWGN信道                                                                              %%
%  OFDM for QPSK：Subcarriers 分别为 64 和 128                                                 %%
%                                                                                             %%
%                                         Yuzhe Yang                                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

%---------------------------------------- 系统参数定义 ----------------------------------------%
nFFT = 128;             % n点fft的长度
nBpsk_Ofdm = 128;       % 每个OFDM符号带有的BPSK符号
nBit_Sym = 128;         % 每个OFDM符号所带的bit数目（对于BPSK而言和 nBpsk_Ofdm 一样）
nSym = 10^4;            % OFDM 符号数目（bit number > 10^5，即symbol > 10^3，此处取10^4保证数目精确度）

EbN0dB = 0:0.5:12;      % Bit to Noise Ratio
EsN0dB = EbN0dB + 10*log10(nBpsk_Ofdm/nFFT) - 10*log10(8);         % BPSK
EsN0dB2 = EbN0dB + 10*log10(2) - 10*log10(8);                      % QPSK

%---------------------------------------- 滤波器设计 及 其幅相频响应 ----------------------------------------%
[n,f,a,w] = firpmord([3.5*10^6 4.5*10^6], [1, 0], [0.005, 0.001], 24*10^6);
Hd = firpm(n,f,a,w);                        % 设计滤波器的系数


%---------------------------------------- OFDM 系统误比特率仿真 ----------------------------------------%
for ii = 1:length(EbN0dB)

    % 发射端
    s_bit = round(rand(1,nBit_Sym*nSym));               % 随机0, 1比特序列
    s_bit_bpsk = 2*s_bit - 1;                           % BPSK调制： 0 --> -1, 1 --> +1
    s_bit_bpsk = reshape(s_bit_bpsk,nBit_Sym,nSym).';   % 将调制后信号重新分组为 nSym*nBit_Sym 矩阵，便于操作
    
    I_bit = 2*(s_bit(1:2:end) - 0.5);
    Q_bit = 2*(s_bit(2:2:end) - 0.5);
    s_bit_qpsk = I_bit + 1j*Q_bit;                      % QPSK调制：1+j, 1-j, -1+j, -1-j
    s_bit_qpsk = reshape(s_bit_qpsk,nBit_Sym,nSym/2).';
    
    % IFFT变换 实现OFDM每个符号
    s_ofdm_bpsk = (nFFT/sqrt(nBpsk_Ofdm))*ifft(fftshift(s_bit_bpsk.')).';       % BPSK
    s_ofdm_qpsk = (nFFT/sqrt(nBpsk_Ofdm/2))*ifft(fftshift(s_bit_qpsk.')).';     % QPSK

    % 插入循环前缀CP，设定其长度为 L=32
    s_cp_bpsk = [s_ofdm_bpsk(:,97:128) s_ofdm_bpsk];            % BPSK
    s_cp_qpsk = [s_ofdm_qpsk(:,97:128) s_ofdm_qpsk];            % QPSK

    % 将多个OFDM符号串联到一个一维向量
    s_cp_bpsk = reshape(s_cp_bpsk.',1,nSym*(32+128));           % BPSK
    s_cp_qpsk = reshape(s_cp_qpsk.',1,nSym/2*(32+128));         % QPSK
    
    % 将实部与虚部信号分别处理
    s_I_bpsk = real(s_cp_bpsk);
    s_Q_bpsk = imag(s_cp_bpsk);
    s_I_qpsk = real(s_cp_qpsk);
    s_Q_qpsk = imag(s_cp_qpsk);
    
    % 8倍上采样
    s_Iup_bpsk = upsample(s_I_bpsk,8);
    s_Qup_bpsk = upsample(s_Q_bpsk,8);
    s_Iup_bpsk = reshape(s_Iup_bpsk, 8*(nBit_Sym+32), nSym).';      % reshape
    s_Qup_bpsk = reshape(s_Qup_bpsk, 8*(nBit_Sym+32), nSym).';
    
    s_Iup_qpsk = upsample(s_I_qpsk,8);
    s_Qup_qpsk = upsample(s_Q_qpsk,8);
    s_Iup_qpsk = reshape(s_Iup_qpsk, 8*(nBit_Sym+32), nSym/2).';      % reshape
    s_Qup_qpsk = reshape(s_Qup_qpsk, 8*(nBit_Sym+32), nSym/2).';
    
    % 经过发射端低通滤波器
    s_I_filter_bpsk = zeros(nSym,160*8 + n);
    s_Q_filter_bpsk = zeros(nSym,160*8 + n);
    for j = 1:nSym
        s_I_filter_bpsk(j,:) = conv(s_Iup_bpsk(j,:),Hd);                  % conv卷积
        s_Q_filter_bpsk(j,:) = conv(s_Qup_bpsk(j,:),Hd);
    end
    
    s_I_filter_qpsk = zeros(nSym/2,160*8 + n);
    s_Q_filter_qpsk = zeros(nSym/2,160*8 + n);
    for j = 1:nSym/2
        s_I_filter_qpsk(j,:) = conv(s_Iup_qpsk(j,:),Hd);                  % conv卷积
        s_Q_filter_qpsk(j,:) = conv(s_Qup_qpsk(j,:),Hd);
    end
    
    % 由于卷积引入了多余项，需要去掉（滤波器阶数 = 点数 = n）
    s_I_filter_bpsk = s_I_filter_bpsk(:,(n/2+1):end-n/2);
    s_Q_filter_bpsk = s_Q_filter_bpsk(:,(n/2+1):end-n/2);
    s_I_filter_bpsk = reshape(s_I_filter_bpsk.', 1, 8*nSym*(32+128));     % reshape
    s_Q_filter_bpsk = reshape(s_Q_filter_bpsk.', 1, 8*nSym*(32+128));
    
    s_I_filter_qpsk = s_I_filter_qpsk(:,(n/2+1):end-n/2);
    s_Q_filter_qpsk = s_Q_filter_qpsk(:,(n/2+1):end-n/2);
    s_I_filter_qpsk = reshape(s_I_filter_qpsk.', 1, 8*nSym/2*(32+128));     % reshape
    s_Q_filter_qpsk = reshape(s_Q_filter_qpsk.', 1, 8*nSym/2*(32+128));
    
    % 两路正交信号经过LPF后合并
    s_filter_bpsk = s_I_filter_bpsk + 1j*s_Q_filter_bpsk;
    s_filter_qpsk = s_I_filter_qpsk + 1j*s_Q_filter_qpsk;
    
    % 通过AWGN信道
    s_awgn_bpsk = awgn(s_filter_bpsk, EsN0dB(ii), 'measured');
    s_awgn_qpsk = awgn(s_filter_qpsk, EsN0dB2(ii), 'measured');
    
    % 接收端将实部与虚部信号分别处理
    r_I_bpsk = real(s_awgn_bpsk);
    r_Q_bpsk = imag(s_awgn_bpsk);
    r_I_bpsk = reshape(r_I_bpsk, 8*(nBit_Sym+32), nSym).';        % reshape
    r_Q_bpsk = reshape(r_Q_bpsk, 8*(nBit_Sym+32), nSym).';
    
    r_I_qpsk = real(s_awgn_qpsk);
    r_Q_qpsk = imag(s_awgn_qpsk);
    r_I_qpsk = reshape(r_I_qpsk, 8*(nBit_Sym+32), nSym/2).';        % reshape
    r_Q_qpsk = reshape(r_Q_qpsk, 8*(nBit_Sym+32), nSym/2).';
    
    % 经过接收端低通滤波器
    r_I_filter_bpsk = zeros(nSym,160*8 + n);
    r_Q_filter_bpsk = zeros(nSym,160*8 + n);
    for j = 1:nSym
        r_I_filter_bpsk(j,:) = conv(r_I_bpsk(j,:),Hd);            % conv卷积
        r_Q_filter_bpsk(j,:) = conv(r_Q_bpsk(j,:),Hd);
    end
    
    r_I_filter_qpsk = zeros(nSym/2,160*8 + n);
    r_Q_filter_qpsk = zeros(nSym/2,160*8 + n);
    for j = 1:nSym/2
        r_I_filter_qpsk(j,:) = conv(r_I_qpsk(j,:),Hd);            % conv卷积
        r_Q_filter_qpsk(j,:) = conv(r_Q_qpsk(j,:),Hd);
    end
    
    % 由于卷积引入了多余项，需要去掉（滤波器阶数 = 点数 = n）
    r_I_filter_bpsk = r_I_filter_bpsk(:,(n/2+1):end-n/2);
    r_Q_filter_bpsk = r_Q_filter_bpsk(:,(n/2+1):end-n/2);
    r_I_filter_bpsk = reshape(r_I_filter_bpsk.', 1, 8*nSym*(32+128));     % reshape
    r_Q_filter_bpsk = reshape(r_Q_filter_bpsk.', 1, 8*nSym*(32+128));
    
    r_I_filter_qpsk = r_I_filter_qpsk(:,(n/2+1):end-n/2);
    r_Q_filter_qpsk = r_Q_filter_qpsk(:,(n/2+1):end-n/2);
    r_I_filter_qpsk = reshape(r_I_filter_qpsk.', 1, 8*nSym/2*(32+128));     % reshape
    r_Q_filter_qpsk = reshape(r_Q_filter_qpsk.', 1, 8*nSym/2*(32+128));
    
    % 8倍下采样
    r_Idown_bpsk = downsample(r_I_filter_bpsk,8);
    r_Qdown_bpsk = downsample(r_Q_filter_bpsk,8);
    
    r_Idown_qpsk = downsample(r_I_filter_qpsk,8);
    r_Qdown_qpsk = downsample(r_Q_filter_qpsk,8);
    
    % 两路正交信号合并
    r_filter_bpsk = r_Idown_bpsk + 1j*r_Qdown_bpsk;
    r_cp_bpsk = reshape(r_filter_bpsk,160,nSym).';        % 将调制后信号重新分组为 nSym*(nBit_Sym+L) 矩阵，便于操作
    
    r_filter_qpsk = r_Idown_qpsk + 1j*r_Qdown_qpsk;
    r_cp_qpsk = reshape(r_filter_qpsk,160,nSym/2).';
    
    % 丢弃循环前缀CP，L=32
    r_fft_bpsk = r_cp_bpsk(:,33:end);
    r_fft_qpsk = r_cp_qpsk(:,33:end);
    
    % 转换到频域，利用FFT将OFDM符号转换为待解调的信息符号
    r_F_bpsk = (sqrt(nBpsk_Ofdm)/nFFT)*fftshift(fft(r_fft_bpsk.')).';
    r_F_qpsk = (sqrt(nBpsk_Ofdm/2)/nFFT)*fftshift(fft(r_fft_qpsk.')).';
    
    % BPSK 解调：将 +value --> 1, 将 -value --> -1
    r_demod_bpsk = real(r_F_bpsk);
    r_demod_bpsk(find(r_demod_bpsk > 0)) = +1;
    r_demod_bpsk(find(r_demod_bpsk < 0)) = -1;    
    r_bit_bpsk = (r_demod_bpsk + 1)/2;
    r_bit_bpsk = reshape(r_bit_bpsk.',nBit_Sym*nSym,1).';
    
    % QPSK 解调：实、虚部分别BPSK解调后重新组合
    r_demod_I = real(r_F_qpsk);
    r_demod_Q = imag(r_F_qpsk);
    r_demod_I(find(r_demod_I > 0)) = 1;
    r_demod_I(find(r_demod_I < 0)) = 0;
    r_demod_Q(find(r_demod_Q > 0)) = 1;
    r_demod_Q(find(r_demod_Q < 0)) = 0;
    r_demod_I = reshape(r_demod_I.', 1, nSym/2*nBit_Sym);
    r_demod_Q = reshape(r_demod_Q.', 1, nSym/2*nBit_Sym);
    r_bit_qpsk = reshape([r_demod_I; r_demod_Q],1,nBit_Sym*nSym);
    
    % 计算错误bits数目
    nErr_bpsk(ii) = size(find(r_bit_bpsk - s_bit),2);
    nErr_qpsk(ii) = size(find(r_bit_qpsk - s_bit),2);
end

simBer_bpsk = nErr_bpsk/(nSym*nBit_Sym);                % 仿真误比特率 BPSK
simBer_qpsk = nErr_qpsk/(nSym*nBit_Sym);                % 仿真误比特率 QPSK
theoryBer = (1/2)*erfc(sqrt(10.^(EbN0dB/10)));          % 理论误比特率 BPSK/QPSK
% load pro.mat simBer_bpsk theoryBer
figure
semilogy(EbN0dB,simBer_bpsk,'--*r','LineWidth',2);      % 仿真误比特率 BPSK
hold on
semilogy(EbN0dB,simBer_qpsk,'--*g','LineWidth',2);      % 仿真误比特率 QPSK
hold on
semilogy(EbN0dB,theoryBer,'--*b','LineWidth',2);        % 理论误比特率 BPSK/QPSK

grid on
legend('BPSK', 'QPSK', 'theory');
xlabel('Eb/No, dB')
ylabel('Bit Error Rate')
title('Bit error probability curve for BPSK and QPSK using OFDM')
