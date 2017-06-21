%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  ��OFDMϵͳ���з��棬�õ������������                                                          %%
%  ����BPSK���е��ƣ�ϵͳ��128��BPSK���Ʒ��Ž���IFFT����                                         %%
%  IFFT����������һ��OFDM����ϵͳ��128��BPSK���Ʒ��Ž���IFFT���㣬IFFT����������һ��OFDM����   %%
%  �շ��ˣ�FIR�˲������Ȳ����˲���                                                              %%
%  �ŵ���AWGN�ŵ�                                                                              %%
%                                                                                             %%
%                                             ���                                           %%
%                                           1400012996                                        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

%---------------------------------------- ϵͳ�������� ----------------------------------------%
nFFT = 128;             % n��fft�ĳ���
nBpsk_Ofdm = 128;       % ÿ��OFDM���Ŵ��е�BPSK����
nBit_Sym = 128;         % ÿ��OFDM����������bit��Ŀ������BPSK���Ժ� nBpsk_Ofdm һ����
nSym = 10^4;            % OFDM ������Ŀ��bit number > 10^5����symbol > 10^3���˴�ȡ10^4��֤��Ŀ��ȷ�ȣ�

EbN0dB = 0:0.5:12;      % Bit to Noise Ratio
EsN0dB = EbN0dB + 10*log10(nBpsk_Ofdm/nFFT) - 10*log10(8);         % ת��Ϊÿ��symbol��SNR, Nsample = 8


%---------------------------------------- �˲������ �� �����Ƶ��Ӧ ----------------------------------------%
[n,f,a,w] = firpmord([3.5*10^6 4.5*10^6], [1, 0], [0.005, 0.001], 24*10^6);
Hd = firpm(n,f,a,w);                        % ����˲�����ϵ��
disp(Hd)                                    % FIRϵͳ����������ʾ

freqz(Hd, 1, 'whole');                      % �˲�����Ƶ\��Ƶ��Ӧ
figure
zplane(Hd, 1); axis([-1.1,1.1,-1.1,1.1]); title('�㼫��ͼ');        % �˲������㼫��ͼ
grpdelay(Hd);  title('Ⱥ�ӳ�ͼ');                                   % �˲�����Ⱥ�ӳ�
stem(Hd); title('�˲���ʱ����Ӧ/�弤��Ӧ');                          % �˲���ʱ����Ӧ

%---------------------------------------- OFDM ϵͳ������ʷ��� ----------------------------------------%
for ii = 1:length(EbN0dB)

    % �����
    s_bit = round(rand(1,nBit_Sym*nSym));               % ���0, 1��������
    s_bit_bpsk = 2*s_bit - 1;                           % BPSK���ƣ� 0 --> -1, 1 --> +1
    s_bit_bpsk = reshape(s_bit_bpsk,nBit_Sym,nSym).';   % �����ƺ��ź����·���Ϊ nSym*nBit_Sym ���󣬱��ڲ���
    % disp(s_bit_bpsk(1:10));
  
    % IFFT�任 ʵ��OFDMÿ������
    % (nFFT/sqrt(nBpsk_Ofdm))�ǽ�ÿ�����ͷ��ŵ�������һ����1
    s_ofdm = (nFFT/sqrt(nBpsk_Ofdm))*ifft(fftshift(s_bit_bpsk.')).';    % ifftǰ��Ҫfftshift�����ģ�������Ϊ�����shift����Ҫ��ת��

    % ����ѭ��ǰ׺CP���趨�䳤��Ϊ L=16
    s_cp = [s_ofdm(:,113:128) s_ofdm];

    % �����OFDM���Ŵ�����һ��һά����
    s_cp = reshape(s_cp.',1,nSym*(16+128));
    
    % ��ʵ�����鲿�źŷֱ���
    s_I = real(s_cp);
    s_Q = imag(s_cp);
    
    % 8���ϲ���
    s_Iup = upsample(s_I,8);
    s_Qup = upsample(s_Q,8);
    
    s_Iup = reshape(s_Iup, 8*(nBit_Sym+16), nSym).';      % reshape
    s_Qup = reshape(s_Qup, 8*(nBit_Sym+16), nSym).';
    
    %{
    % ��������˵�ͨ�˲���
    s_I_filter = zeros(nSym,144*8);
    s_Q_filter = zeros(nSym,144*8);
    for j = 1:nSym
        s_I_filter(j,:) = filter(Hd,1,s_Iup(j,:));      % ���, conv()���������⣿����ʱ��filter
        s_Q_filter(j,:) = filter(Hd,1,s_Qup(j,:));
    end
    
    s_I_filter = reshape(s_I_filter.', 1, 8*nSym*(16+128));     % reshape
    s_Q_filter = reshape(s_Q_filter.', 1, 8*nSym*(16+128));
    %}

    % ��������˵�ͨ�˲���
    % Hd = my_filter;
    s_I_filter = zeros(nSym,144*8 + n);
    s_Q_filter = zeros(nSym,144*8 + n);
    for j = 1:nSym
        s_I_filter(j,:) = conv(s_Iup(j,:),Hd);                  % conv���
        s_Q_filter(j,:) = conv(s_Qup(j,:),Hd);
    end
    
    % ���ھ�������˶������Ҫȥ�����˲������� = ���� = n��
    s_I_filter = s_I_filter(:,(n/2+1):end-n/2);
    s_Q_filter = s_Q_filter(:,(n/2+1):end-n/2);
    
    s_I_filter = reshape(s_I_filter.', 1, 8*nSym*(16+128));     % reshape
    s_Q_filter = reshape(s_Q_filter.', 1, 8*nSym*(16+128));
    
    % ��·�����źž���LPF��ϲ�
    s_filter = s_I_filter + 1j*s_Q_filter;
    
    % ͨ��AWGN�ŵ�
    s_awgn = awgn(s_filter, EsN0dB(ii), 'measured');
    
    % ���ն˽�ʵ�����鲿�źŷֱ���
    r_I = real(s_awgn);
    r_Q = imag(s_awgn);
    
    r_I = reshape(r_I, 8*(nBit_Sym+16), nSym).';        % reshape
    r_Q = reshape(r_Q, 8*(nBit_Sym+16), nSym).';
    %{
    % �������ն˵�ͨ�˲���
    r_I_filter = zeros(nSym,144*8);
    r_Q_filter = zeros(nSym,144*8);
    for j = 1:nSym
        r_I_filter(j,:) = filter(Hd,1,r_I(j,:));      % ���
        r_Q_filter(j,:) = filter(Hd,1,r_Q(j,:));
    end
    
    r_I_filter = reshape(r_I_filter.', 1, 8*nSym*(16+128));     % reshape
    r_Q_filter = reshape(r_Q_filter.', 1, 8*nSym*(16+128));
    %}
    
    % �������ն˵�ͨ�˲���
    r_I_filter = zeros(nSym,144*8 + n);
    r_Q_filter = zeros(nSym,144*8 + n);
    for j = 1:nSym
        r_I_filter(j,:) = conv(r_I(j,:),Hd);            % conv���
        r_Q_filter(j,:) = conv(r_Q(j,:),Hd);
    end
    
    % ���ھ�������˶������Ҫȥ�����˲������� = ���� = n��
    r_I_filter = r_I_filter(:,(n/2+1):end-n/2);
    r_Q_filter = r_Q_filter(:,(n/2+1):end-n/2);
    
    r_I_filter = reshape(r_I_filter.', 1, 8*nSym*(16+128));     % reshape
    r_Q_filter = reshape(r_Q_filter.', 1, 8*nSym*(16+128));
    
    % 8���²���
    r_Idown = downsample(r_I_filter,8);
    r_Qdown = downsample(r_Q_filter,8);
    
    % ��·�����źźϲ�
    r_filter = r_Idown + 1j*r_Qdown;
    
    r_cp = reshape(r_filter,144,nSym).';        % �����ƺ��ź����·���Ϊ nSym*(nBit_Sym+L) ���󣬱��ڲ���
    
    % ����ѭ��ǰ׺CP��L=16
    r_fft = r_cp(:,17:end);
    
    % ת����Ƶ������FFT��OFDM����ת��Ϊ���������Ϣ����
    r_F = (sqrt(nBpsk_Ofdm)/nFFT)*fftshift(fft(r_fft.')).';

    % BPSK ���
    % �� +value --> 1, �� -value --> -1
    r_demod = real(r_F);
    r_demod(find(r_demod > 0)) = +1;
    r_demod(find(r_demod < 0)) = -1;

    % ����������Ϣ����ת��Ϊ������ʽ
    r_bit = (r_demod + 1)/2;
    r_bit = reshape(r_bit.',nBit_Sym*nSym,1).';

    % �������bits��Ŀ
    nErr(ii) = size(find(r_bit - s_bit),2);
        
end

simBer = nErr/(nSym*nBit_Sym);                          % �����������
theoryBer = (1/2)*erfc(sqrt(10.^(EbN0dB/10)));          % ����������� BPSK

figure
semilogy(EbN0dB,simBer,'--*r','LineWidth',2);           % �����������
hold on
semilogy(EbN0dB,theoryBer,'--*b','LineWidth',2);        % ����������� BPSK
% axis([0 10 10^-5 1])
grid on
legend('simulation', 'theory');
xlabel('Eb/No, dB')
ylabel('Bit Error Rate')
title('Bit error probability curve for BPSK using OFDM')
