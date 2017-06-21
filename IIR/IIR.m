%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  IIR�˲�����ƣ�IIR�����˲�����Դ�ź�����ȡ�ض�Ƶ�ʵ��ź�                      %%
%  �ֶ�¼�ƣ��� 8kHz ����Ƶ��¼���Լ������������Ӹ�˹�������������С��20dB       %%
%                                                                            %%
%  ������ֵ�ͨ�˲�����ƣ�������Ӧ���䷨ + �б�ѩ��I�ͣ� ˫���Ա任�� + ������˹ %%
%                                                                            %%
%                                     ���                                  %%
%                                   1400012996                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%------------------------------------- ¼������ ----------------------------------%
% 
% ����ƽ̨��Windows 10 64bit��MATLAB R2014a
% ¼��5���ӣ�����Ƶ��Ĭ��Ϊ 8kHz
%{
recObj = audiorecorder;
disp('Start speaking.')
recordblocking(recObj, 5);
disp('End of Recording.');
play(recObj);       % �ط�¼������
myRecording = getaudiodata(recObj);         % ��ȡ¼������

plot(myRecording);          % ����¼������ʱ����
plot(abs(fftshift(fft(myRecording))));      % ����¼������Ƶ��
%}


%--------------------------------- ���ݴ����� ----------------------------------%
load myrecord.mat myRecording           % ������¼�ƺõ���Ƶ

%------------------------------- ��FFT������ԭ�źŷ�����
N1 = 40000;     %  ����������ʱ�䳤�� 5s
fs = 8000;      %  ����Ƶ�� 8kHz
n1 = 0:N1-1;
t1 = n1/fs;     %  ʱ������

myrecord_F = fftshift(fft(myRecording,N1));    %  ���źŽ���FFT,��ת�Ƶ�����
mag1 = abs(myrecord_F);             %  ������
F1 = n1*fs/N1 - fs/2;       %  Ƶ�����г��ȣ��ҽ���Ƶ�Ƶ�����

figure
plot(t1,myRecording);          %  ����¼������ʱ����
title('Դ�ź�ʱ����');
figure
plot(F1, mag1, '-r');          %  ����¼������Ƶ�ף��������� 2kHz�ڣ�������Ƶ�׷�Χ85~1100Hz��
axis([-2000 2000 0 1500]), xlabel('f/Hz'), title('Դ�ź�Ƶ��');

%------------------------------- ��ԭ�źżӸ�˹�������������5dB��<20dB��
SNR = 5;        %  ����� 5dB
myrecord_noise = awgn(myRecording,SNR,'measured');
mag2 = abs(fftshift(fft(myrecord_noise,N1)));

figure
plot(t1,myrecord_noise);          %  ���Ƽ�������ʱ����
title('���������ź�ʱ����');
figure
plot(F1, mag2, '-r');             %  ���Ƽ�����������Ƶ�ף��������� 2kHz��
axis([-2000 2000 0 1500]), xlabel('f/Hz'), title('���������ź�Ƶ��');


%------------------------------- ����1��������Ӧ���䷨ ��� ��ͨ�б�ѩ��I�� -------------------------------%
%
%   �˲���������ͨ������300Hz, ����ٽ�500Hz������Ƶ��8kHz
%              ͨ�����˥��1dB�������С˥��50dB
wp = 2*pi*300;
ws = 2*pi*500;      %  �����˲�������ת����ģ���˲���Ƶ������
Rp = 1;
Rs = 50;
Fs = 8000;
[N, Wn] = cheb1ord(wp,ws,Rp,Rs,'s');        %  ѡ���˲�������С����
[Z, P, K] = cheb1ap(N,Rp);                  %  ������ͨ�б�ѩ��I�˲���
[A, B, C, D] = zp2ss(Z,P,K);                %  ��һ����ͨ�˲���
[At, Bt, Ct, Dt] = lp2lp(A,B,C,D,Wn);       %  ȥ��һ�����õ���ֹƵ�� Wn ��
[num_ana,den_ana] = ss2tf(At,Bt,Ct,Dt);     %  �õ�ģ���˲������ݺ��������ӡ���ĸ����ʽϵ����

[num_dig,den_dig] = impinvar(num_ana,den_ana,Fs);      %  ������Ӧ���䷨ ת��Ϊ�����˲������ݺ��������ӡ���ĸ����ʽϵ����
disp(num_dig)
disp(den_dig)                %  IIRϵͳ�����㼫��
figure
[H, W] = freqz(num_dig,den_dig,N1,'whole');            %  �˲�����Ƶ��Ӧ
plot(W*Fs/2/pi,20*log10(abs(H)));
axis([0 600 -80 0])
grid; xlabel('f/Hz'); title('����1��������Ӧ���䷨���б�ѩ��I�͵�ͨ�˲���');
figure
subplot(121)
plot(W,angle(H)); xlabel('\omega'); ylabel('\phi'); title('��λ��Ӧ'); axis([0,pi,-pi,pi]);      %  ��Ƶ��Ӧ��ϵͳ�㼫ͼ
subplot(122)
zplane(num_dig,den_dig); axis([-1.1,1.1,-1.1,1.1]); title('�㼶ͼ');

%---------------------------------- �����źž������˲���
out_mag = abs(mag2.*fftshift(H));
figure
plot(F1, out_mag, '-r');             %  �����˲�������Ƶ�ף��������� 2kHz��
axis([-2000 2000 0 1500]), xlabel('f/Hz'), title('�˲����ź�Ƶ��');


%------------------------------- ����2��˫���Ա任�� ��� ��ͨ������˹ --------------------------------%
%
%   �˲���������ͨ������300Hz, ����ٽ�500Hz������Ƶ��8kHz
%              ͨ�����˥��1dB�������С˥��50dB
Wp = 0.075*pi;            %  ͨ���߽�Ƶ��(��һ��)��wp=fp*2*pi/fs
Ws = 0.125*pi;            %  ����߽�Ƶ��(��һ��)��ws=fr*2*pi/fs
Rp = 1;                   %  ͨ������
Rs = 50;                  %  ���˥��
Ts = 0.000125;            %  Fs = 8kHz
%   ת��Ϊģ���˲���ָ�꣨Ԥ���䣩
OmegaP = (2/Ts)*tan(Wp/2);        %  ģ���ͨԭ���˲���ͨ��Ƶ��
OmegaS = (2/Ts)*tan(Ws/2);        %  ģ���ͨԭ���˲������Ƶ��

[N, OmegaC] = buttord(OmegaP,OmegaS,Rp,Rs,'s');     %  ģ�������˹�˲����Ľ�����-3dB��ֹƵ�ʼ���

[Z, P, K] = buttap(N);                              %  ��ƹ�һ�����������˲�����ͨԭ��
num_ana = K * real(poly(Z));
den_ana = real(poly(P));                            %  ԭ�ͣ���һ�����˲���ϵ��
[num_ana, den_ana] = lp2lp(num_ana,den_ana,OmegaC); %  ȥ��һ�������õ���ֹƵ�� OmegaC ��

[num_dig, den_dig] = bilinear(num_ana,den_ana,Fs);  %  ˫���Ա任�� ת��Ϊ�����˲������ݺ���
disp(num_dig)
disp(den_dig)                %  IIRϵͳ�����㼫��
figure
[H, W] = freqz(num_dig,den_dig,N1,'whole');         %  �˲�����Ƶ��Ӧ
plot(W*Fs/2/pi,20*log10(abs(H)));
axis([0 600 -80 0])
grid; xlabel('f/Hz'); title('����2��˫���Ա任����������˹��ͨ�˲���');

%---------------------------------- �����źž������˲���
out_mag = abs(mag2.*fftshift(H));
figure
plot(F1, out_mag, '-r');             %  �����˲�������Ƶ�ף��������� 2kHz��
axis([-2000 2000 0 1500]), xlabel('f/Hz'), title('�˲����ź�Ƶ��');
figure
subplot(121)
plot(W,angle(H)); xlabel('\omega'); ylabel('\phi'); title('��λ��Ӧ'); axis([0,pi,-pi,pi]);      %  ��Ƶ��Ӧ��ϵͳ�㼫ͼ
subplot(122)
zplane(num_dig,den_dig); axis([-1.1,1.1,-1.1,1.1]); title('�㼶ͼ');
