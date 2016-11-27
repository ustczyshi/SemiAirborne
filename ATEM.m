%% ���÷���õ��ķ����������ģ���伤���ľ��ȴ��˲���ų���
% �����о��ضϾֲ��쳣�Ը�Ӧ��ѹ�źŵ�Ӱ��
% by shizongyang
% 2015.11.24
%% 
clc;
close all;
clear all;
%% load the current data
% load Variable_step_discrete.mat;
load ATEM_current.mat;
fs = 1e7;
% fs = Fs;
current = atem_current(1:0.02*fs)';
t = (1:length(current))/fs;

%% the transmitting current
figure(1);
plot(t,current,'r');
title('current');
xlabel('time/s');
ylabel('I/A');

%%  ���ô�ֱ��ż�����ھ��Ȱ�ռ��ر�������Ĵ�ֱ�ų��ļ���ʽ����
% �����Ծ��Ӧ
% �����趨
u0 = 4*pi*1e-7;%H/m
rou =1;% ������
I = 1;
r = 10;
St = pi*r^2;% ������Ȧ�������m^2
Nt = 5;% ������Ȧ����
Qt = St*Nt;% ������Ȧ����Ч���
Sr = 1;% ������Ȧ�������m^2
Nr = 1;% ������Ȧ����
% Srq = Sr*Nr;
Srq = 50;
%% �趨������Ȧ������Ϊ����ԭ�㣬��ֱ����Ϊ+z
xr = 100;yr =0;zr = 0;% ������Ȧ����λ������m��
%%��ֱ��ż��Դ�µľ��ȴ�ؽ�Ծ��ӦHz Eo---
m = I*St;% ������Ȧ�Ĵž�
t_ob = 1/fs:1/fs:0.1;%2��Ĺ۲���
% [Hresponse01,Hresponse10,Himpulse10,Eresponse] = Vertical_Magnetic_Dipole_Step_Response(u0,rou,m,xr,yr,zr,t_ob);
[Hresponse01,Hresponse10,Himpulse10] = Central_Circle_Loop_Step_Response(u0,rou,I,r,t_ob);
%% ��������ͽ�Ծ��Ӧ
% himpulse10 = Himpulse10;
% hstep10 = Hresponse10;
% save('surface_magnetic_dipole_response_jx.mat','himpulse10','hstep10');

%% �����Ծ��Ӧ��������Ӧ
% figure;
% plot(t,Hresponse10(1:length(t)),'r');%����Ծ��Ӧ
% xlabel('time/s');
% ylabel('Hz');
% title('����Ծ��Ӧ');
% figure;
% plot(t(1:end),Himpulse10(1:length(t)),'b');
% hold on;
% title('��������Ӧ');
% xlabel('time/s');
% ylabel('Hz');
%%
% figure;
% plot(t_ob*1000,Hresponse10,'r');%����Ծ��Ӧ
% hold on;
% title('����Ծ��Ӧ');
% xlabel('time/ms');
% ylabel('Hz');
% figure;
% plot(t_ob(1:end-1)*1000,(diff(Hresponse10)*fs),'k');%����Ծ��Ӧ
% hold on;
% legend('����Ծ��Ӧ�ĵ���');
% xlabel('time/ms');
% ylabel('Hz');
% figure;
% plot(t_ob*1000,(Himpulse10),'b');
% hold on;
% legend('��������Ӧ');
% xlabel('time/ms');
% ylabel('Hz');
%% �����Ծ��Ӧ��������Ӧ ������
% figure;
% loglog(t_ob*1000,Hresponse10,'r');%����Ծ��Ӧ
% hold on;
% loglog(t_ob(1:end-1)*1000,abs(diff(Hresponse10)*fs),'k');%����Ծ��Ӧ
% hold on;
% % loglog(t_ob(1:end-1)*1000,abs(diff([Hresponse01])*fs),'m');%����Ծ��Ӧ
% % hold on;
% loglog(t_ob*1000,abs(Himpulse10),'b');
% hold on;
% legend('����Ծ��Ӧ','����Ծ��Ӧ�ĵ���','��������Ӧ');
% xlabel('time/ms');
% ylabel('Hz');
%
%% ������ĵ���
dIdt = Nt*diff([current 0]);%length(t)
response = conv(-dIdt,u0*Hresponse10);% 3*length(t)-1
Response = response(1:length(t_ob)); % ���γ�
%% ֱ������������Ӧ�ı��ʽ���������⣬��ֵƫ��--------������
% response1 = conv(current,-u0*Himpulse10)/fs;% 3*length(t)-1
% Response1 = response1(1:length(t_ob));
%%
response1 = conv(dIdt,u0*Hresponse01);% 3*length(t)-1
Response1 = response1(1:length(t_ob)); %���ܳ�
%%
figure;
plot((1:length(Response))/fs,Response,'r');% ����Ծ
hold on;
plot((1:length(Response1))/fs,Response1,'b');
legend('conv(dIdt,Bz\_step)','conv(I,Bz\_impulse)');
%%
figure;
[hx,h1,h2] = plotyy((1:length(Response))/fs,Response,(1:length(Response1))/fs,Response1);
grid on;
title(['��=' num2str(rou) '����m']);
legend('conv(di,Bstep-):Bz2','conv(di, hstep+):Bz');
xlabel('time/s');
ylabel(hx(1),'Bz2 /T','Color','b');% left y-axis
ylabel(hx(2),'Bz/T','Color','r');% right y-axis
h1.Color = 'b';
h2.Color = 'r';
%
%%
% Vinducted = -Srq*diff([Response Response(end) ])*fs;%3*length(t)-1
%%
Vinducted = -Srq*diff([Response1 Response1(end) ])*fs;%3*length(t)-1
%%
figure;
plot(t,current,'r');
hold on;
title('I');
xlabel('time/s');
ylabel('I/A');
%%
figure;
plot(t,Response1(1:length(t)),'r');
hold on;
plot(t,Nt*current.*(u0*Hresponse01(end)),'b');
hold on;
plot(t,Response1(1:length(t))-u0*Nt*current.*(Hresponse01(end)),'k');
hold on;
grid on;
legend('�ܳ�Bz','һ�γ�Bz1','���γ�Bz2');
xlabel('time/s');
ylabel('Bz/T');
%%
figure;
plot(t,Vinducted(1:length(t)),'r');
grid on;
title('Uz');
xlabel('time/s');
ylabel('U/V');
%
%% 
figure;
plot(t,(Vinducted(1:length(t))),'r');
grid on;
title('the induced voltage ');
xlabel('time/s');
ylabel('U/V');
%%
% figure;
% loglog(t,abs(Vinducted(1:length(t))),'b');
% legend('the induced voltage ');
% xlabel('time/s');
% ylabel('U/V');
%%
figure;
plot(t,(current(1:length(t)))./max(current(1:length(t))),'b');
hold on;
plot(t,Vinducted(1:length(t))./max(abs(Vinducted)),'r');
hold on;
grid on;
title(['��=' num2str(rou) '����m']);
legend('the normalized current','the normalized-induced voltage ');
xlabel('time/s');
ylabel('Uz/max(|Uz|)--I/max(|I|)');
%% ˫������
figure;
[hx,h1,h2] = plotyy(t(1:0.01*fs),current(1:0.01*fs),t(1:0.01*fs),1000*(Vinducted(1:0.01*fs)));
grid on;
title(['��=' num2str(rou) '����m']);
legend('current','the induced voltage ');
xlabel('time/s');
ylabel(hx(1),'the current of  loop /A','Color','b');% left y-axis
ylabel(hx(2),'Uz /mV','Color','r');% right y-axis
h1.Color = 'b';
h2.Color = 'r';
%%
%{
% current1 = current(1:0.01*fs);
% inducted_voltage1 = Vinducted(1:0.01*fs);
% %%
% save('current1','current1');
% save('inducted_voltage1','inducted_voltage1');
%%
% figure;
% loglog(t(1:0.01*fs),abs(Vinducted(1:0.01*fs)),'b');
% % hold on;
% % semilogy(t(1:100000),abs(Vinducted(0.01*fs:0.01*fs+100000-1)),'r');
% legend('the induced voltage ');
% xlabel('time/s');
%%
% figure;
% loglog(t(1e-4*fs:0.01*fs),1000*abs(Vinducted(1e-4*fs:0.01*fs)),'b');%%��Ϊ100us�ض�ʱ��
% legend('the induced voltage ');
% xlabel('time/s');
% ylabel('U/mV');
%}


