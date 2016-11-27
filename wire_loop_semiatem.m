%--------------------------------------------------------------------------
%��������������TEM �뺽��TEM
% ����ˮƽ��żԴ�������������й۲�

%--------------------------------------------------------------------------
%%
format long;
clear all;clc;close all;
%%
u0 = 4*pi*1e-7;
load parameters.txt;
sigma1 = parameters(1,2);%��һ��ĵ絼��
rou = 1./sigma1;
%% �������������
x = 0; % �շ�ˮƽƫ�ƾ࣬��y��
y = 1000;
L= 2500; % �������³��ȣ���x��
I = 40; % �������
m = I*L;
%%  �뺽���շ��߶Ȳ���
z =-100;% �۲������ĸ߶ȣ���������Ϊ��ֵ
h =0;% Դ�����ĸ߶�
n=12;
%% �����ʺ͹۲�ʱ�������
fs = 1e7;% ������
dt = 1./fs;
% t = 1/fs:1/fs:2e-3;% ʱ������
t = logspace(-5,0,100);
%%
[step_ex01,step_ey01,impulse_hx,impulse_hy,impulse_hz] = Horizontal_Electrical_Dipole_Response(u0,rou,m,x,y,z,t);
% [hz_01,hz_10,hz_impulse,hx_01,hx_10,hx_impulse,hy_01,hy_10,hy_impulse,ex_01,ex_impulse,ey_01,ey_impulse]...
% = Calculate_Horizontal_Electrical_Dipole(I,L,h,x,y,z,t);
% tic;
% fine = 0.4;
% [hz_01,hz_10,hz_impulse,hx_01,hx_10,hx_impulse,hy_01,hy_10,hy_impulse,ex_01,ex_10,ex_impulse,ey_01,ey_10,ey_impulse] = ...
%     Calculate_Horizontal_Finite_Electrical_Source(I,L,h,x,y,z,t,fine);% ���ø�˹-���õ»��ֵ�ԭ����ⳤ�ӵص���Դ��˲���ų�
% t0 = toc
tic;
[hz_01,hz_10,hz_impulse,hx_01,hx_10,hx_impulse,hy_01,hy_10,hy_impulse,ex_01,ex_10,ex_impulse,ey_01,ey_10,ey_impulse] = ...
    Calculate_Horizontal_Finite_Electrical_Source_GuaLeg(I,L,h,x,y,z,t,n);% ���ø�˹-���õ»��ֵ�ԭ����ⳤ�ӵص���Դ��˲���ų�
t1 = toc
save('horizontal_electrical_dipole_impulse_shuzhijie','hx_impulse','hy_impulse','hz_impulse');
%% z���Ծ
%{
%    ��ͼ��
hz_10 = hz_01(end) - hz_01;
figure;
plot(t.*10^3,u0.*hz_01,'r','Linewidth',1);
hold on
plot(t.*10^3,u0.*hz_10,'b','Linewidth',1);
grid on;
legend('��ֵ������Ծ��Ӧ','��ֵ�⸺��Ծ��Ӧ');
title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop Bz'])
xlabel('Time/(ms)')
ylabel('Bz/(A/m)');
% x���Ծ
%��ͼ��
% hx_10 = hx_01(end) - hx_01;
figure;
plot(t.*10^3,u0.*(hx_01),'r','Linewidth',1);
hold on
plot(t.*10^3,-u0.*(hx_01),'b','Linewidth',1);
grid on;
legend('��ֵ������Ծ��Ӧ','��ֵ�⸺��Ծ��Ӧ');
title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop Bx'])
xlabel('Time/(ms)')
ylabel('Bx/(A/m)');
% y ���Ծ
hy_10 = hy_01(end) - hy_01;
figure;
plot(t.*10^3,u0.*(hy_01),'r','Linewidth',1);
hold on
plot(t.*10^3,u0.*(hy_10),'b','Linewidth',1);
grid on;
legend('��ֵ������Ծ��Ӧ','��ֵ�⸺��Ծ��Ӧ');
title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop By'])
xlabel('Time/(ms)')
ylabel('By/(A/m)');
%}
%% z ������
jiexijie = ['Horizontal_Electrical_Dipole_Response_rou' num2str(rou) '_xr' num2str(x) '_yr' num2str(y) '_zr' num2str(z) '.mat'];
load(jiexijie);
figure;
loglog(t(1:end).*10^3,abs((ex_10))./L,'r','Linewidth',2);
hold on
loglog(t(1:end).*10^3,abs(step_ex01(end)-step_ex01)./L,'k:','Linewidth',2);
grid on;
legend('��ֵ��ex\_01','������step\_ex01');
title(['source moment' num2str(I) 'A*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop ex step response'])
xlabel('Time/(ms)')
ylabel('Ex/(V/m)');
figure;
loglog(t(1:end).*10^3,abs((ey_10))./L,'r','Linewidth',2);
hold on
loglog(t(1:end).*10^3,abs(step_ey01(end)-step_ey01)./L,'k:','Linewidth',2);
grid on;
legend('��ֵ��ey\_01','������step\_ey01');
title(['source moment' num2str(I) 'A*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop ey step response'])
xlabel('Time/(ms)')
ylabel('Ey/(V/m)');
% ������
% ex_error = abs(abs(ex_01)-abs(step_ex01))./abs(step_ex01);
% ey_error = abs(abs(ey_01)-abs(step_ey01))./abs(step_ey01);
% figure;
% loglog(t.*1e3,ex_error.*100,'r','linewidth',2);
% hold on;
% loglog(t.*1e3,ey_error.*100,'k:','linewidth',2);
% grid on;
% legend('ex','ey');
% title('��ֵ��ͽ������������');
% xlabel('Time/(ms)')
% ylabel('error/(%)');
% ������Ӧ
% figure;
% plot(t(1:end).*10^3,ex_impulse,'r','Linewidth',2);
% hold on
% plot(t(1:end-1).*10^3,diff(step_ex01)./dt,'k:','Linewidth',2);
% grid on;
% legend('��ֵ��ex\_01','������step\_ex01');
% title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop ex step response'])
% xlabel('Time/(ms)')
% ylabel('Ex/(V/m)');
% figure;
% plot(t(1:end).*10^3,abs(ey_01),'r','Linewidth',2);
% hold on
% plot(t(1:end).*10^3,abs(step_ey01),'k:','Linewidth',2);
% grid on;
% legend('��ֵ��ey\_01','������step\_ey01');
% title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop ey step response'])
% xlabel('Time/(ms)')
% ylabel('Ey/(V/m)');
%
%%
figure;
loglog(t(1:end).*10^3,u0.*abs(hz_impulse)./L,'r','Linewidth',2);
hold on
loglog(t(1:end-1).*10^3,u0.*abs(diff(hz_01)./L./diff(t)),'b','Linewidth',2);
grid on;
loglog(t(1:end).*10^3,u0.*abs(impulse_hz)./L,'k:','Linewidth',2);
grid on;
legend('��ֵ��hz\_1\_impulse','��ֵ��diff(hz\_01)./dt','������impulse\_hz');
title(['source moment' num2str(I) 'A*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop Bz impulse response'])
xlabel('Time/(ms)')
ylabel('Bz/(T)');

figure;
loglog(t(1:end).*10^3,u0.*abs(hx_impulse)./L,'r','Linewidth',2);
hold on
loglog(t(1:end-1).*10^3,u0.*abs(diff(hx_01)./L./diff(t)),'b','Linewidth',2);
grid on;
loglog(t(1:end).*10^3,u0.*abs(impulse_hx)./L,'k:','Linewidth',2);
grid on;
legend('��ֵ��hx\_1\_impulse','��ֵ��diff(hx\_01)./dt','������impulse\_hx');
title(['source moment' num2str(I) 'A*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop Bx impulse response'])
xlabel('Time/(ms)')
ylabel('Bx/(T)');
figure;
loglog(t(1:end).*10^3,u0.*abs(hy_impulse)./L,'r','Linewidth',2);
hold on
loglog(t(1:end-1).*10^3,u0.*abs(diff(hy_01)./L./diff(t)),'b','Linewidth',2);
grid on;
loglog(t(1:end).*10^3,u0.*abs(impulse_hy)./L,'k:','Linewidth',2);
grid on;
legend('��ֵ��hy\_1\_impulse','��ֵ��diff(hy\_01)./dt','������impulse\_hy');
title(['source moment' num2str(I) 'A*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop By impulse response'])
xlabel('Time/(ms)')
ylabel('By/(T)');
% ������
hz_error = abs(abs(hz_impulse)-abs(impulse_hz))./abs(impulse_hz);
hx_error = abs(abs(hx_impulse)-abs(impulse_hx))./abs(impulse_hx);
hy_error = abs(abs(hy_impulse)-abs(impulse_hy))./abs(impulse_hy);
figure;
loglog(t.*1e3,hz_error.*100,'r','linewidth',2);
hold on;
loglog(t.*1e3,hx_error.*100,'b','linewidth',2);
hold on;
loglog(t.*1e3,hy_error.*100,'k:','linewidth',2);
grid on;
legend('hz','hx','hy');
title('��ֵ��ͽ������������');
xlabel('Time/(ms)')
ylabel('error/(%)');

%% save data

