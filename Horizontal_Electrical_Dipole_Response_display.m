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
y = 0;
L= 1; % �������³��ȣ���x��
I = 1; % �������
m = I*L;
%%  �뺽���շ��߶Ȳ���
z =-100;% �۲������ĸ߶ȣ���������Ϊ��ֵ
h =0;% Դ�����ĸ߶�

%% �����ʺ͹۲�ʱ�������
fs = 1e7;% ������
dt = 1./fs;
% t = 1/fs:1/fs:2e-3;% ʱ������
t = logspace(-7,2,200);
% [hz_01,hz_10,hz_impulse,hx_01,hx_10,hx_impulse,hy_01,hy_10,hy_impulse,ex_01,ex_10,ex_impulse,ey_01,ey_10,ey_impulse] =...
%     Calculate_Horizontal_Finite_Electrical_Source(I,L,h,x,y,z,t,fine);
 [hz_01,hz_impulse,hx_01,hx_impulse,hy_01,hy_impulse,ex_01,ex_impulse,ey_01,ey_impulse,ez_01,ez_impulse] =...
    Calculate_Horizontal_Electrical_Dipole_SemiAirborne(I,L,h,x,y,z,t);

%% E
% Ex
figure;
loglog(t(1:end).*10^3,abs((ex_01)),'r','Linewidth',2);
hold on
legend('��ֵ��ex\_01');
title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ') ex step response'])
xlabel('Time/(ms)')
ylabel('Ex/(V/m)');
% Ey
figure;
loglog(t(1:end).*10^3,abs((ey_01)),'r','Linewidth',2);
hold on
legend('��ֵ��ey\_01');
title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ') ey step response'])
xlabel('Time/(ms)')
ylabel('Ey/(V/m)');
% Ez
figure;
loglog(t(1:end).*10^3,abs((ez_01)),'r','Linewidth',2);
hold on
legend('��ֵ��ez\_01');
title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')ez step response'])
xlabel('Time/(ms)')
ylabel('Ez/(V/m)');
%% H
% Hx
figure;
loglog(t(1:end).*10^3,abs((hx_01)),'r','Linewidth',2);
hold on
legend('��ֵ��Hx\_01');
title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ') Hx step response'])
xlabel('Time/(ms)')
ylabel('Hx/(V/m)');
% Hy
figure;
loglog(t(1:end).*10^3,abs((hy_01)),'r','Linewidth',2);
hold on
legend('��ֵ��Hy\_01');
title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ') Hy step response'])
xlabel('Time/(ms)')
ylabel('Hy/(V/m)');
% Hz
figure;
loglog(t(1:end).*10^3,abs((hz_01)),'r','Linewidth',2);
hold on
legend('��ֵ��Hz\_01');
title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')Hz step response'])
xlabel('Time/(ms)')
ylabel('Hz/(V/m)');
