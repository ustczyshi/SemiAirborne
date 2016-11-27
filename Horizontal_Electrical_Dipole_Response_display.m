%--------------------------------------------------------------------------
%仿真条件：地面TEM 半航空TEM
% 地面水平电偶源激励，地面或空中观测

%--------------------------------------------------------------------------
%%
format long;
clear all;clc;close all;
%%
u0 = 4*pi*1e-7;
load parameters.txt;
sigma1 = parameters(1,2);%第一层的电导率
rou = 1./sigma1;
%% 发射机参数地面
x = 0; % 收发水平偏移距，沿y轴
y = 0;
L= 1; % 发射线缆长度，沿x轴
I = 1; % 发射电流
m = I*L;
%%  半航空收发高度参数
z =-100;% 观测点距地面的高度，地面以上为负值
h =0;% 源距地面的高度

%% 采样率和观测时间段设置
fs = 1e7;% 采样率
dt = 1./fs;
% t = 1/fs:1/fs:2e-3;% 时间区间
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
legend('数值解ex\_01');
title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ') ex step response'])
xlabel('Time/(ms)')
ylabel('Ex/(V/m)');
% Ey
figure;
loglog(t(1:end).*10^3,abs((ey_01)),'r','Linewidth',2);
hold on
legend('数值解ey\_01');
title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ') ey step response'])
xlabel('Time/(ms)')
ylabel('Ey/(V/m)');
% Ez
figure;
loglog(t(1:end).*10^3,abs((ez_01)),'r','Linewidth',2);
hold on
legend('数值解ez\_01');
title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')ez step response'])
xlabel('Time/(ms)')
ylabel('Ez/(V/m)');
%% H
% Hx
figure;
loglog(t(1:end).*10^3,abs((hx_01)),'r','Linewidth',2);
hold on
legend('数值解Hx\_01');
title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ') Hx step response'])
xlabel('Time/(ms)')
ylabel('Hx/(V/m)');
% Hy
figure;
loglog(t(1:end).*10^3,abs((hy_01)),'r','Linewidth',2);
hold on
legend('数值解Hy\_01');
title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ') Hy step response'])
xlabel('Time/(ms)')
ylabel('Hy/(V/m)');
% Hz
figure;
loglog(t(1:end).*10^3,abs((hz_01)),'r','Linewidth',2);
hold on
legend('数值解Hz\_01');
title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')Hz step response'])
xlabel('Time/(ms)')
ylabel('Hz/(V/m)');
