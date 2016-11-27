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
y = 1000;
L= 2500; % 发射线缆长度，沿x轴
I = 40; % 发射电流
m = I*L;
%%  半航空收发高度参数
z =-100;% 观测点距地面的高度，地面以上为负值
h =0;% 源距地面的高度
n=12;
%% 采样率和观测时间段设置
fs = 1e7;% 采样率
dt = 1./fs;
% t = 1/fs:1/fs:2e-3;% 时间区间
t = logspace(-5,0,100);
%%
[step_ex01,step_ey01,impulse_hx,impulse_hy,impulse_hz] = Horizontal_Electrical_Dipole_Response(u0,rou,m,x,y,z,t);
% [hz_01,hz_10,hz_impulse,hx_01,hx_10,hx_impulse,hy_01,hy_10,hy_impulse,ex_01,ex_impulse,ey_01,ey_impulse]...
% = Calculate_Horizontal_Electrical_Dipole(I,L,h,x,y,z,t);
% tic;
% fine = 0.4;
% [hz_01,hz_10,hz_impulse,hx_01,hx_10,hx_impulse,hy_01,hy_10,hy_impulse,ex_01,ex_10,ex_impulse,ey_01,ey_10,ey_impulse] = ...
%     Calculate_Horizontal_Finite_Electrical_Source(I,L,h,x,y,z,t,fine);% 利用高斯-勒让德积分的原理求解长接地导线源的瞬变电磁场
% t0 = toc
tic;
[hz_01,hz_10,hz_impulse,hx_01,hx_10,hx_impulse,hy_01,hy_10,hy_impulse,ex_01,ex_10,ex_impulse,ey_01,ey_10,ey_impulse] = ...
    Calculate_Horizontal_Finite_Electrical_Source_GuaLeg(I,L,h,x,y,z,t,n);% 利用高斯-勒让德积分的原理求解长接地导线源的瞬变电磁场
t1 = toc
save('horizontal_electrical_dipole_impulse_shuzhijie','hx_impulse','hy_impulse','hz_impulse');
%% z轴阶跃
%{
%    作图。
hz_10 = hz_01(end) - hz_01;
figure;
plot(t.*10^3,u0.*hz_01,'r','Linewidth',1);
hold on
plot(t.*10^3,u0.*hz_10,'b','Linewidth',1);
grid on;
legend('数值解正阶跃响应','数值解负阶跃响应');
title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop Bz'])
xlabel('Time/(ms)')
ylabel('Bz/(A/m)');
% x轴阶跃
%作图。
% hx_10 = hx_01(end) - hx_01;
figure;
plot(t.*10^3,u0.*(hx_01),'r','Linewidth',1);
hold on
plot(t.*10^3,-u0.*(hx_01),'b','Linewidth',1);
grid on;
legend('数值解正阶跃响应','数值解负阶跃响应');
title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop Bx'])
xlabel('Time/(ms)')
ylabel('Bx/(A/m)');
% y 轴阶跃
hy_10 = hy_01(end) - hy_01;
figure;
plot(t.*10^3,u0.*(hy_01),'r','Linewidth',1);
hold on
plot(t.*10^3,u0.*(hy_10),'b','Linewidth',1);
grid on;
legend('数值解正阶跃响应','数值解负阶跃响应');
title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop By'])
xlabel('Time/(ms)')
ylabel('By/(A/m)');
%}
%% z 轴脉冲
jiexijie = ['Horizontal_Electrical_Dipole_Response_rou' num2str(rou) '_xr' num2str(x) '_yr' num2str(y) '_zr' num2str(z) '.mat'];
load(jiexijie);
figure;
loglog(t(1:end).*10^3,abs((ex_10))./L,'r','Linewidth',2);
hold on
loglog(t(1:end).*10^3,abs(step_ex01(end)-step_ex01)./L,'k:','Linewidth',2);
grid on;
legend('数值解ex\_01','解析解step\_ex01');
title(['source moment' num2str(I) 'A*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop ex step response'])
xlabel('Time/(ms)')
ylabel('Ex/(V/m)');
figure;
loglog(t(1:end).*10^3,abs((ey_10))./L,'r','Linewidth',2);
hold on
loglog(t(1:end).*10^3,abs(step_ey01(end)-step_ey01)./L,'k:','Linewidth',2);
grid on;
legend('数值解ey\_01','解析解step\_ey01');
title(['source moment' num2str(I) 'A*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop ey step response'])
xlabel('Time/(ms)')
ylabel('Ey/(V/m)');
% 误差分析
% ex_error = abs(abs(ex_01)-abs(step_ex01))./abs(step_ex01);
% ey_error = abs(abs(ey_01)-abs(step_ey01))./abs(step_ey01);
% figure;
% loglog(t.*1e3,ex_error.*100,'r','linewidth',2);
% hold on;
% loglog(t.*1e3,ey_error.*100,'k:','linewidth',2);
% grid on;
% legend('ex','ey');
% title('数值解和解析解的误差分析');
% xlabel('Time/(ms)')
% ylabel('error/(%)');
% 脉冲响应
% figure;
% plot(t(1:end).*10^3,ex_impulse,'r','Linewidth',2);
% hold on
% plot(t(1:end-1).*10^3,diff(step_ex01)./dt,'k:','Linewidth',2);
% grid on;
% legend('数值解ex\_01','解析解step\_ex01');
% title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop ex step response'])
% xlabel('Time/(ms)')
% ylabel('Ex/(V/m)');
% figure;
% plot(t(1:end).*10^3,abs(ey_01),'r','Linewidth',2);
% hold on
% plot(t(1:end).*10^3,abs(step_ey01),'k:','Linewidth',2);
% grid on;
% legend('数值解ey\_01','解析解step\_ey01');
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
legend('数值解hz\_1\_impulse','数值解diff(hz\_01)./dt','解析解impulse\_hz');
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
legend('数值解hx\_1\_impulse','数值解diff(hx\_01)./dt','解析解impulse\_hx');
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
legend('数值解hy\_1\_impulse','数值解diff(hy\_01)./dt','解析解impulse\_hy');
title(['source moment' num2str(I) 'A*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop By impulse response'])
xlabel('Time/(ms)')
ylabel('By/(T)');
% 误差分析
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
title('数值解和解析解的误差分析');
xlabel('Time/(ms)')
ylabel('error/(%)');

%% save data

