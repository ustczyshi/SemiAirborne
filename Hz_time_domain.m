%--------------------------------------------------------------------------
%仿真条件：仿真瞬断垂直磁偶极的垂直磁场随时间变化规律。
%磁偶极置于电阻率100 Ohm*m 的均匀大地表面，观察点距离偶极100m。
%--------------------------------------------------------------------------
%%
clear all;clc;close all;
format long;
%%
mu_0 = 4*pi*1e-7;
%% atem 收发高度参数
h = 100;% 源距地面的高度
z = 50;% 观测点距地面的高度
%% 地面
r=100; % 收发距
G_S=load ('G_S.txt')';
m2 = 1:length(G_S);
fs = 1e7;
t = 1/fs:1/fs:1e-2;% 时间区间
% t = logspace(-8,-1,1000);
h_z1_t=zeros(1,length(t));
h_z2_t=h_z1_t;
%  for ii=1:length(t)
%     freq = (log(2)*1i/(t(ii)*2*pi))*m2;
    %--------------------------------------------------------------------------
    %第一步：读取已经存储的滤波器系数,表示为行向量；
    load J0_Gupt.txt;       
    J_zero = J0_Gupt( :, 3)'; % 快速汉克尔变换滤波系数
    delta = J0_Gupt( :, 2)'; %  采样点的横坐标偏移量
    %--------------------------------------------------------------------------
    %计算lambda，并将lambda和frequency扩展成二维矩阵
    lambda = (1./r) .*exp(delta); % 如何计算lambda，由采样点的横坐标偏移量转换为积分变量lambda
for ii=1:length(t)
    freq = (log(2)*1i/(t(ii)*2*pi))*m2;
    [lambda_Array,frequency_Array] = Array_trans(lambda,freq);
    %--------------------------------------------------------------------------
    %根据递推公式计算空气-地面的反射系数
    r_TE=calculate_r_TE(lambda,freq);
    %% --------------------------------------------------------------------------
    %计算磁场的垂直分量，并用快速汉克尔变换求积分的数值解。  
    %垂直磁偶源z轴磁场的核函数,h是源位置,z为观测位置
    % 垂直磁偶极源z轴磁场的核函数 u_0 = lambda,z=h=0,地面发射地面接收;
    sum = (1+r_TE).*lambda_Array.^2;  % 正阶跃响应频域核函数
    sum_zeros = (1).*lambda_Array.^2; %% 对应直流成分，此时的核函数中的r_TE=0
    H_vertical = 1./(4*pi) *  Fast_Hankel(r,sum,J_zero);%正阶跃响应频域
    H_vertical_zeros = 1./(4*pi) *  Fast_Hankel(r,sum_zeros,J_zero);%正阶跃响应零频响应
    h_z1_t(ii) = -GS_Trans2(t(ii),H_vertical,G_S);%GS_Trans2(t(ii),H_vertical_zeros); %负脉冲响应时域
    h_z2_t(ii) = +GS_Trans(t(ii),H_vertical_zeros,freq,G_S)-GS_Trans(t(ii),H_vertical,freq,G_S);%负阶跃响应时域
    %% ---------------------------------------------------------------------------
    % 半航空垂直磁偶源z轴磁场的核函数,h是源位置,z为观测位置.对于半航空,z=0,h=-100米。
    % sum = exp(-lambda_Array*h).*(1+r_TE).*lambda_Array.^2;
    %% -----------------------------------atem 垂直方向----------------------------------------
    % 全航空atem垂直磁偶源z轴磁场的核函数,h是源位置,z为观测位置.对于全航空,z=-50,h=-100米。
    sum = (exp(-lambda_Array*(z+h))+r_TE.*exp(lambda_Array*(z-h))).*lambda_Array.^2;
    sum_zeros = (1).*lambda_Array.^2; %% 对应直流成分，此时的核函数中的r_TE=0
    H_vertical = 1./(4*pi) *  Fast_Hankel(r,sum,J_zero);%正阶跃响应频域
    H_vertical_zeros = 1./(4*pi) *  Fast_Hankel(r,sum_zeros,J_zero);%正阶跃响应零频响应
    h_z1_t(ii) = -GS_Trans2(t(ii),H_vertical,G_S);%GS_Trans2(t(ii),H_vertical_zeros); %负脉冲响应时域
    h_z2_t(ii) = +GS_Trans(t(ii),H_vertical_zeros,freq,G_S)-GS_Trans(t(ii),H_vertical,freq,G_S);%负阶跃响应时域
    %% -----------------------------------atem 径向----------------------------------------
    %% -----------------------------------atem phi方向-----------------------------------
 end
%% --------------------------------------------------------------------------
%作图。
loglog(t.*10^3,(h_z1_t),'r','Linewidth',2)
hold on
% axis([10^-5,10^3,10^-11,10^-1])
loglog(t.*10^3,0.5*(abs(h_z1_t)-h_z1_t),'r--','Linewidth',2)
hold on
% axis([10^-5,10^3,10^-11,10^-1])
loglog(t.*10^3,(h_z2_t),'b','Linewidth',2)
hold on 
% axis([10^-5,10^3,10^-11,10^-1])
loglog(t.*10^3,0.5*(abs(h_z2_t)-h_z2_t),'b--','Linewidth',2)
hold on 
% axis([10^-5,10^3,10^-11,10^-1])
legend('负脉冲响应','','负阶跃响应','');
title('瞬断垂直磁偶极子的垂直磁场及其时间导数随时间变化的曲线')
xlabel('时间（ms）')
%% 保存负脉冲响应
himpulse10 = h_z1_t;
hstep10 = h_z2_t;
save('atem_magnetic_dipole_response.mat','himpulse10','hstep10');
%%
