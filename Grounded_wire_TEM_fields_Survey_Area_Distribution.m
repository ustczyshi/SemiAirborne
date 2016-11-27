%--------------------------------------------------------------------------
%仿真条件：地面TEM 半航空TEM
% 地面水平接地长导线源激励，地面或空中观测 B场
% 1. 均匀半空间电阻率设定100欧姆米,先观测不同高度响应的变化，主要是观察Hz,Hx,Hy,高度设定0,30,60,90米
% 2. 改变长导线源中垂线方向偏移距的长度（y=100,500,1000,2000），观察不同偏移距下，响应随高度的变化
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
x = (-200:40:200); % 收发水平偏移距，沿y轴
y = (-200:40:200);
[X,Y] = meshgrid(x,y);% x的值为x轴坐标，y的值为y轴坐标
L= 10; % 发射线缆长度，沿x轴
I = 0.1; % 发射电流
%%  半航空收发高度参数
% z =[0, -30,-60,-90];% 观测点距地面的高度，地面以上为负值
 z =[-50,-100];
h =0;% 源距地面的高度
n = 12;
%% 采样率和观测时间段设置
% fs = 1e5;% 采样率
% dt = 1./fs;
% t = logspace(-3,2,50);% 时间区间
% t = [1e-5 1e-4 1e-3 1e-2 1e-1 1];
t = 1e-4;
%%
%    1-D: 行数（对应坐标是的y轴）;2-D:列数(对应坐标系的x轴);3-D: 高度（对应z轴）；4-D：时间信息(对应不同的观测时间点)
Ex_3D = zeros(length(y),length(x),length(z),length(t));
Ey_3D = zeros(length(y),length(x),length(z),length(t));
Ez_3D = zeros(length(y),length(x),length(z),length(t));

Hz_3D = zeros(length(y),length(x),length(z),length(t));
Hx_3D = zeros(length(y),length(x),length(z),length(t));
Hy_3D = zeros(length(y),length(x),length(z),length(t));

Uz_3D = zeros(length(y),length(x),length(z),length(t));
Uy_3D = zeros(length(y),length(x),length(z),length(t));
Ux_3D = zeros(length(y),length(x),length(z),length(t));

tic;
for kt = 1:length(t)
    for kz = 1:length(z) %the third Dimension
        for ky = 1:length(y) % row
            for kx = 1:length(x) % col
        
%         [hz_01,hz_10,hz_impulse,hx_01,hx_10,hx_impulse,hy_01,hy_10,hy_impulse,ex_01,ex_10,ex_impulse,ey_01,ey_10,ey_impulse] =...
%             Calculate_Horizontal_Finite_Electrical_Source_unsave(I,L,h,x(kx),y(ky),z(kz),t,fine);
 [hz_01,hz_10,hz_impulse,hx_01,hx_10,hx_impulse,hy_01,hy_10,hy_impulse,ex_01,ex_10,ex_impulse,ey_01,ey_10,ey_impulse,ez_01,ez_impulse] = ...
    Calculate_Horizontal_Finite_Electrical_Source_GuaLeg(I,L,h,x(kx),y(ky),z(kz),t(kt),n)
        %
        Hz_3D(ky,kx,kz,kt) = hz_01;
        Hy_3D(ky,kx,kz,kt) = hy_01;
        Hx_3D(ky,kx,kz,kt) = hx_01;
        
        Ex_3D(ky,kx,kz,kt) = ex_01;
        Ey_3D(ky,kx,kz,kt) = ey_01;
        Ez_3D(ky,kx,kz,kt) = ez_01;
        
        Uz_3D(ky,kx,kz,kt) = hz_impulse;
        Uy_3D(ky,kx,kz,kt) = hy_impulse;
        Ux_3D(ky,kx,kz,kt) = hx_impulse;
            end
        end
    end
end
t = toc
%% U
%
% 
U_max  = u0.*max( max( max( max([Ux_3D Uy_3D Uz_3D]) ) ) ).*1e9 ;% 转化为nT/s
U_min  = u0.*min( min( min( min([Ux_3D Uy_3D Uz_3D]) ) ) ).*1e9  ;% 转化为nT/s
figure; 
subplot(1,3,1);
C_Ux = contourf(X,Y,u0.*(Ux_3D(:,:,2,1)).*1e9 ,10,'ShowText','on');% 标记等高线的数值
xlabel('X position / (m)');
ylabel('Y position / (m)');
title('the contour of the Ux');
% clabel(C_Hx);% 标记等高线数值
set(gca,'CLim',[U_min,U_max]);% 统一不同子图的colorbar
c = colorbar;
c.Label.String = 'nT/s/m^2';
% Uy
subplot(1,3,2);
C_Uy = contourf(X,Y,u0.*(Uy_3D(:,:,2,1)).*1e9 ,10,'ShowText','on');% 标记等高线的数值
xlabel('X position / (m)');
ylabel('Y position / (m)');
title('the contour of the Uy');
% clabel(C_Hy);% 标记等高线数值
set(gca,'CLim',[U_min,U_max]);
c = colorbar;
c.Label.String = 'nT/s/m^2';
% Uz
subplot(1,3,3);
C_Uz = contourf(X,Y,u0.*(Uz_3D(:,:,2,1)).*1e9 ,10,'ShowText','on');% 标记等高线的数值
xlabel('X position / (m)');
ylabel('Y position / (m)');
title('the contour of the Uz');
% clabel(C_Hz);% 标记等高线数值
set(gca,'CLim',[U_min,U_max]);
c = colorbar;% 显示图片旁边的颜色棒
c.Label.String = 'nT/s/m^2'; % 标注颜色棒值单位
%
%% EM contourf
% Ex
E_max  = max( max( max( max( abs([Ex_3D Ey_3D Ez_3D]) )) ) ).*1e3 ;% 转化为mV/m
E_min  = min( min( min( min( abs([Ex_3D Ey_3D Ez_3D]) ) ) ) ).*1e3   ;% 转化为mV/m
figure; 
subplot(1,3,1);
C_Ex = contourf(X,Y,abs(Ex_3D(:,:,2,1).*1e3 ),6,'ShowText','on');% 标记等高线的数值
xlabel('X position / (m)');
ylabel('Y position / (m)');
title('the contour of the Ex');
% clabel(C_Ex);% 标记等高线数值
set(gca,'CLim',[E_min,E_max]);
c = colorbar;% 显示图片旁边的颜色棒
c.Label.String = 'mV/m'; % 标注颜色棒值单位
% Ey
subplot(1,3,2);
C_Ey = contourf(X,Y,abs(Ey_3D(:,:,2,1)).*1e3 ,6,'ShowText','on');% 标记等高线的数值
xlabel('X position / (m)');
ylabel('Y position / (m)');
title('the contour of the Ey');
% clabel(C_Ey);% 标记等高线数值
set(gca,'CLim',[E_min,E_max]);
c = colorbar;% 显示图片旁边的颜色棒
c.Label.String = 'mV/m'; % 标注颜色棒值单位
% Ez
subplot(1,3,3);
C_Ez = contourf(X,Y,abs(Ez_3D(:,:,2,1)).*1e3 ,6,'ShowText','on');% 标记等高线的数值
xlabel('X position / (m)');
ylabel('Y position / (m)');
title('the contour of the Ez');
% clabel(C_Ey);% 标记等高线数值
set(gca,'CLim',[E_min,E_max]);
c = colorbar;% 显示图片旁边的颜色棒
c.Label.String = 'mV/m'; % 标注颜色棒值单位
%%
% Hx

H_max  = u0.*max( max( max( max( abs([Hx_3D Hy_3D Hz_3D]) )) ) ).*1e9 ;% 转化为nT
H_min  = u0.*min( min( min( min( abs([Hx_3D Hy_3D Hz_3D]) ) ) ) ).*1e9   ;% 转化为nT
figure; 
subplot(1,3,1);
C_Hx = contourf(X,Y,u0.*(Hx_3D(:,:,2,1)).*1e9,10,'ShowText','on');% 标记等高线的数值
xlabel('X position / (m)');
ylabel('Y position / (m)');
title('the contour of the Bx');
% clabel(C_Hx);% 标记等高线数值
set(gca,'CLim',[H_min,H_max]);
c = colorbar;% 显示图片旁边的颜色棒
c.Label.String = 'nT'; % 标注颜色棒值单位
%  Hy
subplot(1,3,2);
C_Hy = contourf(X,Y,u0.*(Hy_3D(:,:,2,1)).*1e9,10,'ShowText','on');% 标记等高线的数值
xlabel('X position / (m)');
ylabel('Y position / (m)');
title('the contour of the By');
% clabel(C_Hy);% 标记等高线数值
set(gca,'CLim',[H_min,H_max]);
c = colorbar;% 显示图片旁边的颜色棒
c.Label.String = 'nT'; % 标注颜色棒值单位
% Hz
subplot(1,3,3);
C_Hz = contourf(X,Y,u0.*(Hz_3D(:,:,2,1)).*1e9,10,'ShowText','on');% 标记等高线的数值
xlabel('X position / (m)');
ylabel('Y position / (m)');
title('the contour of the Bz');
% clabel(C_Hz);% 标记等高线数值
set(gca,'CLim',[H_min,H_max]);
c = colorbar;% 显示图片旁边的颜色棒
c.Label.String = 'nT'; % 标注颜色棒值单位

%% EM quiver
% E
figure; 
Q_E = quiver(X,Y,abs(Ex_3D(:,:,2,1)),abs(Ey_3D(:,:,2,1)) );% 标记等高线的数值
hold on;
streamslice(X,Y,abs( Ex_3D(:,:,2,1) ), abs( Ey_3D(:,:,2,1)) );
xlabel('X position / (m)');
ylabel('Y position / (m)');
title('the quiver of the E');
% colorbar;% 显示图片旁边的颜色棒
hold off;
%% H quiver
figure; 
Q_H = quiver(X,Y,Hx_3D(:,:,2,1),Hy_3D(:,:,2,1));% 标记等高线的数值
% streamline(X,Y,Hx_3D(:,:,1,1),Hy_3D(:,:,1,1));
% C_Ex = contour(X,Y,(Hx_3D(:,:,2,1).^2+Hy_3D(:,:,2,1).^2).^0.5,20,'ShowText','on');% 标记等高线的数值
% C_Ex = contour(X,Y,(Hx_3D(:,:,2,1).^2+Hy_3D(:,:,2,1).^2).^0.5,20);% 标记等高线的数值
hold on;
streamslice(X,Y,Hx_3D(:,:,2,1),Hy_3D(:,:,2,1));
xlabel('X position / (m)');
ylabel('Y position / (m)');
title('the quiver of the H');
% colorbar;% 显示图片旁边的颜色棒
hold off;



