%% 绘制等高线图和矢量图
clc;close all;
x = (-200:3:200); % 收发水平偏移距，沿y轴
y = (-400:3:400);
[X,Y] = meshgrid(x,y);% x的值为x轴坐标，y的值为y轴坐标
%% EM contourf
% Ex
figure; 
C_Ex = contourf(X,Y,Ex_3D(:,:,1,1),'ShowText','on');% 标记等高线的数值
xlabel('X position / (m)');
ylabel('Y position / (m)');
title('the contour of the Ex');
clabel(C_Ex);% 标记等高线数值
colorbar;% 显示图片旁边的颜色棒
% Ey
figure; 
C_Ey = contourf(X,Y,Ey_3D(:,:,1,1),'ShowText','on');% 标记等高线的数值
xlabel('X position / (m)');
ylabel('Y position / (m)');
title('the contour of the Ey');
clabel(C_Ey);% 标记等高线数值
colorbar;% 显示图片旁边的颜色棒
%
% Hx
figure; 
C_Hx = contourf(X,Y,Hx_3D(:,:,1,1),'ShowText','on');% 标记等高线的数值
xlabel('X position / (m)');
ylabel('Y position / (m)');
title('the contour of the Hx');
clabel(C_Hx);% 标记等高线数值
colorbar;% 显示图片旁边的颜色棒
% Hy
figure; 
C_Hy = contourf(X,Y,Hy_3D(:,:,1,1),'ShowText','on');% 标记等高线的数值
xlabel('X position / (m)');
ylabel('Y position / (m)');
title('the contour of the Hy');
clabel(C_Hy);% 标记等高线数值
colorbar;% 显示图片旁边的颜色棒
% Hz
figure; 
C_Hz = contourf(X,Y,Hz_3D(:,:,1,1),'ShowText','on');% 标记等高线的数值
xlabel('X position / (m)');
ylabel('Y position / (m)');
title('the contour of the Hz');
clabel(C_Hz);% 标记等高线数值
colorbar;% 显示图片旁边的颜色棒

%% EM quiver
% E
figure; 
Q_E = quiver(X,Y,Ex_3D(:,:,1,1),Ey_3D(:,:,1,1));% 标记等高线的数值
xlabel('X position / (m)');
ylabel('Y position / (m)');
title('the quiver of the E');
colorbar;% 显示图片旁边的颜色棒
% E
figure; 
Q_H = quiver(X,Y,Hx_3D(:,:,1,1),Hy_3D(:,:,1,1));% 标记等高线的数值
xlabel('X position / (m)');
ylabel('Y position / (m)');
title('the quiver of the H');
colorbar;% 显示图片旁边的颜色棒