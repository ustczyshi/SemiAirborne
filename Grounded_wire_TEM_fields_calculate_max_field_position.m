%--------------------------------------------------------------------------
%仿真条件：地面TEM 半航空TEM
% 地面水平接地长导线源激励，地面或空中观测 B场和dB/dt
% 1. 均匀半空间电阻率设定100欧姆米,主要是观察Hz,Hx,Hy,高度设定0,50,100,150米
% 2. 改变长导线源中垂线方向偏移距的长度（y=100,500,1000,2000），观察不同偏移距下，响应随高度的变化
% 3.得到第一个采样时窗时刻dBz/dt最大值所在的偏移距位置
% 4.半航空系统的采样率为48kHz(电流记录和响应记录)
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
% x = (-2500:200:2500); % 收发水平偏移距，沿y轴
x = 0;
% y = (0:50:2500);
 y = (0:25:1000);
% x = 100;
% y =  1500;
[X,Y] = meshgrid(x,y);% x的值为x轴坐标，y的值为y轴坐标
% y = [500,1000,2000];
L= 2500; % 发射线缆长度，沿x轴
I = 40; % 发射电流
fine = 0.4;% 细化因子，决定网格划分大小,一般选择0.5或1
n = 12;
%%  半航空收发高度参数
% z =[0, -30,-60,-90];% 观测点距地面的高度，地面以上为负值
z = [0 -50 -100 -150 ] ;
h =0;% 源距地面的高度

%% 采样率和观测时间段设置
% fs = 1e5;% 采样率
% dt = 1./fs;
% t = logspace(-3,2,50);% 时间区间
 t = [1/48e3 1e-1];
% t = 0.1;

%%
%    1-D: 行数（对应坐标是的y轴）;2-D:列数(对应坐标系的x轴);3-D: 高度（对应z轴）；4-D：时间信息(对应不同的观测时间点)
Ex_3D = zeros(length(y),length(x),length(z),length(t));
Ey_3D = zeros(length(y),length(x),length(z),length(t));

Hz_3D = zeros(length(y),length(x),length(z),length(t));
Hx_3D = zeros(length(y),length(x),length(z),length(t));
Hy_3D = zeros(length(y),length(x),length(z),length(t));

Uz_3D = zeros(length(y),length(x),length(z),length(t));
Uy_3D = zeros(length(y),length(x),length(z),length(t));
Ux_3D = zeros(length(y),length(x),length(z),length(t));


for kx = 1:length(x) % col
for kt = 1:length(t)
    
for kz = 1:length(z) %the third Dimension
    parpool;
    parfor ky = 1:length(y) % row
 
%         [hz_01,hz_10,hz_impulse,hx_01,hx_10,hx_impulse,hy_01,hy_10,hy_impulse,ex_01,ex_10,ex_impulse,ey_01,ey_10,ey_impulse] =...
%             Calculate_Horizontal_Finite_Electrical_Source_unsave(I,L,h,x(kx),y(ky),z(kz),t,fine);
        [hz_01,hz_10,hz_impulse,hx_01,hx_10,hx_impulse,hy_01,hy_10,hy_impulse,ex_01,ex_10,ex_impulse,ey_01,ey_10,ey_impulse] =...
            Calculate_Horizontal_Finite_Electrical_Source_GuaLeg(I,L,h,x(kx),y(ky),z(kz),t(kt),n);
        %
        Hz_3D(ky,kx,kz,kt) = hz_10;
        Hy_3D(ky,kx,kz,kt) = hy_10;
        
        Hx_3D(ky,kx,kz,kt) = hx_10;
        Ex_3D(ky,kx,kz,kt) = -ex_01;
        Ey_3D(ky,kx,kz,kt) = -ey_01;
        
        Uz_3D(ky,kx,kz,kt) = hz_impulse;
        Uy_3D(ky,kx,kz,kt) = hy_impulse;
        Ux_3D(ky,kx,kz,kt) = hx_impulse;
            %
            %{
        Hz_3D(ky,kx,kz,kt) = hz_impulse;
        Hy_3D(ky,kx,kz,kt) = hy_impulse;
        Hx_3D(ky,kx,kz,kt) = hx_impulse;
        Ex_3D(ky,kx,kz,kt) = -ex_01;
        Ey_3D(ky,kx,kz,kt) = -ey_01;
          %}  
    end
delete(gcp);   
end

end

end

%% U
% 绘制Uz vs 偏移距的曲线@ first time window
figure;
plot(y,abs(u0.*Uz_3D(:,1,1,1)),'r','linewidth',2);
hold on;
plot(y,abs(u0.*Uz_3D(:,1,2,1)),'b','linewidth',2);
hold on;
plot(y,abs(u0.*Uz_3D(:,1,3,1)),'k','linewidth',2);
hold on;
plot(y,abs(u0.*Uz_3D(:,1,4,1)),'m','linewidth',2);
hold on;
grid on;
title(['dBz/dt field @first time window('  num2str(t(1)*1e3) 'ms )']);
legend(['h=' num2str(-z(1)) 'm'],['h=' num2str(-z(2)) 'm'],['h=' num2str(-z(3)) 'm'],['h=' num2str(-z(4)) 'm'] );
legend boxoff;
xlabel('offset(m)');
ylabel('dBz/dt(V/m^2)');
% 绘制Uz vs 偏移距的曲线@@ first time window
figure;
plot(y,abs(u0.*Uy_3D(:,1,1,1)),'r','linewidth',2);
hold on;
plot(y,abs(u0.*Uy_3D(:,1,2,1)),'b','linewidth',2);
hold on;
plot(y,abs(u0.*Uy_3D(:,1,3,1)),'k','linewidth',2);
hold on;
plot(y,abs(u0.*Uy_3D(:,1,4,1)),'m','linewidth',2);
hold on;
grid on;
title(['dBy/dt field @first time window('  num2str(t(1)*1e3) 'ms )']);
legend(['h=' num2str(-z(1)) 'm'],['h=' num2str(-z(2)) 'm'],['h=' num2str(-z(3)) 'm'],['h=' num2str(-z(4)) 'm'] );
legend boxoff;
xlabel('offset(m)');
ylabel('dBy/dt(V/m^2)');

% 绘制Uz vs 偏移距的曲线@ last time window
figure;
plot(y,abs(u0.*Uz_3D(:,1,1,2)),'r','linewidth',2);
hold on;
plot(y,abs(u0.*Uz_3D(:,1,2,2)),'b','linewidth',2);
hold on;
plot(y,abs(u0.*Uz_3D(:,1,3,2)),'k','linewidth',2);
hold on;
plot(y,abs(u0.*Uz_3D(:,1,4,2)),'m','linewidth',2);
hold on;
grid on;
title(['dBz/dt field @last time window('  num2str(t(2)*1e3) 'ms )']);
legend(['h=' num2str(-z(1)) 'm'],['h=' num2str(-z(2)) 'm'],['h=' num2str(-z(3)) 'm'],['h=' num2str(-z(4)) 'm'] );
legend boxoff;
xlabel('offset(m)');
ylabel('dBz/dt(V/m^2)');
% 绘制Uz vs 偏移距的曲线@@ last time window
figure;
plot(y,abs(u0.*Uy_3D(:,1,1,2)),'r','linewidth',2);
hold on;
plot(y,abs(u0.*Uy_3D(:,1,2,2)),'b','linewidth',2);
hold on;
plot(y,abs(u0.*Uy_3D(:,1,3,2)),'k','linewidth',2);
hold on;
plot(y,abs(u0.*Uy_3D(:,1,4,2)),'m','linewidth',2);
hold on;
grid on;
title(['dBy/dt field @last time window('  num2str(t(2)*1e3) 'ms )']);
legend(['h=' num2str(-z(1)) 'm'],['h=' num2str(-z(2)) 'm'],['h=' num2str(-z(3)) 'm'],['h=' num2str(-z(4)) 'm'] );
legend boxoff;
xlabel('offset(m)');
ylabel('dBy/dt(V/m^2)');
%%  B
%% B
% 绘制Uz vs 偏移距的曲线@ first time window
figure;
plot(y,abs(u0.*Hz_3D(:,1,1,1)),'r','linewidth',2);
hold on;
plot(y,abs(u0.*Hz_3D(:,1,2,1)),'b','linewidth',2);
hold on;
plot(y,abs(u0.*Hz_3D(:,1,3,1)),'k','linewidth',2);
hold on;
plot(y,abs(u0.*Hz_3D(:,1,4,1)),'m','linewidth',2);
hold on;
grid on;
title(['Bz field @first time window('  num2str(t(1)*1e3) 'ms )']);
legend(['h=' num2str(-z(1)) 'm'],['h=' num2str(-z(2)) 'm'],['h=' num2str(-z(3)) 'm'],['h=' num2str(-z(4)) 'm'] );
legend boxoff;
xlabel('offset(m)');
ylabel('Bz(T)');
% 绘制Uz vs 偏移距的曲线@@ first time window
figure;
plot(y,abs(u0.*Hy_3D(:,1,1,1)),'r','linewidth',2);
hold on;
plot(y,abs(u0.*Hy_3D(:,1,2,1)),'b','linewidth',2);
hold on;
plot(y,abs(u0.*Hy_3D(:,1,3,1)),'k','linewidth',2);
hold on;
plot(y,abs(u0.*Hy_3D(:,1,4,1)),'m','linewidth',2);
hold on;
grid on;
title(['By field @first time window('  num2str(t(1)*1e3) 'ms )']);
legend(['h=' num2str(-z(1)) 'm'],['h=' num2str(-z(2)) 'm'],['h=' num2str(-z(3)) 'm'],['h=' num2str(-z(4)) 'm'] );
legend boxoff;
xlabel('offset(m)');
ylabel('By(T)');

% 绘制Uz vs 偏移距的曲线@ last time window
figure;
plot(y,abs(u0.*Hz_3D(:,1,1,2)),'r','linewidth',2);
hold on;
plot(y,abs(u0.*Hz_3D(:,1,2,2)),'b','linewidth',2);
hold on;
plot(y,abs(u0.*Hz_3D(:,1,3,2)),'k','linewidth',2);
hold on;
plot(y,abs(u0.*Hz_3D(:,1,4,2)),'m','linewidth',2);
hold on;
grid on;
title(['Bz field @last time window('  num2str(t(2)*1e3) 'ms )']);
legend(['h=' num2str(-z(1)) 'm'],['h=' num2str(-z(2)) 'm'],['h=' num2str(-z(3)) 'm'],['h=' num2str(-z(4)) 'm'] );
legend boxoff;
xlabel('offset(m)');
ylabel('Bz(T)');
% 绘制Uz vs 偏移距的曲线@@ last time window
figure;
plot(y,abs(u0.*Hy_3D(:,1,1,2)),'r','linewidth',2);
hold on;
plot(y,abs(u0.*Hy_3D(:,1,2,2)),'b','linewidth',2);
hold on;
plot(y,abs(u0.*Hy_3D(:,1,3,2)),'k','linewidth',2);
hold on;
plot(y,abs(u0.*Hy_3D(:,1,4,2)),'m','linewidth',2);
hold on;
grid on;
title(['By field @last time window('  num2str(t(2)*1e3) 'ms )']);
legend(['h=' num2str(-z(1)) 'm'],['h=' num2str(-z(2)) 'm'],['h=' num2str(-z(3)) 'm'],['h=' num2str(-z(4)) 'm'] );
legend boxoff;
xlabel('offset(m)');
ylabel('By(T)');





%{
% Ux
figure; 
C_Ux = contourf(X,Y,u0.*(Ux_3D(:,:,2,1)),20,'ShowText','on');% 标记等高线的数值
xlabel('X position / (m)');
ylabel('Y position / (m)');
title('the contour of the Ux');
% clabel(C_Hx);% 标记等高线数值
colorbar;% 显示图片旁边的颜色棒
% Uy
figure; 
C_Uy = contourf(X,Y,u0.*(Uy_3D(:,:,2,1)),20,'ShowText','on');% 标记等高线的数值
xlabel('X position / (m)');
ylabel('Y position / (m)');
title('the contour of the Uy');
% clabel(C_Hy);% 标记等高线数值
colorbar;% 显示图片旁边的颜色棒
% Uz
figure; 
C_Uz = contourf(X,Y,u0.*(Uz_3D(:,:,2,1)),20,'ShowText','on');% 标记等高线的数值
xlabel('X position / (m)');
ylabel('Y position / (m)');
title('the contour of the Uz');
% clabel(C_Hz);% 标记等高线数值
colorbar;% 显示图片旁边的颜色棒
%
%% EM contourf
% Ex
figure; 
C_Ex = contourf(X,Y,Ex_3D(:,:,2,1),20,'ShowText','on');% 标记等高线的数值
xlabel('X position / (m)');
ylabel('Y position / (m)');
title('the contour of the Ex');
% clabel(C_Ex);% 标记等高线数值
colorbar;% 显示图片旁边的颜色棒
% Ey
figure; 
C_Ey = contourf(X,Y,Ey_3D(:,:,2,1),20,'ShowText','on');% 标记等高线的数值
xlabel('X position / (m)');
ylabel('Y position / (m)');
title('the contour of the Ey');
% clabel(C_Ey);% 标记等高线数值
colorbar;% 显示图片旁边的颜色棒
%%
% Hx
figure; 
C_Hx = contourf(X,Y,(Hx_3D(:,:,2,1)),20,'ShowText','on');% 标记等高线的数值
xlabel('X position / (m)');
ylabel('Y position / (m)');
title('the contour of the Hx');
% clabel(C_Hx);% 标记等高线数值
colorbar;% 显示图片旁边的颜色棒
%  Hy
figure; 
C_Hy = contourf(X,Y,(Hy_3D(:,:,2,1)),20,'ShowText','on');% 标记等高线的数值
xlabel('X position / (m)');
ylabel('Y position / (m)');
title('the contour of the Hy');
% clabel(C_Hy);% 标记等高线数值
colorbar;% 显示图片旁边的颜色棒
% Hz
figure; 
C_Hz = contourf(X,Y,(Hz_3D(:,:,2,1)),20,'ShowText','on');% 标记等高线的数值
xlabel('X position / (m)');
ylabel('Y position / (m)');
title('the contour of the Hz');
% clabel(C_Hz);% 标记等高线数值
colorbar;% 显示图片旁边的颜色棒

%% EM quiver
% E
figure; 
Q_E = quiver(X,Y,Ex_3D(:,:,2,1),Ey_3D(:,:,2,1));% 标记等高线的数值
hold on;
streamslice(X,Y,Ex_3D(:,:,2,1),Ey_3D(:,:,2,1));
xlabel('X position / (m)');
ylabel('Y position / (m)');
title('the quiver of the E');
colorbar;% 显示图片旁边的颜色棒
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
colorbar;% 显示图片旁边的颜色棒
hold off;
%}


