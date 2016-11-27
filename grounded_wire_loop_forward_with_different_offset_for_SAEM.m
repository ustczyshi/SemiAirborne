%% 考察 不同电阻率均匀大地的响应变化
%--------------------------------------------------------------------------
%仿真条件：地面TEM 半航空TEM
% 地面水平接地长导线源激励，地面或空中观测 B场
%--------------------------------------------------------------------------
%%
format long;
clear all;clc;close all;
%--------------------------------------------------------------------------
%仿真条件：地面TEM 半航空TEM
% 地面水平接地长导线源激励，空中观测 B场
%--------------------------------------------------------------------------
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
x = 0; % 收发水平偏移距，沿y轴
y = [100,500,1000,1500,2000,2500];
[X,Y] = meshgrid(x,y);% x的值为x轴坐标，y的值为y轴坐标
L= 2500; % 发射线缆长度，沿x轴
I = 40; % 发射电流
%%  半航空收发高度参数
% z =[0, -30,-60,-90];% 观测点距地面的高度，地面以上为负值
 z =[-50,-100];
h =0;% 源距地面的高度
n = 8;
%% 采样率和观测时间段设置
% fs = 1e5;% 采样率
% dt = 1./fs;
t = logspace(-5,-1,50);% 时间区间
% t = [1e-5 1e-4 1e-3 1e-2 1e-1 1];
% t = 1e-3;
tol = 1e-8;
% 确定长导线的分点位置
[Ak,xk,dxk] = GuaLeg_DiscreteSource_Out(L,n,tol);
% Ak : 积分系数，列向量；
% xk : 积分节点，实际是n阶勒让德多项式的n各节点，为列向量；
% dxk : 列向量，存放分成的长导线各段的长度；
[row,col ] = size(X);% row = length(y)
offset = zeros(row,col,n);
for k = 1:n
    offset(:,:,k) =sqrt((X-xk(k)).^2+Y.^2); 
end

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
% for kt = 1:length(t)  
    for kz = 1:length(z) %the third Dimension
        for ky = 1:length(y) % row
            for kx = 1:length(x) % col
%         [hz_01,hz_10,hz_impulse,hx_01,hx_10,hx_impulse,hy_01,hy_10,hy_impulse,ex_01,ex_10,ex_impulse,ey_01,ey_10,ey_impulse] =...
%             Calculate_Horizontal_Finite_Electrical_Source_unsave(I,L,h,x(kx),y(ky),z(kz),t,fine);
 [hz_01,hz_10,hz_impulse,hx_01,hx_10,hx_impulse,hy_01,hy_10,hy_impulse,ex_01,ex_10,ex_impulse,ey_01,ey_10,ey_impulse,ez_01,ez_impulse] = ...
    Calculate_Horizontal_Finite_Electrical_Source_GuaLeg_out(I,L,h,x(kx),y(ky),z(kz),offset(ky,kx,:),Ak,t,n);
        %
        Hz_3D(ky,kx,kz,:) = hz_10;
        Hy_3D(ky,kx,kz,:) = hy_10;
        Hx_3D(ky,kx,kz,:) = hx_10;
        
        Ex_3D(ky,kx,kz,:) = ex_01;
        Ey_3D(ky,kx,kz,:) = ey_01;
        Ez_3D(ky,kx,kz,:) = ez_01;
        
        Uz_3D(ky,kx,kz,:) = hz_impulse;
        Uy_3D(ky,kx,kz,:) = hy_impulse;
        Ux_3D(ky,kx,kz,:) = hx_impulse;
            end
%     delete(gcp('nocreate'));
        end
    end
% end
tt = toc
%% display;
% Hz
figure;
Hz  = [reshape(Hz_3D(1,1,2,:),[],1) reshape(Hz_3D(2,1,2,:),[],1) reshape(Hz_3D(3,1,2,:),[],1) reshape(Hz_3D(4,1,2,:),[],1)  ...
    reshape(Hz_3D(5,1,2,:),[],1) reshape(Hz_3D(6,1,2,:),[],1)] ;
loglog(t,u0.*abs(Hz),'Linewidth',2);
hold on;
title('Bz');
% title(['Grounded Finite wire L=' num2str(L) ' offset=' num2str(y) ' rou =' num2str(rou)]);
legend([  num2str(y(1)) 'm'],[ num2str(y(2)) 'm'],...
    [ num2str(y(3)) 'm'],[ num2str(y(4)) 'm'],[ num2str(y(5)) 'm'],[ num2str(y(6)) 'm']);
legend boxoff;
grid on;
xlabel('Time/s');
ylabel('Bz/(T)');
figure;
Hy  = [reshape(Hy_3D(1,1,2,:),[],1) reshape(Hy_3D(2,1,2,:),[],1) reshape(Hy_3D(3,1,2,:),[],1) reshape(Hy_3D(4,1,2,:),[],1) ...
    reshape(Hy_3D(5,1,2,:),[],1) reshape(Hy_3D(6,1,2,:),[],1)] ;
loglog(t,u0.*abs(Hy),'Linewidth',2);
hold on;
title('By');
% title(['Grounded Finite wire L=' num2str(L) ' offset=' num2str(y) ' rou =' num2str(rou)]);
legend([  num2str(y(1)) 'm'],[ num2str(y(2)) 'm'],...
    [ num2str(y(3)) 'm'],[ num2str(y(4)) 'm'],[ num2str(y(5)) 'm'],[ num2str(y(6)) 'm']);
legend boxoff;
grid on;
xlabel('Time/s');
ylabel('By/(T)');

figure;
Hx  = [reshape(Hx_3D(1,1,2,:),[],1) reshape(Hx_3D(2,1,2,:),[],1) reshape(Hx_3D(3,1,2,:),[],1) reshape(Hx_3D(4,1,2,:),[],1) ...
    reshape(Hx_3D(5,1,2,:),[],1) reshape(Hx_3D(6,1,2,:),[],1)] ;
loglog(t,u0.*abs(Hx),'Linewidth',2);
hold on;
title('Bx');
% title(['Grounded Finite wire L=' num2str(L) ' offset=' num2str(y) ' rou =' num2str(rou)]);
legend([  num2str(y(1)) 'm'],[ num2str(y(2)) 'm'],...
    [ num2str(y(3)) 'm'],[ num2str(y(4)) 'm'],[ num2str(y(5)) 'm'],[ num2str(y(6)) 'm']);
legend boxoff;
grid on;
xlabel('Time/s');
ylabel('Bx/(T)');
%%  Uz
figure;
Uz  = [reshape(Uz_3D(1,1,2,:),[],1) reshape(Uz_3D(2,1,2,:),[],1) reshape(Uz_3D(3,1,2,:),[],1) reshape(Uz_3D(4,1,2,:),[],1)  ...
    reshape(Uz_3D(5,1,2,:),[],1) reshape(Uz_3D(6,1,2,:),[],1)] ;
loglog(t,u0.*abs(Uz),'Linewidth',2);%./repmat(diff(t),length(y),1)
hold on;
hold on;
title('dBz/dt');
% title(['Grounded Finite wire L=' num2str(L) ' offset=' num2str(y) ' rou =' num2str(rou)]);
legend([  num2str(y(1)) 'm'],[ num2str(y(2)) 'm'],...
    [ num2str(y(3)) 'm'],[ num2str(y(4)) 'm'],[ num2str(y(5)) 'm'],[ num2str(y(6)) 'm']);
legend boxoff;
grid on;
xlabel('Time/s');
ylabel('dBz/dt(V/m^2)');

figure;
Uy  = [reshape(Uy_3D(1,1,2,:),[],1) reshape(Uy_3D(2,1,2,:),[],1) reshape(Uy_3D(3,1,2,:),[],1) reshape(Uy_3D(4,1,2,:),[],1)  ...
    reshape(Uy_3D(5,1,2,:),[],1) reshape(Uy_3D(6,1,2,:),[],1)] ;
loglog(t,u0.*abs(Uy),'Linewidth',2);%./repmat(diff(t),length(y),1)
hold on;
hold on;
title('dBy/dt');
% title(['Grounded Finite wire L=' num2str(L) ' offset=' num2str(y) ' rou =' num2str(rou)]);
legend([  num2str(y(1)) 'm'],[ num2str(y(2)) 'm'],...
    [ num2str(y(3)) 'm'],[ num2str(y(4)) 'm'],[ num2str(y(5)) 'm'],[ num2str(y(6)) 'm']);
legend boxoff;
grid on;
xlabel('Time/s');
ylabel('dBy/dt(V/m^2)');

figure;
Ux  = [reshape(Ux_3D(1,1,2,:),[],1) reshape(Ux_3D(2,1,2,:),[],1) reshape(Ux_3D(3,1,2,:),[],1) reshape(Ux_3D(4,1,2,:),[],1)  ...
    reshape(Ux_3D(5,1,2,:),[],1) reshape(Ux_3D(6,1,2,:),[],1)] ;
loglog(t,u0.*abs(Ux),'Linewidth',2);%./repmat(diff(t),length(y),1)
hold on;
hold on;
title('dBx/dt');
% title(['Grounded Finite wire L=' num2str(L) ' offset=' num2str(y) ' rou =' num2str(rou)]);
legend([  num2str(y(1)) 'm'],[ num2str(y(2)) 'm'],...
    [ num2str(y(3)) 'm'],[ num2str(y(4)) 'm'],[ num2str(y(5)) 'm'],[ num2str(y(6)) 'm']);
legend boxoff;
grid on;
xlabel('Time/s');
ylabel('dBx/dt(V/m^2)');
