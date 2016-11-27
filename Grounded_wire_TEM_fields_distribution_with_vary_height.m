%--------------------------------------------------------------------------
%��������������TEM �뺽��TEM
% ����ˮƽ�ӵس�����Դ�������������й۲� B��
% 1. ���Ȱ�ռ�������趨100ŷķ��,�ȹ۲ⲻͬ�߶���Ӧ�ı仯����Ҫ�ǹ۲�Hz,Hx,Hy,�߶��趨0,30,60,90��
% 2. �ı䳤����Դ�д��߷���ƫ�ƾ�ĳ��ȣ�y=100,500,1000,2000�����۲첻ͬƫ�ƾ��£���Ӧ��߶ȵı仯
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
x = (-2500:200:2500); % �շ�ˮƽƫ�ƾ࣬��y��
y = (-2500:200:2500);
% x = 100;
% y =  1500;
[X,Y] = meshgrid(x,y);% x��ֵΪx�����꣬y��ֵΪy������
% y = [500,1000,2000];
L= 2500; % �������³��ȣ���x��
I = 40; % �������
fine = 0.4;% ϸ�����ӣ��������񻮷ִ�С,һ��ѡ��0.5��1
n = 12;
%%  �뺽���շ��߶Ȳ���
% z =[0, -30,-60,-90];% �۲������ĸ߶ȣ���������Ϊ��ֵ
z = [0 -50 -100 -150] ;
h =0;% Դ�����ĸ߶�

%% �����ʺ͹۲�ʱ�������
% fs = 1e5;% ������
% dt = 1./fs;
% t = logspace(-3,2,50);% ʱ������
 t = 1e-1;

%%
%    1-D: ��������Ӧ�����ǵ�y�ᣩ;2-D:����(��Ӧ����ϵ��x��);3-D: �߶ȣ���Ӧz�ᣩ��4-D��ʱ����Ϣ(��Ӧ��ͬ�Ĺ۲�ʱ���)
Ex_3D = zeros(length(y),length(x),length(z),length(t));
Ey_3D = zeros(length(y),length(x),length(z),length(t));

Hz_3D = zeros(length(y),length(x),length(z),length(t));
Hx_3D = zeros(length(y),length(x),length(z),length(t));
Hy_3D = zeros(length(y),length(x),length(z),length(t));

Uz_3D = zeros(length(y),length(x),length(z),length(t));
Uy_3D = zeros(length(y),length(x),length(z),length(t));
Ux_3D = zeros(length(y),length(x),length(z),length(t));
for kt = 1:length(t)
    
for kz = 1:length(z) %the third Dimension
    
for ky = 1:length(y) % row
    
    for kx = 1:length(x) % col
        
%         [hz_01,hz_10,hz_impulse,hx_01,hx_10,hx_impulse,hy_01,hy_10,hy_impulse,ex_01,ex_10,ex_impulse,ey_01,ey_10,ey_impulse] =...
%             Calculate_Horizontal_Finite_Electrical_Source_unsave(I,L,h,x(kx),y(ky),z(kz),t,fine);
        [hz_01,hz_10,hz_impulse,hx_01,hx_10,hx_impulse,hy_01,hy_10,hy_impulse,ex_01,ex_10,ex_impulse,ey_01,ey_10,ey_impulse] =...
            Calculate_Horizontal_Finite_Electrical_Source_GuaLeg(I,L,h,x(kx),y(ky),z(kz),t,n);
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
    
end

end

end
%% U
%
% Ux
figure; 
C_Ux = contourf(X,Y,u0.*(Ux_3D(:,:,2,1)),20,'ShowText','on');% ��ǵȸ��ߵ���ֵ
xlabel('X position / (m)');
ylabel('Y position / (m)');
title('the contour of the Ux');
% clabel(C_Hx);% ��ǵȸ�����ֵ
colorbar;% ��ʾͼƬ�Աߵ���ɫ��
% Uy
figure; 
C_Uy = contourf(X,Y,u0.*(Uy_3D(:,:,2,1)),20,'ShowText','on');% ��ǵȸ��ߵ���ֵ
xlabel('X position / (m)');
ylabel('Y position / (m)');
title('the contour of the Uy');
% clabel(C_Hy);% ��ǵȸ�����ֵ
colorbar;% ��ʾͼƬ�Աߵ���ɫ��
% Uz
figure; 
C_Uz = contourf(X,Y,u0.*(Uz_3D(:,:,2,1)),20,'ShowText','on');% ��ǵȸ��ߵ���ֵ
xlabel('X position / (m)');
ylabel('Y position / (m)');
title('the contour of the Uz');
% clabel(C_Hz);% ��ǵȸ�����ֵ
colorbar;% ��ʾͼƬ�Աߵ���ɫ��
%
%% EM contourf
% Ex
figure; 
C_Ex = contourf(X,Y,Ex_3D(:,:,2,1),20,'ShowText','on');% ��ǵȸ��ߵ���ֵ
xlabel('X position / (m)');
ylabel('Y position / (m)');
title('the contour of the Ex');
% clabel(C_Ex);% ��ǵȸ�����ֵ
colorbar;% ��ʾͼƬ�Աߵ���ɫ��
% Ey
figure; 
C_Ey = contourf(X,Y,Ey_3D(:,:,2,1),20,'ShowText','on');% ��ǵȸ��ߵ���ֵ
xlabel('X position / (m)');
ylabel('Y position / (m)');
title('the contour of the Ey');
% clabel(C_Ey);% ��ǵȸ�����ֵ
colorbar;% ��ʾͼƬ�Աߵ���ɫ��
%%
% Hx
figure; 
C_Hx = contourf(X,Y,(Hx_3D(:,:,2,1)),20,'ShowText','on');% ��ǵȸ��ߵ���ֵ
xlabel('X position / (m)');
ylabel('Y position / (m)');
title('the contour of the Hx');
% clabel(C_Hx);% ��ǵȸ�����ֵ
colorbar;% ��ʾͼƬ�Աߵ���ɫ��
%  Hy
figure; 
C_Hy = contourf(X,Y,(Hy_3D(:,:,2,1)),20,'ShowText','on');% ��ǵȸ��ߵ���ֵ
xlabel('X position / (m)');
ylabel('Y position / (m)');
title('the contour of the Hy');
% clabel(C_Hy);% ��ǵȸ�����ֵ
colorbar;% ��ʾͼƬ�Աߵ���ɫ��
% Hz
figure; 
C_Hz = contourf(X,Y,(Hz_3D(:,:,2,1)),20,'ShowText','on');% ��ǵȸ��ߵ���ֵ
xlabel('X position / (m)');
ylabel('Y position / (m)');
title('the contour of the Hz');
% clabel(C_Hz);% ��ǵȸ�����ֵ
colorbar;% ��ʾͼƬ�Աߵ���ɫ��

%% EM quiver
% E
figure; 
Q_E = quiver(X,Y,Ex_3D(:,:,2,1),Ey_3D(:,:,2,1));% ��ǵȸ��ߵ���ֵ
hold on;
streamslice(X,Y,Ex_3D(:,:,2,1),Ey_3D(:,:,2,1));
xlabel('X position / (m)');
ylabel('Y position / (m)');
title('the quiver of the E');
colorbar;% ��ʾͼƬ�Աߵ���ɫ��
hold off;
%% H quiver
figure; 
Q_H = quiver(X,Y,Hx_3D(:,:,2,1),Hy_3D(:,:,2,1));% ��ǵȸ��ߵ���ֵ
% streamline(X,Y,Hx_3D(:,:,1,1),Hy_3D(:,:,1,1));
% C_Ex = contour(X,Y,(Hx_3D(:,:,2,1).^2+Hy_3D(:,:,2,1).^2).^0.5,20,'ShowText','on');% ��ǵȸ��ߵ���ֵ
% C_Ex = contour(X,Y,(Hx_3D(:,:,2,1).^2+Hy_3D(:,:,2,1).^2).^0.5,20);% ��ǵȸ��ߵ���ֵ
hold on;
streamslice(X,Y,Hx_3D(:,:,2,1),Hy_3D(:,:,2,1));
xlabel('X position / (m)');
ylabel('Y position / (m)');
title('the quiver of the H');
colorbar;% ��ʾͼƬ�Աߵ���ɫ��
hold off;
%


