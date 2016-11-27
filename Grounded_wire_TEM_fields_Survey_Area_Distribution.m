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
x = (-200:40:200); % �շ�ˮƽƫ�ƾ࣬��y��
y = (-200:40:200);
[X,Y] = meshgrid(x,y);% x��ֵΪx�����꣬y��ֵΪy������
L= 10; % �������³��ȣ���x��
I = 0.1; % �������
%%  �뺽���շ��߶Ȳ���
% z =[0, -30,-60,-90];% �۲������ĸ߶ȣ���������Ϊ��ֵ
 z =[-50,-100];
h =0;% Դ�����ĸ߶�
n = 12;
%% �����ʺ͹۲�ʱ�������
% fs = 1e5;% ������
% dt = 1./fs;
% t = logspace(-3,2,50);% ʱ������
% t = [1e-5 1e-4 1e-3 1e-2 1e-1 1];
t = 1e-4;
%%
%    1-D: ��������Ӧ�����ǵ�y�ᣩ;2-D:����(��Ӧ����ϵ��x��);3-D: �߶ȣ���Ӧz�ᣩ��4-D��ʱ����Ϣ(��Ӧ��ͬ�Ĺ۲�ʱ���)
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
U_max  = u0.*max( max( max( max([Ux_3D Uy_3D Uz_3D]) ) ) ).*1e9 ;% ת��ΪnT/s
U_min  = u0.*min( min( min( min([Ux_3D Uy_3D Uz_3D]) ) ) ).*1e9  ;% ת��ΪnT/s
figure; 
subplot(1,3,1);
C_Ux = contourf(X,Y,u0.*(Ux_3D(:,:,2,1)).*1e9 ,10,'ShowText','on');% ��ǵȸ��ߵ���ֵ
xlabel('X position / (m)');
ylabel('Y position / (m)');
title('the contour of the Ux');
% clabel(C_Hx);% ��ǵȸ�����ֵ
set(gca,'CLim',[U_min,U_max]);% ͳһ��ͬ��ͼ��colorbar
c = colorbar;
c.Label.String = 'nT/s/m^2';
% Uy
subplot(1,3,2);
C_Uy = contourf(X,Y,u0.*(Uy_3D(:,:,2,1)).*1e9 ,10,'ShowText','on');% ��ǵȸ��ߵ���ֵ
xlabel('X position / (m)');
ylabel('Y position / (m)');
title('the contour of the Uy');
% clabel(C_Hy);% ��ǵȸ�����ֵ
set(gca,'CLim',[U_min,U_max]);
c = colorbar;
c.Label.String = 'nT/s/m^2';
% Uz
subplot(1,3,3);
C_Uz = contourf(X,Y,u0.*(Uz_3D(:,:,2,1)).*1e9 ,10,'ShowText','on');% ��ǵȸ��ߵ���ֵ
xlabel('X position / (m)');
ylabel('Y position / (m)');
title('the contour of the Uz');
% clabel(C_Hz);% ��ǵȸ�����ֵ
set(gca,'CLim',[U_min,U_max]);
c = colorbar;% ��ʾͼƬ�Աߵ���ɫ��
c.Label.String = 'nT/s/m^2'; % ��ע��ɫ��ֵ��λ
%
%% EM contourf
% Ex
E_max  = max( max( max( max( abs([Ex_3D Ey_3D Ez_3D]) )) ) ).*1e3 ;% ת��ΪmV/m
E_min  = min( min( min( min( abs([Ex_3D Ey_3D Ez_3D]) ) ) ) ).*1e3   ;% ת��ΪmV/m
figure; 
subplot(1,3,1);
C_Ex = contourf(X,Y,abs(Ex_3D(:,:,2,1).*1e3 ),6,'ShowText','on');% ��ǵȸ��ߵ���ֵ
xlabel('X position / (m)');
ylabel('Y position / (m)');
title('the contour of the Ex');
% clabel(C_Ex);% ��ǵȸ�����ֵ
set(gca,'CLim',[E_min,E_max]);
c = colorbar;% ��ʾͼƬ�Աߵ���ɫ��
c.Label.String = 'mV/m'; % ��ע��ɫ��ֵ��λ
% Ey
subplot(1,3,2);
C_Ey = contourf(X,Y,abs(Ey_3D(:,:,2,1)).*1e3 ,6,'ShowText','on');% ��ǵȸ��ߵ���ֵ
xlabel('X position / (m)');
ylabel('Y position / (m)');
title('the contour of the Ey');
% clabel(C_Ey);% ��ǵȸ�����ֵ
set(gca,'CLim',[E_min,E_max]);
c = colorbar;% ��ʾͼƬ�Աߵ���ɫ��
c.Label.String = 'mV/m'; % ��ע��ɫ��ֵ��λ
% Ez
subplot(1,3,3);
C_Ez = contourf(X,Y,abs(Ez_3D(:,:,2,1)).*1e3 ,6,'ShowText','on');% ��ǵȸ��ߵ���ֵ
xlabel('X position / (m)');
ylabel('Y position / (m)');
title('the contour of the Ez');
% clabel(C_Ey);% ��ǵȸ�����ֵ
set(gca,'CLim',[E_min,E_max]);
c = colorbar;% ��ʾͼƬ�Աߵ���ɫ��
c.Label.String = 'mV/m'; % ��ע��ɫ��ֵ��λ
%%
% Hx

H_max  = u0.*max( max( max( max( abs([Hx_3D Hy_3D Hz_3D]) )) ) ).*1e9 ;% ת��ΪnT
H_min  = u0.*min( min( min( min( abs([Hx_3D Hy_3D Hz_3D]) ) ) ) ).*1e9   ;% ת��ΪnT
figure; 
subplot(1,3,1);
C_Hx = contourf(X,Y,u0.*(Hx_3D(:,:,2,1)).*1e9,10,'ShowText','on');% ��ǵȸ��ߵ���ֵ
xlabel('X position / (m)');
ylabel('Y position / (m)');
title('the contour of the Bx');
% clabel(C_Hx);% ��ǵȸ�����ֵ
set(gca,'CLim',[H_min,H_max]);
c = colorbar;% ��ʾͼƬ�Աߵ���ɫ��
c.Label.String = 'nT'; % ��ע��ɫ��ֵ��λ
%  Hy
subplot(1,3,2);
C_Hy = contourf(X,Y,u0.*(Hy_3D(:,:,2,1)).*1e9,10,'ShowText','on');% ��ǵȸ��ߵ���ֵ
xlabel('X position / (m)');
ylabel('Y position / (m)');
title('the contour of the By');
% clabel(C_Hy);% ��ǵȸ�����ֵ
set(gca,'CLim',[H_min,H_max]);
c = colorbar;% ��ʾͼƬ�Աߵ���ɫ��
c.Label.String = 'nT'; % ��ע��ɫ��ֵ��λ
% Hz
subplot(1,3,3);
C_Hz = contourf(X,Y,u0.*(Hz_3D(:,:,2,1)).*1e9,10,'ShowText','on');% ��ǵȸ��ߵ���ֵ
xlabel('X position / (m)');
ylabel('Y position / (m)');
title('the contour of the Bz');
% clabel(C_Hz);% ��ǵȸ�����ֵ
set(gca,'CLim',[H_min,H_max]);
c = colorbar;% ��ʾͼƬ�Աߵ���ɫ��
c.Label.String = 'nT'; % ��ע��ɫ��ֵ��λ

%% EM quiver
% E
figure; 
Q_E = quiver(X,Y,abs(Ex_3D(:,:,2,1)),abs(Ey_3D(:,:,2,1)) );% ��ǵȸ��ߵ���ֵ
hold on;
streamslice(X,Y,abs( Ex_3D(:,:,2,1) ), abs( Ey_3D(:,:,2,1)) );
xlabel('X position / (m)');
ylabel('Y position / (m)');
title('the quiver of the E');
% colorbar;% ��ʾͼƬ�Աߵ���ɫ��
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
% colorbar;% ��ʾͼƬ�Աߵ���ɫ��
hold off;



