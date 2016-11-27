%% ���Ƶȸ���ͼ��ʸ��ͼ
clc;close all;
x = (-200:3:200); % �շ�ˮƽƫ�ƾ࣬��y��
y = (-400:3:400);
[X,Y] = meshgrid(x,y);% x��ֵΪx�����꣬y��ֵΪy������
%% EM contourf
% Ex
figure; 
C_Ex = contourf(X,Y,Ex_3D(:,:,1,1),'ShowText','on');% ��ǵȸ��ߵ���ֵ
xlabel('X position / (m)');
ylabel('Y position / (m)');
title('the contour of the Ex');
clabel(C_Ex);% ��ǵȸ�����ֵ
colorbar;% ��ʾͼƬ�Աߵ���ɫ��
% Ey
figure; 
C_Ey = contourf(X,Y,Ey_3D(:,:,1,1),'ShowText','on');% ��ǵȸ��ߵ���ֵ
xlabel('X position / (m)');
ylabel('Y position / (m)');
title('the contour of the Ey');
clabel(C_Ey);% ��ǵȸ�����ֵ
colorbar;% ��ʾͼƬ�Աߵ���ɫ��
%
% Hx
figure; 
C_Hx = contourf(X,Y,Hx_3D(:,:,1,1),'ShowText','on');% ��ǵȸ��ߵ���ֵ
xlabel('X position / (m)');
ylabel('Y position / (m)');
title('the contour of the Hx');
clabel(C_Hx);% ��ǵȸ�����ֵ
colorbar;% ��ʾͼƬ�Աߵ���ɫ��
% Hy
figure; 
C_Hy = contourf(X,Y,Hy_3D(:,:,1,1),'ShowText','on');% ��ǵȸ��ߵ���ֵ
xlabel('X position / (m)');
ylabel('Y position / (m)');
title('the contour of the Hy');
clabel(C_Hy);% ��ǵȸ�����ֵ
colorbar;% ��ʾͼƬ�Աߵ���ɫ��
% Hz
figure; 
C_Hz = contourf(X,Y,Hz_3D(:,:,1,1),'ShowText','on');% ��ǵȸ��ߵ���ֵ
xlabel('X position / (m)');
ylabel('Y position / (m)');
title('the contour of the Hz');
clabel(C_Hz);% ��ǵȸ�����ֵ
colorbar;% ��ʾͼƬ�Աߵ���ɫ��

%% EM quiver
% E
figure; 
Q_E = quiver(X,Y,Ex_3D(:,:,1,1),Ey_3D(:,:,1,1));% ��ǵȸ��ߵ���ֵ
xlabel('X position / (m)');
ylabel('Y position / (m)');
title('the quiver of the E');
colorbar;% ��ʾͼƬ�Աߵ���ɫ��
% E
figure; 
Q_H = quiver(X,Y,Hx_3D(:,:,1,1),Hy_3D(:,:,1,1));% ��ǵȸ��ߵ���ֵ
xlabel('X position / (m)');
ylabel('Y position / (m)');
title('the quiver of the H');
colorbar;% ��ʾͼƬ�Աߵ���ɫ��