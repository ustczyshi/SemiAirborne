%% ���� ��ͬ�����ʾ��ȴ�ص���Ӧ�仯
%--------------------------------------------------------------------------
%��������������TEM �뺽��TEM
% ����ˮƽ�ӵس�����Դ�������������й۲� B��
%--------------------------------------------------------------------------
%%
format long;
clear all;clc;close all;
%--------------------------------------------------------------------------
%��������������TEM �뺽��TEM
% ����ˮƽ�ӵس�����Դ���������й۲� B��
%--------------------------------------------------------------------------
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
x = 0; % �շ�ˮƽƫ�ƾ࣬��y��
y = [100,500,1000,1500,2000,2500];
[X,Y] = meshgrid(x,y);% x��ֵΪx�����꣬y��ֵΪy������
L= 2500; % �������³��ȣ���x��
I = 40; % �������
%%  �뺽���շ��߶Ȳ���
% z =[0, -30,-60,-90];% �۲������ĸ߶ȣ���������Ϊ��ֵ
 z =[-50,-100];
h =0;% Դ�����ĸ߶�
n = 8;
%% �����ʺ͹۲�ʱ�������
% fs = 1e5;% ������
% dt = 1./fs;
t = logspace(-5,-1,50);% ʱ������
% t = [1e-5 1e-4 1e-3 1e-2 1e-1 1];
% t = 1e-3;
tol = 1e-8;
% ȷ�������ߵķֵ�λ��
[Ak,xk,dxk] = GuaLeg_DiscreteSource_Out(L,n,tol);
% Ak : ����ϵ������������
% xk : ���ֽڵ㣬ʵ����n�����õ¶���ʽ��n���ڵ㣬Ϊ��������
% dxk : ����������ŷֳɵĳ����߸��εĳ��ȣ�
[row,col ] = size(X);% row = length(y)
offset = zeros(row,col,n);
for k = 1:n
    offset(:,:,k) =sqrt((X-xk(k)).^2+Y.^2); 
end

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
