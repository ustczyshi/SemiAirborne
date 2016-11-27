%%  �о�ˮƽ��żԴ dBz/dt ��dBy/dt�ڲ�ͬ�絼���µ���Ӧǿ��
%   ���ø���Ծ��Ӧ����
%  ������2�� x�᣺��ͬ�ĵ絼��
%                        y�᣺��ͬ��ƫ�ƾ�
%                        ��ɫ��Uz��һ������ʱ������Ӧǿ��   
%  ������1�� x�᣺��ͬ�ĵ絼��
%                      y�᣺dB/dt
%%
clc;close all;clear all;
%
format long;
%%
mu0 = 4*pi*1e-7;
sigma = logspace(-3,2,100);
rou = 1./sigma;
I = 40;
ds = 1;
M = I*ds;

t = [1/48e3,1e-4,1e-3,1e-2];
x = 0;
y = linspace(50,300,20);
r = (x.^2+y.^2).^(0.5);
% �м�������� results1
n1 = length(sigma);
n2 = length(r);
Hz_1 = zeros(n2,n1);
parpool; % �������г�
parfor k1 = 1:n1  % �絼�ʣ�x��
    for k2 = 1:n2
        u = r(k2).*(mu0.*sigma(k1)./4./t(1));
        Hz_1(k2,k1) =M./(2.*pi.*sigma(k1).*mu0).*y(k2)./(r(k2).^5).*(3.*erf(u) -2./(sqrt(pi)) .*u.*(3+2.*u^2).*exp(-u.^2) );  
    end
end
%%
delete(gcp); % �رղ��г�


%% display the result1
[X,Y] =meshgrid(rou,r); 
figure;
surf(X,Y,mu0.*abs(Hz_1));
title([ '$\frac{dBz}{dt}$' '@first time window'],'interpreter','latex');
xlabel('resistivity(\Omega\cdotm)');
ylabel('offset(m)');
% colorbar;
% colormap summer;
%%
figure;
image((rou),(r),(abs(Hz_1)));
title([ '$\frac{dHz}{dy}$' '@first time window'],'interpreter','latex');
xlabel('resistivity(\Omega\cdotm)');
ylabel('offset(m)');
colorbar;
%%
figure;
image((r),(rou),(abs(Hz_1)));
title([ '$\frac{dHz}{dy}$' '@first time window'],'interpreter','latex');
xlabel('resistivity(\Omega\cdotm)');
ylabel('offset(m)');
colorbar;
%% result2
nt = length(t);
Uz_2 = zeros(nt,n1);
parpool;
parfor k1 = 1:n1
    for kt = 1:nt
        u = r(8).*(mu0.*sigma(k1)./4./t(kt));
        Uz_2(kt,k1) =M./(2.*pi.*sigma(k1)).*y(8)./(r(8).^5).*(3.*erf(u) -2./(sqrt(pi)) .*u.*(3+2.*u^2).*exp(-u.^2) );  
    end
end
delete(gcp);
%% display
figure;
loglog(rou,abs(Uz_2));
xlabel('resistivity(\Omega\cdotm)');
ylabel('dBz/dt');






