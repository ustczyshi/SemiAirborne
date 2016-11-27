%% the sensitivity of the hx-hy-hz
%%
%  分析磁场各分量对电阻率及层厚的敏感性
%  求误差的rms值
clc;close all;clear all;
%% load parameters and set parameters
u0 = 4*pi*1e-7;
load parameters.txt; % add the variable  to the workplace
parameters0 = parameters;
parametersD = [1 0.01 100;0 1/50 0];
sigma1 = parameters(1,2);%第一层的电导率
rou_0 = 1./sigma1;
h0 = parameters(1,3);
% set resistivity of  the model
rate_r  =  logspace(-5,2,35);
rou_a = rou_0.*rate_r;
sigma_a = 1./rou_a;
% set the thickness of the anomal object
rate_h = logspace(-2,1,10);
h_a = h0.*rate_h;
%% set parameters
% Tx
L= 2000; % 发射线缆长度，沿x轴
I = 40; % 发射电流
% mesh
fine = 0.5;% 细化因子，决定网格划分大小,一般选择0.5或1
% time 
t  = logspace(-5,1,50);
%%  receive position set
x = 500; % 收发水平偏移距，沿y轴
y = 1000;
z =0;% 观测点距地面的高度，地面以上为负值
h =0;% 源距地面的高度
%% rebuild the parameters.m
% parfor parallel compute
rlen = length(rou_a);
hlen = length(h_a);
tlen = length(t);
hzh_10_mat = zeros(hlen,tlen);
hyh_10_mat = zeros(hlen,tlen);
hxh_10_mat = zeros(hlen,tlen);
% change the thickness
for kh = 1:hlen
        parameters(2,3) = h_a(kh) ;  
        save('parameters.txt','parameters','-ascii');
        
         [hz_01,hzh_10_mat(kh,:),hz_impulse,hx_01,hxh_10_mat(kh,:),hx_impulse,hy_01,hyh_10_mat(kh,:),hy_impulse,ex_01,ex_10,ex_impulse,ey_01,ey_10,ey_impulse] =...
            Calculate_Horizontal_Finite_Electrical_Source_unsave(I,L,h,x,y,z,t,fine);      
end
% change the resistivity
hzr_10_mat = zeros(rlen,tlen);
hyr_10_mat = zeros(rlen,tlen);
hxr_10_mat = zeros(rlen,tlen);
% save('parameters.txt','parameters0','-ascii');
parameters = parameters0;
for kr = 1:rlen
    parameters(2,2) = sigma_a(kr) ; 
    save('parameters.txt','parameters','-ascii');
    
    [hz_01,hzr_10_mat(kr,:),hz_impulse,hx_01,hxr_10_mat(kr,:),hx_impulse,hy_01,hyr_10_mat(kr,:),hy_impulse,ex_01,ex_10,ex_impulse,ey_01,ey_10,ey_impulse] =...
            Calculate_Horizontal_Finite_Electrical_Source_unsave(I,L,h,x,y,z,t,fine);      
end 
%% reference response
% D model
save('parameters.txt','parametersD','-ascii');
[hz_01,hzD_10,hz_impulse,hx_01,hxD_10,hx_impulse,hy_01,hyD_10,hy_impulse,ex_01,ex_10,ex_impulse,ey_01,ey_10,ey_impulse] =...
            Calculate_Horizontal_Finite_Electrical_Source_unsave(I,L,h,x,y,z,t,fine);    
%% Calculate the rms 
% the sensitivity to the thickness
rms_hzh = rms((hzh_10_mat - repmat(hzD_10,hlen,1))'./repmat(hzD_10,hlen,1)');
rms_hyh = rms((hyh_10_mat - repmat(hyD_10,hlen,1))'./repmat(hyD_10,hlen,1)');
rms_hxh = rms((hxh_10_mat - repmat(hxD_10,hlen,1))'./repmat(hxD_10,hlen,1)');
figure;
loglog(rate_h,rms_hzh,'r','linewidth',2);
hold on;
loglog(rate_h,rms_hyh,'b','linewidth',2');
hold on;
loglog(rate_h,rms_hxh,'k','linewidth',2);
grid on;
legend('hz','hy','hx');
title('sensitivity to the thickness of high resistivity layer');
xlabel('h_a/h_1');
ylabel('RMS Difference');
%%
% the sensitivity to the resistivity
rms_hzr = rms((hzr_10_mat - repmat(hzD_10,rlen,1))'./repmat(hzD_10,rlen,1)');
rms_hyr = rms((hyr_10_mat - repmat(hyD_10,rlen,1))'./repmat(hyD_10,rlen,1)');
rms_hxr = rms((hxr_10_mat - repmat(hxD_10,rlen,1))'./repmat(hxD_10,rlen,1)');
%
figure;
loglog(rate_r,rms_hzr,'r','linewidth',2);
hold on;
loglog(rate_r,rms_hyr,'b','linewidth',2);
hold on;
loglog(rate_r,rms_hxr,'k','linewidth',2);
legend('hz','hy','hx');
grid on;
title('sensitivity to the resistivity ');
xlabel('rho_a/rho_1');
ylabel('RMS Difference');
%% recovery the parameters
save('parameters.txt','parameters0','-ascii');
load parameters.txt;


