%--------------------------------------------------------------------------
%仿真条件：仿真垂直磁偶极的水平磁场随频率变化规律。
%磁偶极置于电阻率100 Ohm*m 的均匀大地表面，观察点距离偶极100m。
%--------------------------------------------------------------------------
format long;
mu_0 = 4*pi*1e-7;
freq= logspace(-1,5,250);
r = 100;
%--------------------------------------------------------------------------
%第一步：读取已经存储的滤波器系数,表示为行向量；
load FHT_J1_filter_Kong.txt;       
J_one = FHT_J1_filter_Kong( :, 3)';
delta = FHT_J1_filter_Kong( :, 2)';
%--------------------------------------------------------------------------
%计算lambda，并将lambda和frequency扩展成二维矩阵
lambda=(1./r) .*exp(delta);
[lambda_Array,frequency_Array] = Array_trans(lambda,freq);
%--------------------------------------------------------------------------
% %根据递推公式计算第一层的反射系数
r_TE = calculate_r_TE(lambda,freq);
%--------------------------------------------------------------------------
%计算磁场的水平分量，并用快速汉克尔变换求积分的数值解。  
sum = (1-r_TE).*lambda_Array.^2 ;
H_horizontal = 1/(4*pi) *  Fast_Hankel(r,sum,J_one);
%--------------------------------------------------------------------------
%作图。红色表示实分量，蓝色表示虚分量。虚线为负值，实线为正值。
loglog(freq,0.5*(abs(real(H_horizontal))-real(H_horizontal)),'r--','Linewidth',2)
hold on
loglog(freq,(real(H_horizontal)),'r','Linewidth',2)
hold on
loglog(freq,0.5*(abs(imag(H_horizontal))-imag(H_horizontal)),'b--','Linewidth',2)
hold on 
loglog(freq,(imag(H_horizontal)),'b','Linewidth',2)
xlabel('频率（Hz）')
ylabel('Hz(A/m)')
title('垂直磁偶极的水平磁场随频率的变化曲线')
axis([10^-1,10^5,10^-12,10^-6])