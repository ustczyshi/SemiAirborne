%% 通过高斯-勒让德积分原理，先离散化长接地导线源
% ================================
% 调用子函数  GuaLegendIntegral_Coef.m
% 先得到勒让德各分节点的坐标值xk
% 及对应的积分系数Ak
% 接地长导线源位于地面沿x轴

function [Ak,xk,dxk,r] = GuaLeg_DiscreteSource(obs_point,L,n,tol)
% ================================
%                           input parameters
% obs_point: 向量，[x,y,z],存放观测点坐标值；
% L: 长接地导线长度；
% n: GuaLegend 高斯-勒让德积分节点数，一般取7~12；
%     节点数过多不一定积分精度高
% -----------------------------------------------------------
%                          output parameters
% Ak : 积分系数，列向量；
% xk : 积分节点，实际是n阶勒让德多项式的n各节点，为列向量；
% dxk : 列向量，存放分成的长导线各段的长度；
% r : 各分点距观测点的水平距离,列向量；
% ================================
switch nargin 
    case 4
        disp('------------输入数据正确！----------');
    case 3
        tol = 1e-10;
    case 2
        n = 8;% 默认 8各分点
    otherwise
        error('--!!! 输入错误，请重新输入参数!!!--');
end
% 求的n阶勒让德多项式的各根，作为分隔点（积分节点）
[C,ZP]=Legendre_Roots(n); % ZP n阶勒让德多项式的零点xk
% 计算各零点（积分节点）对应的积分系数Ak
[Ak ,roots]=GuaLegendIntegral_Coef(ZP,tol);
a = -L/2;
b = L/2;
% 积分变量代换，将[a,b]变换到[-1,1]
xk=(b-a)/2*roots+(b+a)/2;% 列向量
% 长导线各段的长度；
dxk = [ b xk']'-[xk' a]';
% 各分点距观测点坐标
x = obs_point(1);% 第一列为x轴坐标，列向量
y = obs_point(2);% 第二列为y轴坐标，列向量
z = obs_point(3);% 第三列为z轴坐标，列向量
% 各分点距观测点的水平距离
r = sqrt((x-xk).^2+y.^2);
r1 = ((x+L/2)^2+y^2)^(0.5);
r2 = ((x-L/2)^2+y^2)^(0.5);
end