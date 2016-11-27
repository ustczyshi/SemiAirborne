%% 比较相同偏移距下，不同接收高度响应的衰减情况
% clc;close all;
% %
% u0 = 4*pi*1e-7;
% fs = 1e5;
% dt = 1/fs;
% L = 2000;
% I = 40;
% h = 0;
% load parameters.txt;
% sigma1 = parameters(1,2);%第一层的电导率
% rou = 1./sigma1;
% x = 500; % 收发水平偏移距，沿y轴
% y = 500;
% z = [0,-30,-60,-90,-150];
% % N = 4e-2*fs;
% % t = (1:N).*dt;
% t = logspace(-5,1,100);% 时间区间
N = length(t);
%% dataname
hz_step_y100 = zeros(length(z),N);
hz_impulse_y100 = zeros(length(z)+1,N);

hy_step_y100 = zeros(length(z),N);
hy_impulse_y100 = zeros(length(z),N);

hx_step_y100 = zeros(length(z),N);
hx_impulse_y100 = zeros(length(z),N);

ex_step_y100 = zeros(length(z),N);
ex_impulse_y100 = zeros(length(z),N);

ey_step_y100 = zeros(length(z),N);
ey_impulse_y100 = zeros(length(z),N);
for k = 1:length(y)
%     dataname0 = ['ground_finite_wire_source_jiexi_rou' num2str(rou) '_L' num2str(L) '_y' num2str(y(k)) '.mat'];
%     load(dataname0);
%     hz_impulse_y100(length(z)+1,:) = -hz_impulse_jiexijie;
    for kk = 1:length(z)
        dataname =['SemiAtem_Horizontal_Finite_Electrical_Source_height_varying_L' num2str(L) '_h' num2str(h) '_z' num2str(z(kk)) '_x' num2str(x) '_y' num2str(y(k)) '.mat'];
        load(dataname);
        hz_step_y100(kk,:) = hz_10;
        hz_impulse_y100(kk,:) = hz_impulse;
        hx_step_y100(kk,:) = hx_10;
        hx_impulse_y100(kk,:) = hx_impulse;
        hy_step_y100(kk,:) = hy_10;
        hy_impulse_y100(kk,:) = hy_impulse;
        ex_step_y100(kk,:) = ex_10;
        ex_impulse_y100(kk,:) = ex_impulse;
        ey_step_y100(kk,:) = ey_10;
        ey_impulse_y100(kk,:) = ey_impulse;
    end
end
%% display 
%   H场
figure;
loglog(t,abs(hz_step_y100(3,:) ),'r','linewidth',2);
hold on;
loglog(t,abs( hx_step_y100(3,:) ),'b','linewidth',2);
hold on;
loglog(t,abs( hy_step_y100(3,:) ),'k','linewidth',2);
hold on;
title(['L=' num2str(L) 'm, Magnetic@ (' num2str(x) 'm,' num2str(y) 'm,' num2str(z(3)) 'm)' ]);
legend('Hz','Hx','Hy');
legend boxoff;
grid on;
xlabel('Time/s');
ylabel('Magnetic/(A/m)');
% dB场
figure;
loglog(t,u0.*abs(hz_impulse_y100(3,:) ),'r','linewidth',2);
hold on;
loglog(t,u0.*abs( hx_impulse_y100(3,:) ),'b','linewidth',2);
hold on;
loglog(t,u0.*abs( hy_impulse_y100(3,:) ),'k','linewidth',2);
hold on;
title(['L=' num2str(L) 'm, dB/dt@ (' num2str(x) 'm,' num2str(y) 'm,' num2str(z(3)) 'm)' ]);
legend('dBz/dt','dBx/dt','dBy/dt');
legend boxoff;
grid on;
xlabel('Time/s');
ylabel('dB/dt(T/s^2)');
%
%% Hz
figure;
loglog(t,hz_step_y100,'linewidth',2);
% title('负阶跃响应hz随接收高度的变化');
title(['L=' num2str(L) 'm, dB/dt@ (' num2str(x) 'm,' num2str(y) 'm,' num2str(z(3)) 'm)' ]);
% title(['Grounded Finite wire L=' num2str(L) ' offset=' num2str(y) ' rou =' num2str(rou)]);
legend([ 'z=' num2str(z(1)) 'm'],[ 'z=' num2str(z(2)) 'm'],...
    [ 'z=' num2str(z(3)) 'm'],[ 'z=' num2str(z(4)) 'm'],[ 'z=' num2str(z(5)) 'm']);
legend boxoff;
grid on;
xlabel('Time/s');
ylabel('Hz/(A/m)');
% percentage difference
error_hz =(1- abs( (hz_step_y100(1:length(z),:))./repmat(hz_step_y100(1,:),length(z),1) ));
figure;
semilogx(t,error_hz.*100,'linewidth',2);
title(['L=' num2str(L) 'm, dB/dt@ (' num2str(x) 'm,' num2str(y) 'm,' num2str(z(3)) 'm)' ]);
% title(['Grounded Finite wire L=' num2str(L) ' offset=' num2str(y) ' rou =' num2str(rou)]);
legend([ 'z=' num2str(z(1)) 'm'],[ 'z=' num2str(z(2)) 'm'],...
    [ 'z=' num2str(z(3)) 'm'],[ 'z=' num2str(z(4)) 'm'],[ 'z=' num2str(z(5)) 'm']);
legend boxoff;
grid on;
xlabel('Time/s');
ylabel('Hz percentage difference/(%)');
%% 计算最大误差
error_hz_mat  = max(error_hz');
figure;
plot(abs(z),error_hz_mat);
hold on;
grid on;
title('负阶跃响应误差随高度的变化');
xlabel('heght/m');
ylabel('percentage difference /(%)');
%% 脉冲响应
figure;
loglog(t,u0.*(hz_impulse_y100),'linewidth',2);
% title('脉冲响应hz随接收高度的变化');
title(['L=' num2str(L) 'm, dB/dt@ (' num2str(x) 'm,' num2str(y) 'm,' num2str(z(3)) 'm)' ]);
% title(['Grounded Finite wire L=' num2str(L) ' offset=' num2str(y) ' rou =' num2str(rou)]);
legend([ 'z=' num2str(z(1)) 'm'],[ 'z=' num2str(z(2)) 'm'],...
    [ 'z=' num2str(z(3)) 'm'],[ 'z=' num2str(z(4)) 'm'],[ 'z=' num2str(z(5)) 'm']);
legend boxoff;
grid on;
xlabel('Time/s');
ylabel('dBz/dt/(V/m^2)');
% percentage difference
error_vz =(1-abs( (hz_impulse_y100(1:length(z),:))./repmat(hz_impulse_y100(1,:),length(z),1) ));
figure;
semilogx(t,error_vz.*100,'linewidth',2);
title(['L=' num2str(L) 'm, dB/dt@ (' num2str(x) 'm,' num2str(y) 'm,' num2str(z(3)) 'm)' ]);
% title(['Grounded Finite wire L=' num2str(L) ' offset=' num2str(y) ' rou =' num2str(rou)]);
legend([ 'z=' num2str(z(1)) 'm'],[ 'z=' num2str(z(2)) 'm'],...
    [ 'z=' num2str(z(3)) 'm'],[ 'z=' num2str(z(4)) 'm'],[ 'z=' num2str(z(5)) 'm']);
legend boxoff;
grid on;
xlabel('Time/s');
ylabel('Uz percentage difference/(%)');
%% 计算最大误差
error_mat  = max(abs(error_vz'));
figure;
plot(abs(z),error_mat);
hold on;
title('响应误差随高度的变化');
xlabel('heght/m');
ylabel('percentage difference /(%)');
%%
% Hx
figure;
loglog(t,abs(hx_step_y100),'linewidth',2);
title('hx with different height');
legend([ 'z=' num2str(z(1)) 'm'],[ 'z=' num2str(z(2)) 'm'],...
    [ 'z=' num2str(z(3)) 'm'],[ 'z=' num2str(z(4)) 'm'],[ 'z=' num2str(z(5)) 'm']);
legend boxoff;
grid on;
xlabel('Time/s');
ylabel('Hx/(A/m)');
figure;
loglog(t,u0.*abs(hx_impulse_y100),'linewidth',2);
title('dBx/dt with different height');
legend([ 'z=' num2str(z(1)) 'm'],[ 'z=' num2str(z(2)) 'm'],...
    [ 'z=' num2str(z(3)) 'm'],[ 'z=' num2str(z(4)) 'm'],[ 'z=' num2str(z(5)) 'm']);
legend boxoff;
grid on;
xlabel('Time/s');
ylabel('dBx/dt/(V/m^2)');
%% Hy
figure;
loglog(t,abs( hy_step_y100),'linewidth',2);
% loglog(t,abs(hy_step_y100(end) - hy_step_y100)./L,'linewidth',2);
title('Hy with different height');
legend([ 'z=' num2str(z(1)) 'm'],[ 'z=' num2str(z(2)) 'm'],...
    [ 'z=' num2str(z(3)) 'm'],[ 'z=' num2str(z(4)) 'm'],[ 'z=' num2str(z(5)) 'm']);
legend boxoff;
grid on;
xlabel('Time/s');
ylabel('Hy/(A/m)');

% percentage difference
error_hy =(1-  abs( (hy_step_y100(1:length(z),:))./repmat(hy_step_y100(1,:),length(z),1) ));
figure;
semilogx(t,error_hy.*100,'linewidth',2);
title(['L=' num2str(L) 'm, dB/dt@ (' num2str(x) 'm,' num2str(y) 'm,' num2str(z(3)) 'm)' ]);
% title(['Grounded Finite wire L=' num2str(L) ' offset=' num2str(y) ' rou =' num2str(rou)]);
legend([ 'z=' num2str(z(1)) 'm'],[ 'z=' num2str(z(2)) 'm'],...
    [ 'z=' num2str(z(3)) 'm'],[ 'z=' num2str(z(4)) 'm'],[ 'z=' num2str(z(5)) 'm']);
legend boxoff;
grid on;
xlabel('Time/s');
ylabel('Hy percentage difference/(%)');
%%
figure;
loglog(t,u0.*abs(hy_impulse_y100),'linewidth',2);
title(['L=' num2str(L) 'm, dBy/dt@ (' num2str(x) 'm,' num2str(y) 'm,' num2str(z(3)) 'm)' ]);
legend([ 'z=' num2str(z(1)) 'm'],[ 'z=' num2str(z(2)) 'm'],...
    [ 'z=' num2str(z(3)) 'm'],[ 'z=' num2str(z(4)) 'm'],[ 'z=' num2str(z(5)) 'm']);
legend boxoff;
grid on;
xlabel('Time/s');
ylabel('dBy/dt/(V/m^2)');

%% percentage difference
error_vy =(1-  abs( (hy_impulse_y100(1:length(z),:))./repmat(hy_impulse_y100(1,:),length(z),1) ));
figure;
semilogx(t,error_vy.*100,'linewidth',2);
title(['Grounded Finite wire L=' num2str(L) ' offset=' num2str(y) ' rou =' num2str(rou)]);
legend([ 'z=' num2str(z(1)) 'm'],[ 'z=' num2str(z(2)) 'm'],...
    [ 'z=' num2str(z(3)) 'm'],[ 'z=' num2str(z(4)) 'm'],[ 'z=' num2str(z(5)) 'm']);
legend boxoff;
grid on;
xlabel('Time/s');
ylabel('Uy percentage difference/(%)');
%
%% Ex
% figure;
% loglog(t,abs( ex_step_y100),'linewidth',2);
% title('负阶跃响应ex随接收高度的变化');
% legend([ 'z=' num2str(z(1)) 'm'],[ 'z=' num2str(z(2)) 'm'],...
%     [ 'z=' num2str(z(3)) 'm'],[ 'z=' num2str(z(4)) 'm'],[ 'z=' num2str(z(5)) 'm']);
% grid on;
% xlabel('Time/s');
% ylabel('ex/(V/m)');
%% Ey
%
% figure;
% loglog(t,abs(ey_step_y100),'linewidth',2);
% title('负阶跃响应ey随接收高度的变化');
% legend([ 'z=' num2str(z(1)) 'm'],[ 'z=' num2str(z(2)) 'm'],...
%     [ 'z=' num2str(z(3)) 'm'],[ 'z=' num2str(z(4)) 'm'],[ 'z=' num2str(z(5)) 'm']);
% grid on;
% xlabel('Time/s');
% ylabel('ey/(V/m)');
%