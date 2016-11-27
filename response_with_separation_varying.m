%% 比较相同偏移距下，不同接收高度响应的衰减情况
% clc;close all;
% %
% load parameters.txt;
% sigma1 = parameters(1,2);%第一层的电导率
% rou = 1./sigma1;
% u0 = 4*pi*1e-7;
% fs = 1e5;
% dt = 1/fs;
% x = 500;
% L = 2000;
% I=40;
% h = 0;
% y = [100,500,1000,2000,10000];
% % y = 10;
% z =-30;
% N = 4e-2*fs;
% % t = (1:N).*dt;
% clear t;
% t = logspace(-5,0,100);% 时间区间
N = length(t);
%% dataname


hz_step_mat = zeros(length(y),N);
hz_impulse_mat = zeros(length(y),N);
hz_impulse_jiexi = zeros(length(y),N);
hy_step_mat = zeros(length(y),N);
hy_impulse_mat = zeros(length(y),N);

hx_step_mat = zeros(length(y),N);
hx_impulse_mat = zeros(length(y),N);
for kk = 1:length(y)
    dataname0 = ['ground_finite_wire_source_jiexi_rou' num2str(rou) '_L' num2str(L) '_y' num2str(y(kk)) '.mat'];
    load(dataname0);
    hz_impulse_jiexi(kk,:) = -hz_impulse_jiexijie;
%     for kk = 1:length(z)
%         dataname =['SemiAtem_Horizontal_Finite_Electrical_Source_L' num2str(L) '_h' num2str(h) '_z' num2str(z(kk)) '_x' num2str(x) '_y' num2str(y(k)) '.mat'];
        dataname = ['SemiAtem_Horizontal_Finite_Electrical_Source_separation_varying_L'...
            num2str(L) '_h' num2str(h) '_z' num2str(z) '_x' num2str(x) '_y' num2str(y(kk)) '.mat'];
%         dataname = ['SemiAtem_Horizontal_Electrical_Dipole_separation_varying_L'...
%             num2str(L) '_h' num2str(h) '_z' num2str(z) '_x' num2str(x) '_y' num2str(y(kk)) '.mat'];
        load(dataname);
        hz_step_mat(kk,:) = hz_10;
        hz_impulse_mat(kk,:) = hz_impulse;
        hy_step_mat(kk,:) = hy_10;
        hy_impulse_mat(kk,:) = hy_impulse;
        hx_step_mat(kk,:) = hx_10;
        hx_impulse_mat(kk,:) = hx_impulse;
%     end
end
%% display;
% Hz
figure;
loglog(t,hz_step_mat,'Linewidth',2);
title('Hz');
% title(['Grounded Finite wire L=' num2str(L) ' offset=' num2str(y) ' rou =' num2str(rou)]);
legend([  num2str(y(1)) 'm'],[ num2str(y(2)) 'm'],...
    [ num2str(y(3)) 'm'],[ num2str(y(4)) 'm'],[ num2str(y(5)) 'm'],[ num2str(y(6)) 'm']);
grid on;
xlabel('Time/s');
ylabel('Hz/(A/m)');
%%  Uz
figure;
loglog(t(1:end),u0.*hz_impulse_mat,'Linewidth',2);%./repmat(diff(t),length(y),1)
hold on;
% loglog(t(1:end),u0.*hz_impulse_jiexi,'--','Linewidth',2);
% hold on;
title('dBz/dt');
% title(['Grounded Finite wire L=' num2str(L) ' offset=' num2str(y) ' rou =' num2str(rou)]);
legend([  num2str(y(1)) 'm'],[ num2str(y(2)) 'm'],...
    [ num2str(y(3)) 'm'],[ num2str(y(4)) 'm'],[ num2str(y(5)) 'm'],[ num2str(y(6)) 'm']);
grid on;
xlabel('Time/s');
ylabel('dBz/dt(V/m^2)');

%%
figure;
% loglog(t(1:end),u0.*hz_impulse,'Linewidth',2);%./repmat(diff(t),length(y),1)
% hold on;
loglog(t(1:end),u0.*hz_impulse_mat,'Linewidth',2);
hold on;
title('Uz with different offset');
% title(['Grounded Finite wire L=' num2str(L) ' offset=' num2str(y) ' rou =' num2str(rou)]);
legend([  num2str(y(1)) 'm'],[ num2str(y(2)) 'm'],...
    [ num2str(y(3)) 'm'],[ num2str(y(4)) 'm'],[ num2str(y(5)) 'm'],[ num2str(y(1)) 'm']);
grid on;
xlabel('Time/s');
ylabel('dBz/dt/(V/m^2)');

%
%%
% Hy
figure;
loglog(t,abs(hy_step_mat),'Linewidth',2);
title('Hy');
legend([  num2str(y(1)) 'm'],[  num2str(y(2)) 'm'],...
    [ num2str(y(3)) 'm'],[ num2str(y(4)) 'm'],[  num2str(y(5)) 'm'],[  num2str(y(6)) 'm']);
grid on;
xlabel('Time/s');
ylabel('Hy/(A/m)');

figure;
loglog(t,u0.*abs(hy_impulse_mat),'Linewidth',2);
% title('脉冲响应hy随接收高度的变化');
title('dBy/dt');
legend([ num2str(y(1)) 'm'],[ num2str(y(2)) 'm'],...
    [  num2str(y(3)) 'm'],[  num2str(y(4)) 'm'],[ num2str(y(5)) 'm'],[ num2str(y(6)) 'm']);
grid on;
xlabel('Time/s');
ylabel('dBy/dt/(V/m^2)');
%%
% Hx
figure;
loglog(t,abs(hx_step_mat),'Linewidth',2);
title('Hx with different offset');
legend([  num2str(y(1)) 'm'],[  num2str(y(2)) 'm'],...
    [ num2str(y(3)) 'm'],[ num2str(y(4)) 'm'],[  num2str(y(5)) 'm'],[ num2str(y(6)) 'm']);
grid on;
xlabel('Time/s');
ylabel('Hx/(A/m)');
figure;
loglog(t,u0.*abs(hx_impulse_mat),'Linewidth',2);
% title('脉冲响应hy随接收高度的变化');
title('Ux with different offset');
legend([ num2str(y(1)) 'm'],[ num2str(y(2)) 'm'],...
    [  num2str(y(3)) 'm'],[  num2str(y(4)) 'm'],[ num2str(y(5)) 'm'],[ num2str(y(6)) 'm']);
grid on;
xlabel('Time/s');
ylabel('dBx/dt/(V/m^2)');
%%
% Uz-Uy-Ux
figure;
loglog(t,hz_impulse_mat(4,:),'r','Linewidth',2);
hold on;
loglog(t,abs(hy_impulse_mat(4,:)),'b','Linewidth',2);
hold on;
loglog(t,abs(hx_impulse_mat(4,:)),'k','Linewidth',2);
hold on;
title(['(dBx/dt,dBy/dt,dBz/dt) at('  num2str(x)  ',' num2str(y(4))  ',' num2str(z) ')' ]);
% title(['Grounded Finite wire L=' num2str(L) ' offset=' num2str(y) ' rou =' num2str(rou)]);
legend( 'dBz/dt', 'dBy/dt', 'dBx/dt');
grid on;
xlabel('Time/s');
ylabel('dB/dt(V/m^2)');
% Hz-Hy-Hx
figure;
loglog(t,hz_step_mat(4,:),'r','Linewidth',2);
hold on;
loglog(t,abs(hy_step_mat(4,:)),'b','Linewidth',2);
hold on;
loglog(t,abs(hx_step_mat(4,:)),'k','Linewidth',2);
hold on;
title(['(hx,hy,hz) at('  num2str(x)  ',' num2str(y(4))  ',' num2str(z) ')' ]);
% title(['Grounded Finite wire L=' num2str(L) ' offset=' num2str(y) ' rou =' num2str(rou)]);
legend( 'hz', 'hy', 'hx');
grid on;
xlabel('Time/s');
ylabel('H/(A/m)');

%{
%% Ex
figure;
loglog(t,abs(ex_step_y100));
title('正阶跃响应ex随接收高度的变化');
legend([ '接收高度' num2str(z(1)) 'm'],[ '接收高度' num2str(z(2)) 'm'],...
    [ '接收高度' num2str(z(3)) 'm'],[ '接收高度' num2str(z(4)) 'm'],[ '接收高度' num2str(z(5)) 'm']);
grid on;
xlabel('Time/s');
ylabel('ex/(V/m)');
figure;
loglog(t,abs(ex_impulse_y100));
title('脉冲响应ex随接收高度的变化');
legend([ '接收高度' num2str(z(1)) 'm'],[ '接收高度' num2str(z(2)) 'm'],...
    [ '接收高度' num2str(z(3)) 'm'],[ '接收高度' num2str(z(4)) 'm'],[ '接收高度' num2str(z(5)) 'm']);
grid on;
xlabel('Time/s');
ylabel('ex/(V/m)');
%}