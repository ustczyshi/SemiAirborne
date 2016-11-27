%% 利用仿真得到的发射电流数据模拟其激发的均匀大地瞬变电磁场，
% 以期研究关断局部异常对感应电压信号的影响
% by shizongyang
% 2015.11.24
%% 
clc;
close all;
clear all;
%% load the current data
% load Variable_step_discrete.mat;
load ATEM_current.mat;
fs = 1e7;
% fs = Fs;
current = atem_current(1:0.02*fs)';
t = (1:length(current))/fs;

%% the transmitting current
figure(1);
plot(t,current,'r');
title('current');
xlabel('time/s');
ylabel('I/A');

%%  利用垂直磁偶极子在均匀半空间大地表面产生的垂直磁场的计算式仿真
% 理想阶跃响应
% 参数设定
u0 = 4*pi*1e-7;%H/m
rou =1;% 电阻率
I = 1;
r = 10;
St = pi*r^2;% 发射线圈单匝面积m^2
Nt = 5;% 发射线圈匝数
Qt = St*Nt;% 发射线圈的有效面积
Sr = 1;% 接受线圈单匝面积m^2
Nr = 1;% 接受线圈匝数
% Srq = Sr*Nr;
Srq = 50;
%% 设定发射线圈的中心为坐标原点，垂直向下为+z
xr = 100;yr =0;zr = 0;% 接受线圈中心位置坐标m；
%%垂直磁偶极源下的均匀大地阶跃响应Hz Eo---
m = I*St;% 发射线圈的磁矩
t_ob = 1/fs:1/fs:0.1;%2秒的观测期
% [Hresponse01,Hresponse10,Himpulse10,Eresponse] = Vertical_Magnetic_Dipole_Step_Response(u0,rou,m,xr,yr,zr,t_ob);
[Hresponse01,Hresponse10,Himpulse10] = Central_Circle_Loop_Step_Response(u0,rou,I,r,t_ob);
%% 保存脉冲和阶跃响应
% himpulse10 = Himpulse10;
% hstep10 = Hresponse10;
% save('surface_magnetic_dipole_response_jx.mat','himpulse10','hstep10');

%% 考察阶跃响应和脉冲响应
% figure;
% plot(t,Hresponse10(1:length(t)),'r');%负阶跃响应
% xlabel('time/s');
% ylabel('Hz');
% title('负阶跃响应');
% figure;
% plot(t(1:end),Himpulse10(1:length(t)),'b');
% hold on;
% title('负脉冲响应');
% xlabel('time/s');
% ylabel('Hz');
%%
% figure;
% plot(t_ob*1000,Hresponse10,'r');%负阶跃响应
% hold on;
% title('负阶跃响应');
% xlabel('time/ms');
% ylabel('Hz');
% figure;
% plot(t_ob(1:end-1)*1000,(diff(Hresponse10)*fs),'k');%负阶跃响应
% hold on;
% legend('负阶跃响应的导数');
% xlabel('time/ms');
% ylabel('Hz');
% figure;
% plot(t_ob*1000,(Himpulse10),'b');
% hold on;
% legend('负脉冲响应');
% xlabel('time/ms');
% ylabel('Hz');
%% 考察阶跃响应和脉冲响应 对数域
% figure;
% loglog(t_ob*1000,Hresponse10,'r');%负阶跃响应
% hold on;
% loglog(t_ob(1:end-1)*1000,abs(diff(Hresponse10)*fs),'k');%负阶跃响应
% hold on;
% % loglog(t_ob(1:end-1)*1000,abs(diff([Hresponse01])*fs),'m');%负阶跃响应
% % hold on;
% loglog(t_ob*1000,abs(Himpulse10),'b');
% hold on;
% legend('负阶跃响应','负阶跃响应的导数','负脉冲响应');
% xlabel('time/ms');
% ylabel('Hz');
%
%% 求电流的导数
dIdt = Nt*diff([current 0]);%length(t)
response = conv(-dIdt,u0*Hresponse10);% 3*length(t)-1
Response = response(1:length(t_ob)); % 二次场
%% 直接利用脉冲响应的表达式计算有问题，数值偏大--------？？？
% response1 = conv(current,-u0*Himpulse10)/fs;% 3*length(t)-1
% Response1 = response1(1:length(t_ob));
%%
response1 = conv(dIdt,u0*Hresponse01);% 3*length(t)-1
Response1 = response1(1:length(t_ob)); %　总场
%%
figure;
plot((1:length(Response))/fs,Response,'r');% 负阶跃
hold on;
plot((1:length(Response1))/fs,Response1,'b');
legend('conv(dIdt,Bz\_step)','conv(I,Bz\_impulse)');
%%
figure;
[hx,h1,h2] = plotyy((1:length(Response))/fs,Response,(1:length(Response1))/fs,Response1);
grid on;
title(['ρ=' num2str(rou) 'Ω·m']);
legend('conv(di,Bstep-):Bz2','conv(di, hstep+):Bz');
xlabel('time/s');
ylabel(hx(1),'Bz2 /T','Color','b');% left y-axis
ylabel(hx(2),'Bz/T','Color','r');% right y-axis
h1.Color = 'b';
h2.Color = 'r';
%
%%
% Vinducted = -Srq*diff([Response Response(end) ])*fs;%3*length(t)-1
%%
Vinducted = -Srq*diff([Response1 Response1(end) ])*fs;%3*length(t)-1
%%
figure;
plot(t,current,'r');
hold on;
title('I');
xlabel('time/s');
ylabel('I/A');
%%
figure;
plot(t,Response1(1:length(t)),'r');
hold on;
plot(t,Nt*current.*(u0*Hresponse01(end)),'b');
hold on;
plot(t,Response1(1:length(t))-u0*Nt*current.*(Hresponse01(end)),'k');
hold on;
grid on;
legend('总场Bz','一次场Bz1','二次场Bz2');
xlabel('time/s');
ylabel('Bz/T');
%%
figure;
plot(t,Vinducted(1:length(t)),'r');
grid on;
title('Uz');
xlabel('time/s');
ylabel('U/V');
%
%% 
figure;
plot(t,(Vinducted(1:length(t))),'r');
grid on;
title('the induced voltage ');
xlabel('time/s');
ylabel('U/V');
%%
% figure;
% loglog(t,abs(Vinducted(1:length(t))),'b');
% legend('the induced voltage ');
% xlabel('time/s');
% ylabel('U/V');
%%
figure;
plot(t,(current(1:length(t)))./max(current(1:length(t))),'b');
hold on;
plot(t,Vinducted(1:length(t))./max(abs(Vinducted)),'r');
hold on;
grid on;
title(['ρ=' num2str(rou) 'Ω·m']);
legend('the normalized current','the normalized-induced voltage ');
xlabel('time/s');
ylabel('Uz/max(|Uz|)--I/max(|I|)');
%% 双坐标轴
figure;
[hx,h1,h2] = plotyy(t(1:0.01*fs),current(1:0.01*fs),t(1:0.01*fs),1000*(Vinducted(1:0.01*fs)));
grid on;
title(['ρ=' num2str(rou) 'Ω·m']);
legend('current','the induced voltage ');
xlabel('time/s');
ylabel(hx(1),'the current of  loop /A','Color','b');% left y-axis
ylabel(hx(2),'Uz /mV','Color','r');% right y-axis
h1.Color = 'b';
h2.Color = 'r';
%%
%{
% current1 = current(1:0.01*fs);
% inducted_voltage1 = Vinducted(1:0.01*fs);
% %%
% save('current1','current1');
% save('inducted_voltage1','inducted_voltage1');
%%
% figure;
% loglog(t(1:0.01*fs),abs(Vinducted(1:0.01*fs)),'b');
% % hold on;
% % semilogy(t(1:100000),abs(Vinducted(0.01*fs:0.01*fs+100000-1)),'r');
% legend('the induced voltage ');
% xlabel('time/s');
%%
% figure;
% loglog(t(1e-4*fs:0.01*fs),1000*abs(Vinducted(1e-4*fs:0.01*fs)),'b');%%认为100us关断时间
% legend('the induced voltage ');
% xlabel('time/s');
% ylabel('U/mV');
%}


