%--------------------------------------------------------------------------
%��������������TEM �뺽��TEM
% ����ˮƽ�ӵس�����Դ�������������й۲� B��
% 1. ���Ȱ�ռ�������趨100ŷķ��,�ȹ۲ⲻͬ�߶���Ӧ�ı仯����Ҫ�ǹ۲�Hz,Hx,Hy,�߶��趨0,30,60,90��
% 
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
x = 100; % �շ�ˮƽƫ�ƾ࣬��y��
y = 1000;
L= 2500; % �������³��ȣ���x��
I = 40; % �������
m = I*L;
%%  �뺽���շ��߶Ȳ���
z =0;% �۲������ĸ߶ȣ���������Ϊ��ֵ
h =0;% Դ�����ĸ߶�

%% �����ʺ͹۲�ʱ�������
fs = 1e7;% ������
dt = 1./fs;
% t = 1/fs:1/fs:2e-3;% ʱ������
 t = logspace(-5,0,100);
%%

%% ������
% ���޳��ӵص����д��ߵ����ϵĽ�����
% [ hz_impulse_jiexijie] = ground_finite_wire_source_jiexi(u0,rou,I,L,y,t);
% save(['ground_finite_wire_source_jiexi_y' num2str(y) '.mat'],'hz_impulse_jiexijie');
% ˮƽ��żԴ
[step_ex01,step_ey01,impulse_hx,impulse_hy,impulse_hz] = Horizontal_Electrical_Dipole_Response(u0,rou,m,x,y,z,t);
% fine = 0.5;
% [hz_01,hz_10,hz_1_impulse,hx_01,hx_10,hx_1_impulse,hy_01,hy_10,hy_1_impulse,ex_01,ex_10,ex_impulse,ey_01,ey_10,ey_impulse] =...
%             Calculate_Horizontal_Finite_Electrical_Source(I,L,h,x,y,z,t,fine);
%% ���ø�˹-���õ»���ʵ���س����ߵĻ�������
n = 10;
[hz_01,hz_10,hz_1_impulse,hx_01,hx_10,hx_1_impulse,hy_01,hy_10,hy_1_impulse,ex_01,ex_10,ex_impulse,ey_01,ey_10,ey_impulse] =...
            Calculate_Horizontal_Finite_Electrical_Source_GuaLeg(I,L,h,x,y,z,t,n);
%         [ hz_impulse_jiexijie] = ground_finite_wire_source_jiexi(u0,rou,I,L,y,t);
%% Ex-Ey 
%
% ���ӵس����߳ߴ�ԴС��ƫ�ƾ��ǣ�������˲����������ˮƽ��żԴ��Ӧ�ĶԱ�
% Ex  step response
figure;
loglog(t,abs(ex_01) ,'r','linewidth',2);
hold on;
loglog(t,abs(step_ex01) ,'b','linewidth',1);
hold on;
grid on;
legend('��ֵ��ex\_1\_step','������step\_ex01');
title(['source moment' num2str(I) 'A*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop  Ex step response'])
xlabel('Time/(ms)')
ylabel('E/(V/m)');
% Ey impulse response
figure;
loglog(t,abs(ey_01) ,'r','linewidth',2);
hold on;
loglog(t,abs(step_ey01) ,'b','linewidth',1);
hold on;
grid on;
legend('��ֵ��ey_01','������step_ey01');
title(['source moment' num2str(I) 'Ax' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop Ey step response'])
xlabel('Time/(ms)')
ylabel('E/(V/m)');
%
% Hz impulse
figure;
loglog(t(1:end).*10^3,u0.*abs(hz_1_impulse),'r','Linewidth',2);
hold on
loglog(t(1:end).*10^3,u0.*abs(impulse_hz),'b','Linewidth',1);
grid on;
legend('��ֵ��hz\_1\_impulse','������impulse\_hz');
title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop Bz impulse response'])
xlabel('Time/(ms)')
ylabel('Bz/(T)');
% Hx impulse
figure;
loglog(t(1:end).*10^3,u0.*abs(hx_1_impulse),'r','Linewidth',2);
hold on
loglog(t(1:end).*10^3,u0.*abs(impulse_hx),'b','Linewidth',1);
grid on;
legend('��ֵ��hx\_1\_impulse','������impluse\_hx');
title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop Bx step response'])
xlabel('Time/(ms)')
ylabel('Bx/(T)');
% Hy impulse
figure;
loglog(t(1:end).*10^3,u0.*abs(hy_1_impulse),'r','Linewidth',2);
hold on
loglog(t(1:end).*10^3,u0.*abs(impulse_hy),'b','Linewidth',1);
grid on;
legend('��ֵ��hy\_1\_impulse','������impulse\_hy');
title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop By step response'])
xlabel('Time/(ms)')
ylabel('By/(T)');
%
%%
% hz impulse response
% figure;
% loglog(t(1:end).*10^3,u0.*abs(hz_1_impulse),'r','Linewidth',2);
% hold on
% loglog(t(1:end-1).*10^3,u0.*abs(diff(hz_01)./dt),'b','Linewidth',2);
% grid on;
% legend('��ֵ��hz\_1\_impulse','��ֵ��diff(hz\_01)./dt');
% title(['source moment' num2str(I) 'A*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop Bz impulse response'])
% xlabel('Time/(ms)')
% ylabel('Bz/(T)');
%{
%% �ӵس����ߵ�hz��ֵ��ͽ�����Ա�
figure;
loglog(t(1:end).*10^3,u0.*abs(hz_1_impulse),'r','Linewidth',2);
hold on
loglog(t(1:end).*10^3,u0.*abs(hz_impulse_jiexijie),'b','Linewidth',2);
grid on;
legend('��ֵ��hz\_1\_impulse','������hz\_impulse\_jiexijie');
title(['source moment' num2str(I) 'A*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop Bz impulse response'])
xlabel('Time/(ms)')
ylabel('Bz/(T)');
%%
% ��֤�ӵس����ߵ�Bz��Ӧ
figure;
plot(t(1:end).*1e3,abs(hz_impulse_jiexijie)./abs(hz_1_impulse));
title('Bz��ֵ��ͽ���������')
xlabel('Time/(ms)')
ylabel('error/(T)');

figure;
plot(t(1:end).*1e3,(abs(hz_impulse_jiexijie)./abs(hz_1_impulse)-1).*100);
title('Bz��ֵ��ͽ���������')
xlabel('Time/(ms)')
ylabel('error/(%)');
%}



