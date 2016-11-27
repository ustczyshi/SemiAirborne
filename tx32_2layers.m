%%  ����TX32 TEM���Ļ��ߣ�����ģ�ͶԴ��ڹ���ĵ�����������Ӧ���յ�ѹ�ź�
%----------------------------       -------------------------------    
%   rou = 100  h = 20              rou = 0.1  h = 20 
%----------------------------       --------------------------------
%   rou = 0.1 h = 50                    rou = 100  h = 50
%----------------------------       --------------------------------
%   rou = 100  h = 0              rou = 100  h = 0
clc;close all;clear all;
%%
mu_0 = 4*pi*1e-7;
%% load the current data
% load Variable_step_discrete.mat;
% load ATEM_current.mat;% �й���
load TX41_with_overshoot.mat;% �޹���
fs = 1e7;
current = tx41_current_1(1:0.01*fs);% row vector
Ns = length(current);
t = (1:Ns)/fs;
% ��������
figure;
plot(t,current,'r');
title('the ATEM current')
xlabel('time/s');
ylabel('current/A');
%%
%  h = [50 100 150 200 250 300 350 400 450 500 ];
% h = [0 50 150 250 350 450 ];
h = [ 100 200 300 400 500];
rou = 0.001;%% ʵ��
r = 0;
L =30;
z = -(h - L);
% z(1) = 0;
color = ['r','b','g','k','c','m'];
tob = (1/fs:1/fs:1e-2);
Hzimpulse01 =  zeros(length(h),length(tob));
Hzstep01 =  zeros(length(h),length(tob));
Hrimpulse01 =  zeros(length(h),length(tob));
Hrstep01 =  zeros(length(h),length(tob));
epimpulse01 = zeros(length(h),length(tob));
epstep01 =  zeros(length(h),length(tob));

Hzimpulse10 =  zeros(length(h),length(tob));
Hzstep10 =  zeros(length(h),length(tob));
Hrimpulse10 =  zeros(length(h),length(tob));
Hrstep10 =  zeros(length(h),length(tob));
epimpulse10 = zeros(length(h),length(tob));
epstep10 =  zeros(length(h),length(tob));


%%
for k = 1:length(h)
    filename = ['atem_circle_loop_response_h' num2str(h(k)) 'z' num2str(z(k)) '_r' num2str(r) '_rou' num2str(rou) '.mat'];
    load(filename);
    Hzimpulse01(k,:) = hzimpulse01;
    Hzstep01(k,:) = hzstep01;
    Hrimpulse01(k,:) = hrimpulse01;
    Hrstep01(k,:) = hrstep01;
    epimpulse01(k,:) = epimpulse01;
    epstep01(k,:) = epstep01;
    
    Hzimpulse10(k,:) = hzimpulse10;
    Hzstep10(k,:) = hzstep10;
    Hrimpulse10(k,:) = hrimpulse10;
    Hrstep10(k,:) = hrstep10;
    epimpulse10(k,:) = epimpulse10;
    epstep10(k,:) = epstep10;
end
%%
%{
figure;
loglog(t*1e3,abs(Hzimpulse));
hold on;
grid on;
title('Hz-impulse');
legend(['h=' num2str(h(1)) 'm'],['h=' num2str(h(2)) 'm'],['h=' num2str(h(3)) 'm'],['h=' num2str(h(4)) 'm'],['h=' num2str(h(5)) 'm'],...
    ['h=' num2str(h(6)) 'm']);
xlabel('time/(ms)');
ylabel('Hz/(A/m)');
%
figure;
loglog(t*1e3,(Hzstep));
hold on;
grid on;
title('Hz-step');
legend(['h=' num2str(h(1)) 'm'],['h=' num2str(h(2)) 'm'],['h=' num2str(h(3)) 'm'],['h=' num2str(h(4)) 'm'],['h=' num2str(h(5)) 'm'],...
    ['h=' num2str(h(6)) 'm']);
xlabel('time/(ms)');
ylabel('Hz/(A/m)');
%%
figure;
loglog(t*1e3,abs(Hrimpulse));
hold on;
grid on;
title('Hr-impulse');
legend(['h=' num2str(h(1)) 'm'],['h=' num2str(h(2)) 'm'],['h=' num2str(h(3)) 'm'],['h=' num2str(h(4)) 'm'],['h=' num2str(h(5)) 'm'],...
    ['h=' num2str(h(6)) 'm']);
xlabel('time/(ms)');
ylabel('Hr/(A/m)');
%
figure;
loglog(t*1e3,abs(Hrstep));
hold on;
grid on;
title('Hr-step');
legend(['h=' num2str(h(1)) 'm'],['h=' num2str(h(2)) 'm'],['h=' num2str(h(3)) 'm'],['h=' num2str(h(4)) 'm'],['h=' num2str(h(5)) 'm'],...
    ['h=' num2str(h(6)) 'm']);
xlabel('time/(ms)');
ylabel('Hr/(A/m)');
%%
figure;
loglog(t*1e3,abs(epimpulse));
hold on;
grid on;
title('Ephi-impulse');
legend(['h=' num2str(h(1)) 'm'],['h=' num2str(h(2)) 'm'],['h=' num2str(h(3)) 'm'],['h=' num2str(h(4)) 'm'],['h=' num2str(h(5)) 'm'],...
    ['h=' num2str(h(6)) 'm']);
xlabel('time/(ms)');
ylabel('Ephi/(V/m)');
%
figure;
loglog(t*1e3,(epstep));
hold on;
grid on;
title('Ephi-step');
legend(['h=' num2str(h(1)) 'm'],['h=' num2str(h(2)) 'm'],['h=' num2str(h(3)) 'm'],['h=' num2str(h(4)) 'm'],['h=' num2str(h(5)) 'm'],...
    ['h=' num2str(h(6)) 'm']);
xlabel('time/(ms)');
ylabel('Ephi/(V/m)');
%}
%%
%%  ���Ʋ�ͬ�߶ȣ�������Ȧ����Ծ��Ӧ
%% didt
n = 5;%������Ȧ����
didt =n*diff([current 0 ]);
figure;
plot(t,didt*fs,'r');
title('the ATEM current  derivation')
xlabel('time/s');
ylabel('dIdt');
%% Bz ��о��ȦSz = 50;
Sz = 50;% ��о��Ȧ���
% Sz = 490;% ��о��Ȧ���
%% ************************���ø���Ծ��Ӧ����**************************************
Bz2 = zeros(length(h),length(t)+length(tob)-1);
for  k = 1:length(h)
    Bz2(k,:) = mu_0*conv(-didt,Hzstep10(k,:));
end
Bz1 = zeros(length(h),length(t));
for  k = 1:length(h)
    Bz1(k,:) = mu_0*Hzstep01(k,end)*n*current;
end
Bz = Bz1+Bz2(:,1:length(t));
Uz2 = -Sz.*diff(Bz2')'*fs;
Uz = -Sz.*diff(Bz')'*fs;
%% hr
Br2 = zeros(length(h),length(t)+length(tob)-1);
for  k = 1:length(h)
    Br2(k,:) = mu_0*conv(-didt,Hrstep10(k,:));
end
Br1 = zeros(length(h),length(t));
for  k = 1:length(h)
    Br1(k,:) = mu_0*Hrstep01(k,end)*n*current;
end
Br = Br1+Br2(:,1:length(t));
Ur2 = -Sz.*diff(Br2')'*fs;
Ur = -Sz.*diff(Br')'*fs;
%% phi
Ephi2 = zeros(length(h),length(t)+length(tob)-1);
for  k = 1:length(h)
    Ephi2(k,:) = conv(-didt,epstep10(k,:));
end
Ephi1 = zeros(length(h),length(t));
for  k = 1:length(h)
    Ephi1(k,:) = epstep01(k,end)*n*current;
end
sensor_len = 1;
Ephi = Ephi1+Ephi2(:,1:length(t));
Uphi2 = Ephi2*sensor_len;
Uphi = Ephi*sensor_len;
%%  ������ų���phi��ų�ת����x-y�����£����辶����x����н�Ϊsita
sita =pi/6;% ����sita = 30��
Ux = Ur(:,1:length(t)-1).*cos(sita)-Uphi(:,1:length(t)-1).*sin(sita);
Uy = Ur(:,1:length(t)-1).*sin(sita)+Uphi(:,1:length(t)-1).*cos(sita);

Ux2 = Ur2(:,1:length(t)-1).*cos(sita)-Uphi2(:,1:length(t)-1).*sin(sita);
Uy2 = Ur2(:,1:length(t)-1).*sin(sita)+Uphi2(:,1:length(t)-1).*cos(sita);


%% ��ֱ��Ӧ��ѹ
figure;
plot(t(1:Ns-1),1e3*Uz(:,1:Ns-1),'Linewidth',1.5);
title(['rou =' num2str(rou) ' �շ���Ȧ��ֱ����' num2str(L) '��-Uz']);
grid on;
% legend('h=50 Uz','h=100 Uz ','h=150 Uz ','h=200 Uz','h=250 Uz' ,'h=300 Uz','h=350 Uz','h=400 Uz','h=450 Uz','h=500 Uz');
legend('h=100','h=200' ,'h=300','h=400','h=500');
xlabel('time/s');
ylabel('Uz/mV');
figure;
plot(t(1:Ns-1),1e3*Uz2(:,1:Ns-1),'Linewidth',1.5);
grid on;
title(['rou =' num2str(rou) ' �շ���Ȧ��ֱ����' num2str(L) '��-Uz2']);
legend('h=100','h=200' ,'h=300','h=400','h=500');
% legend('h=50 Uz2','h=100 Uz2 ','h=150 Uz2 ','h=200 Uz2','h=250 Uz2' ,'h=300 Uz2','h=350 Uz2','h=400 Uz2','h=450 Uz2','h=500 Uz2');
xlabel('time/s');
ylabel('Uz2/mV');
%% �����ѹ
figure;
plot(t(1:Ns-1),1e3*Ur(:,1:Ns-1),'Linewidth',1.5);
grid on;
title(['rou =' num2str(rou) ' �շ���Ȧ��ֱ����' num2str(L) '��-Ur']);
legend('h=100','h=200' ,'h=300','h=400','h=500');
% legend('h=50 Ur','h=100 Ur ','h=150 Ur ','h=200 Ur','h=250 Ur' ,'h=300 Ur','h=350 Ur','h=400 Ur','h=450 Ur','h=500 Ur');
xlabel('time/s');
ylabel('Ur/mV');
figure;
plot(t(1:Ns-1),1e3*Ur2(:,1:Ns-1),'Linewidth',1.5);
grid on;
title(['rou =' num2str(rou) ' �շ���Ȧ��ֱ����' num2str(L) '��-Ur2']);
legend('h=100','h=200' ,'h=300','h=400','h=500');
% legend('h=50 Ur2','h=100 Ur2 ','h=150 Ur2 ','h=200 Ur2','h=250 Ur2' ,'h=300 Ur2','h=350 Ur2','h=400 Ur2','h=450 Ur2','h=500 Ur2');
xlabel('time/s');
ylabel('Ur2/mV');
%% phi ���ѹ
figure;
plot(t(1:Ns-1),1e3*Uphi(:,1:Ns-1),'Linewidth',1.5);
grid on;
title(['rou =' num2str(rou) ' �շ���Ȧ��ֱ����' num2str(L) '��-Uphi']);
legend('h=100','h=200' ,'h=300','h=400','h=500');
% legend('h=50 Uphi','h=100 Uphi ','h=150 Uphi','h=200 Uphi','h=250 Uphi' ,'h=300 Uphi','h=350 Uphi','h=400 Uphi','h=450 Uphi','h=500 Uphi');
xlabel('time/s');
ylabel('Uphi/mV');
figure;
plot(t(1:Ns-1),1e3*Uphi2(:,1:Ns-1),'Linewidth',1.5);
grid on;
title(['rou =' num2str(rou) ' �շ���Ȧ��ֱ����' num2str(L) '��-Uphi2']);
legend('h=100','h=200' ,'h=300','h=400','h=500');
% legend('h=50 Uphi2','h=100 Uphi2 ','h=150 Uphi2','h=200 Uphi2','h=250 Uphi2' ,'h=300 Uphi2','h=350 Uphi2','h=400 Uphi2','h=450 Uphi2','h=500 Uphi2');
xlabel('time/s');
ylabel('Uphi2/mV');
%% x-y�����Ӧ��ѹ
figure;
plot(t(1:Ns-1),1e3*Ux(:,1:Ns-1),'Linewidth',1.5);
grid on;
title(['rou =' num2str(rou) ' �շ���Ȧ��ֱ����' num2str(L) '��-Uphi']);
legend('h=100','h=200' ,'h=300','h=400','h=500');
% legend('h=50 Uphi','h=100 Uphi ','h=150 Uphi','h=200 Uphi','h=250 Uphi' ,'h=300 Uphi','h=350 Uphi','h=400 Uphi','h=450 Uphi','h=500 Uphi');
xlabel('time/s');
ylabel('Uphi/mV');
figure;
plot(t(1:Ns-1),1e3*Ux2(:,1:Ns-1),'Linewidth',1.5);
grid on;
title(['rou =' num2str(rou) ' �շ���Ȧ��ֱ����' num2str(L) '��-Uphi2']);
legend('h=100','h=200' ,'h=300','h=400','h=500');
% legend('h=50 Uphi2','h=100 Uphi2 ','h=150 Uphi2','h=200 Uphi2','h=250 Uphi2' ,'h=300 Uphi2','h=350 Uphi2','h=400 Uphi2','h=450 Uphi2','h=500 Uphi2');
xlabel('time/s');
ylabel('Uphi2/mV');
%% x-y�����Ӧ��ѹ
figure;
plot(t(1:Ns-1),1e3*Uy(:,1:Ns-1),'Linewidth',1.5);
grid on;
title(['rou =' num2str(rou) ' �շ���Ȧ��ֱ����' num2str(L) '��-Uphi']);
legend('h=100','h=200' ,'h=300','h=400','h=500');
% legend('h=50 Uphi','h=100 Uphi ','h=150 Uphi','h=200 Uphi','h=250 Uphi' ,'h=300 Uphi','h=350 Uphi','h=400 Uphi','h=450 Uphi','h=500 Uphi');
xlabel('time/s');
ylabel('Uphi/mV');
figure;
plot(t(1:Ns-1),1e3*Uy2(:,1:Ns-1),'Linewidth',1.5);
grid on;
title(['rou =' num2str(rou) ' �շ���Ȧ��ֱ����' num2str(L) '��-Uphi2']);
legend('h=100','h=200' ,'h=300','h=400','h=500');
% legend('h=50 Uphi2','h=100 Uphi2 ','h=150 Uphi2','h=200 Uphi2','h=250 Uphi2' ,'h=300 Uphi2','h=350 Uphi2','h=400 Uphi2','h=450 Uphi2','h=500 Uphi2');
xlabel('time/s');
ylabel('Uphi2/mV');
