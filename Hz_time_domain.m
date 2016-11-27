%--------------------------------------------------------------------------
%��������������˲�ϴ�ֱ��ż���Ĵ�ֱ�ų���ʱ��仯���ɡ�
%��ż�����ڵ�����100 Ohm*m �ľ��ȴ�ر��棬�۲�����ż��100m��
%--------------------------------------------------------------------------
%%
clear all;clc;close all;
format long;
%%
mu_0 = 4*pi*1e-7;
%% atem �շ��߶Ȳ���
h = 100;% Դ�����ĸ߶�
z = 50;% �۲������ĸ߶�
%% ����
r=100; % �շ���
G_S=load ('G_S.txt')';
m2 = 1:length(G_S);
fs = 1e7;
t = 1/fs:1/fs:1e-2;% ʱ������
% t = logspace(-8,-1,1000);
h_z1_t=zeros(1,length(t));
h_z2_t=h_z1_t;
%  for ii=1:length(t)
%     freq = (log(2)*1i/(t(ii)*2*pi))*m2;
    %--------------------------------------------------------------------------
    %��һ������ȡ�Ѿ��洢���˲���ϵ��,��ʾΪ��������
    load J0_Gupt.txt;       
    J_zero = J0_Gupt( :, 3)'; % ���ٺ��˶��任�˲�ϵ��
    delta = J0_Gupt( :, 2)'; %  ������ĺ�����ƫ����
    %--------------------------------------------------------------------------
    %����lambda������lambda��frequency��չ�ɶ�ά����
    lambda = (1./r) .*exp(delta); % ��μ���lambda���ɲ�����ĺ�����ƫ����ת��Ϊ���ֱ���lambda
for ii=1:length(t)
    freq = (log(2)*1i/(t(ii)*2*pi))*m2;
    [lambda_Array,frequency_Array] = Array_trans(lambda,freq);
    %--------------------------------------------------------------------------
    %���ݵ��ƹ�ʽ�������-����ķ���ϵ��
    r_TE=calculate_r_TE(lambda,freq);
    %% --------------------------------------------------------------------------
    %����ų��Ĵ�ֱ���������ÿ��ٺ��˶��任����ֵ���ֵ�⡣  
    %��ֱ��żԴz��ų��ĺ˺���,h��Դλ��,zΪ�۲�λ��
    % ��ֱ��ż��Դz��ų��ĺ˺��� u_0 = lambda,z=h=0,���淢��������;
    sum = (1+r_TE).*lambda_Array.^2;  % ����Ծ��ӦƵ��˺���
    sum_zeros = (1).*lambda_Array.^2; %% ��Ӧֱ���ɷ֣���ʱ�ĺ˺����е�r_TE=0
    H_vertical = 1./(4*pi) *  Fast_Hankel(r,sum,J_zero);%����Ծ��ӦƵ��
    H_vertical_zeros = 1./(4*pi) *  Fast_Hankel(r,sum_zeros,J_zero);%����Ծ��Ӧ��Ƶ��Ӧ
    h_z1_t(ii) = -GS_Trans2(t(ii),H_vertical,G_S);%GS_Trans2(t(ii),H_vertical_zeros); %��������Ӧʱ��
    h_z2_t(ii) = +GS_Trans(t(ii),H_vertical_zeros,freq,G_S)-GS_Trans(t(ii),H_vertical,freq,G_S);%����Ծ��Ӧʱ��
    %% ---------------------------------------------------------------------------
    % �뺽�մ�ֱ��żԴz��ų��ĺ˺���,h��Դλ��,zΪ�۲�λ��.���ڰ뺽��,z=0,h=-100�ס�
    % sum = exp(-lambda_Array*h).*(1+r_TE).*lambda_Array.^2;
    %% -----------------------------------atem ��ֱ����----------------------------------------
    % ȫ����atem��ֱ��żԴz��ų��ĺ˺���,h��Դλ��,zΪ�۲�λ��.����ȫ����,z=-50,h=-100�ס�
    sum = (exp(-lambda_Array*(z+h))+r_TE.*exp(lambda_Array*(z-h))).*lambda_Array.^2;
    sum_zeros = (1).*lambda_Array.^2; %% ��Ӧֱ���ɷ֣���ʱ�ĺ˺����е�r_TE=0
    H_vertical = 1./(4*pi) *  Fast_Hankel(r,sum,J_zero);%����Ծ��ӦƵ��
    H_vertical_zeros = 1./(4*pi) *  Fast_Hankel(r,sum_zeros,J_zero);%����Ծ��Ӧ��Ƶ��Ӧ
    h_z1_t(ii) = -GS_Trans2(t(ii),H_vertical,G_S);%GS_Trans2(t(ii),H_vertical_zeros); %��������Ӧʱ��
    h_z2_t(ii) = +GS_Trans(t(ii),H_vertical_zeros,freq,G_S)-GS_Trans(t(ii),H_vertical,freq,G_S);%����Ծ��Ӧʱ��
    %% -----------------------------------atem ����----------------------------------------
    %% -----------------------------------atem phi����-----------------------------------
 end
%% --------------------------------------------------------------------------
%��ͼ��
loglog(t.*10^3,(h_z1_t),'r','Linewidth',2)
hold on
% axis([10^-5,10^3,10^-11,10^-1])
loglog(t.*10^3,0.5*(abs(h_z1_t)-h_z1_t),'r--','Linewidth',2)
hold on
% axis([10^-5,10^3,10^-11,10^-1])
loglog(t.*10^3,(h_z2_t),'b','Linewidth',2)
hold on 
% axis([10^-5,10^3,10^-11,10^-1])
loglog(t.*10^3,0.5*(abs(h_z2_t)-h_z2_t),'b--','Linewidth',2)
hold on 
% axis([10^-5,10^3,10^-11,10^-1])
legend('��������Ӧ','','����Ծ��Ӧ','');
title('˲�ϴ�ֱ��ż���ӵĴ�ֱ�ų�����ʱ�䵼����ʱ��仯������')
xlabel('ʱ�䣨ms��')
%% ���渺������Ӧ
himpulse10 = h_z1_t;
hstep10 = h_z2_t;
save('atem_magnetic_dipole_response.mat','himpulse10','hstep10');
%%
