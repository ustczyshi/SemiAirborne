%--------------------------------------------------------------------------
%��������������TEM �뺽��TEM
% ����ˮƽ��żԴ�������������й۲�
%--------------------------------------------------------------------------
%%
format long;
clear all;clc;close all;

%   for k2 = 1:5
% d = [0 25 30 35 50];
%%
mu_0 = 4*pi*1e-7;
rou = 1;
%% ����
r=50; % �շ�ˮƽƫ�ƾ�
a = 10; % Բ�λ��߰뾶
% �뺽���շ��߶Ȳ���
z =0;% �۲������ĸ߶ȣ���������Ϊ��ֵ
h =0;% Դ�����ĸ߶�


%% ȫ�����շ��߶Ȳ���
%% ���泡��ѡ��
%{
disp('************************************************');
disp('����ģ��ѡ��');
disp('0�����ȴ�ص���TEM----------');
disp('1�����ȴ��ȫ����TEM----------');
disp('2�����ȴ�ذ뺽��TEM����----------');
tig = input('��ѡ��');
switch tig
    case 0
        r  = input('�շ�ˮƽ�ࣺ');
%         a = input('Բ�λ��߰뾶 ��');
        h = 0;
        z = 0;
    case 1  %  ���ȴ��ȫ����TEM
%         a = input('Բ�λ��߰뾶 ��');
        h = input('����Դ�����߶ȣ�');
        z = input('�۲���Ȧ�����߶ȣ�'); 
    case 2   %  ���ȴ�ذ뺽��TEM
        h = 0;
        z = input('�۲���Ȧ�����߶ȣ�'); 
    case -1
        h = 0;
        z = 0;
        r = 100;
    otherwise 
         error('����ѡ�����Ĭ�Ϸ������TEM ');
         tig = input('��ѡ��');
        return;
end
disp('************************************************');
%}
%% Դ����ģ��ѡ��: 1--��ֱ��ż��Դ��0��Բ�λ���Դ

% fig = 1;% ��ֱ��ż��Դ
 fig = 0;% Բ�λ���
%% 
G_S=load ('G_S.txt')';
m2 = 1:length(G_S);
fs = 1e7;
t = 1/fs:1/fs:1e-2;% ʱ������
% t = logspace(-8,1,10000);
%%  z�����
h_z1_t=zeros(1,length(t));% ��������Ӧʱ��
h_z2_t=h_z1_t;%����Ծ
h01_z1_t=zeros(1,length(t));% ��������Ӧʱ��
h01_z2_t=h_z1_t;%����Ծ
%%  x�����
h_x1_t=zeros(1,length(t));
h_x2_t=h_x1_t;
h01_x1_t=zeros(1,length(t));% ��������Ӧʱ��
h01_x2_t=h_x1_t;%����Ծ
%%  y����
h_y1_t=zeros(1,length(t));
h_y2_t=h_y1_t;
h01_y1_t=zeros(1,length(t));% ��������Ӧʱ��
h01_y2_t=h_y1_t;%����Ծ
%--------------------------------------------------------------------------
%��һ������ȡ�Ѿ��洢���˲���ϵ��,��ʾΪ��������
%%  y�������
load J0_Gupt.txt;       
J_zero = J0_Gupt( :, 3)'; % ���ٺ��˶��任�˲�ϵ��
delta = J0_Gupt( :, 2)'; %  ������ĺ�����ƫ����
%% x-z�������
 load J1_Gupt.txt;       
J_1 = J1_Gupt( :, 3)'; % ���ٺ��˶��任�˲�ϵ��
delta_1 = J1_Gupt( :, 2)'; %  ������ĺ�����ƫ����
%--------------------------------------------------------------------------
%����lambda������lambda��frequency��չ�ɶ�ά����-----��Դ�ֱ���򺺿˶��任
if fig == 1; % ˮƽ��ż��Դ
    lambda = (1./r) .*exp(delta); % ��μ���lambda���ɲ�����ĺ�����ƫ����ת��Ϊ���ֱ���lambda�����J0
    lambda_r = (1./r) .*exp(delta_1); % ��μ���lambda���ɲ�����ĺ�����ƫ����ת��Ϊ���ֱ���lambda,���J1���˶��任
end
%{
if fig ==0 ;% Բ�λ��ߣ���ֱ����ų����˶��任��lambdaֵ 
    lambda = (1./a) .*exp(delta_1); 
end
%}
%--------------------------------------------------------------------------

    for ii=1:length(t)
        freq = (log(2)*1i/(t(ii)*2*pi))*m2;
        %--------------------------------------------------------------------------
        %���ݵ��ƹ�ʽ�������-����ķ���ϵ��
        if fig == 1;
            [lambda_Array,frequency_Array] = Array_trans(lambda,freq); % ��Դ�ֱ���򺺿˶��任
            [lambdar_Array,frequencyr_Array] = Array_trans(lambda_r,freq); % ���ˮƽ���򺺿˶��任   
            r_TE=calculate_r_TE(lambda,freq); % ��Դ�ֱ���򺺿˶��任
            r_TEr=calculate_r_TE(lambda_r,freq); % ��Դ�ֱ���򺺿˶��任
        end
        if fig ==0;
            [lambda_Array,frequency_Array] = Array_trans(lambda,freq); % ��Դ�ֱ���򺺿˶��任
            r_TE=calculate_r_TE(lambda,freq); % ��Դ�ֱ���򺺿˶��任
        end
        %% --------------------------------------------------------------------------
        %{
        %����ų��Ĵ�ֱ���������ÿ��ٺ��˶��任����ֵ���ֵ�⡣  
        %��ֱ��żԴz��ų��ĺ˺���,h��Դλ��,zΪ�۲�λ��
        % ��ֱ��ż��Դz��ų��ĺ˺��� u_0 = lambda,z=h=0,���淢��������;
    %     sum = (1+r_TE).*lambda_Array.^2;  % ����Ծ��ӦƵ��˺���
    %     sum_zeros = (1).*lambda_Array.^2; %% ��Ӧֱ���ɷ֣���ʱ�ĺ˺����е�r_TE=0
    %     H_vertical = 1./(4*pi) *  Fast_Hankel(r,sum,J_zero);%����Ծ��ӦƵ��
    %     H_vertical_zeros = 1./(4*pi) *  Fast_Hankel(r,sum_zeros,J_zero);%����Ծ��Ӧ��Ƶ��Ӧ
    %     h_z1_t(ii) = -GS_Trans2(t(ii),H_vertical,G_S);%GS_Trans2(t(ii),H_vertical_zeros); %��������Ӧʱ��
    %     h_z2_t(ii) = +GS_Trans(t(ii),H_vertical_zeros,freq,G_S)-GS_Trans(t(ii),H_vertical,freq,G_S);%����Ծ��Ӧʱ��
        %}
        %% ---------------------------------------------------------------------------
        if fig == 1; % ��ֱ��ż��Դ
            %% -----------------------------------�뺽��tem ��ֱ����----------------------------------------
            % ȫ����atem��ֱ��żԴz��ų��ĺ˺���,h��Դλ��,zΪ�۲�λ��.����ȫ����,z=-50,h=-100�ס�
            sum = (exp(-lambda_Array*(z+h))+r_TE.*exp(lambda_Array*(z-h))).*lambda_Array.^2;
            sum_zeros = (1).*lambda_Array.^2.*exp(-lambda_Array*(z+h)); %% ��Ӧֱ���ɷ֣���ʱ�ĺ˺����е�r_TE=0
            H_vertical = 1./(4*pi) *  Fast_Hankel(r,sum,J_zero);%����Ծ��ӦƵ��
            H_vertical_zeros = 1./(4*pi) *  Fast_Hankel(r,sum_zeros,J_zero);%����Ծ��Ӧ��Ƶ��Ӧ
            h_z1_t(ii) = -GS_Trans2(t(ii),H_vertical,G_S);%GS_Trans2(t(ii),H_vertical_zeros); %��������Ӧʱ��
            h_z2_t(ii) = +GS_Trans(t(ii),H_vertical_zeros,freq,G_S)-GS_Trans(t(ii),H_vertical,freq,G_S);%����Ծ��Ӧʱ��
            h01_z1_t(ii) = GS_Trans2(t(ii),H_vertical,G_S);%GS_Trans2(t(ii),H_vertical_zeros); %��������Ӧʱ��
            h01_z2_t(ii) = GS_Trans(t(ii),H_vertical,freq,G_S);%����Ծ��Ӧʱ��
            %% -----------------------------------�뺽��tem ����----------------------------------------
            sumr = (exp(-lambdar_Array*(z+h))-r_TEr.*exp(lambdar_Array*(z-h))).*lambdar_Array.^2;
            sumr_zeros = (1).*lambdar_Array.^2.*exp(-lambdar_Array*(z+h)); %% ��Ӧֱ���ɷ֣���ʱ�ĺ˺����е�r_TE=0
            Hr_vertical = 1./(4*pi) *  Fast_Hankel(r,sumr,J_1);%����Ծ��ӦƵ��
            Hr_vertical_zeros = 1./(4*pi) *  Fast_Hankel(r,sumr_zeros,J_1);%����Ծ��Ӧ��Ƶ��Ӧ
            h_r1_t(ii) = -GS_Trans2(t(ii),Hr_vertical,G_S);%GS_Trans2(t(ii),H_vertical_zeros); %��������Ӧʱ��
            h_r2_t(ii) = +GS_Trans(t(ii),Hr_vertical_zeros,freq,G_S)-GS_Trans(t(ii),Hr_vertical,freq,G_S);%����Ծ��Ӧʱ��
            h01_r1_t(ii) = GS_Trans2(t(ii),Hr_vertical,G_S);%GS_Trans2(t(ii),H_vertical_zeros); %��������Ӧʱ��
            h01_r2_t(ii) = GS_Trans(t(ii),Hr_vertical,freq,G_S);%����Ծ��Ӧʱ��
            %% -----------------------------------�뺽��tem phi����糡-----------------------------------
            sump = (exp(-lambdar_Array*(z+h))+r_TEr.*exp(lambdar_Array*(z-h))).*lambdar_Array;
            sump_zeros = (1).*lambdar_Array.*exp(-lambdar_Array*(z+h)); %% ��Ӧֱ���ɷ֣���ʱ�ĺ˺����е�r_TE=0
            Hp_vertical = -1i*2*pi*mu_0*freq'./(4*pi) .* Fast_Hankel(r,sump,J_1);%����Ծ��ӦƵ��
            Hp_vertical_zeros = -1i*2*pi*mu_0*freq'./(4*pi) .* Fast_Hankel(r,sump_zeros,J_1);%����Ծ��Ӧ��Ƶ��Ӧ
            h_p1_t(ii) = -GS_Trans2(t(ii),Hp_vertical,G_S);%GS_Trans2(t(ii),H_vertical_zeros); %��������Ӧʱ��
            h_p2_t(ii) = +GS_Trans(t(ii),Hp_vertical_zeros,freq,G_S)-GS_Trans(t(ii),Hp_vertical,freq,G_S);%����Ծ��Ӧʱ��
            h01_p1_t(ii) = GS_Trans2(t(ii),Hp_vertical,G_S);%GS_Trans2(t(ii),H_vertical_zeros); %��������Ӧʱ��
            h01_p2_t(ii) = GS_Trans(t(ii),Hp_vertical,freq,G_S);%����Ծ��Ӧʱ��
        end
        if fig == 0; % Բ�λ���
            %% ***************************��ֱ���� z ��*********************************
            sum = (exp(-lambda_Array*(z+h))+r_TE.*exp(lambda_Array*(z-h))).*lambda_Array.*besselj(0,lambda_Array.*r);
            sum_zeros = exp(-lambda_Array*(z+h)).*lambda_Array.*besselj(0,lambda_Array.*r);%% ��Ӧֱ���ɷ֣���ʱ�ĺ˺����е�r_TE=0
            H_z = 0.5*a * Fast_Hankel(a,sum,J_1);%����Ծ��ӦƵ��
            H_z_zeros = 0.5*a *  Fast_Hankel(a,sum_zeros,J_1);
            h_z1_t(ii) = -GS_Trans2(t(ii),H_z,G_S);%GS_Trans2(t(ii),H_vertical_zeros); %��������Ӧʱ��
            h_z2_t(ii) = +GS_Trans(t(ii),H_z_zeros,freq,G_S)-GS_Trans(t(ii),H_z,freq,G_S);%����Ծ��Ӧʱ��
            h01_z1_t(ii) = GS_Trans2(t(ii),H_z,G_S);%GS_Trans2(t(ii),H_vertical_zeros); %��������Ӧʱ��
            h01_z2_t(ii) = GS_Trans(t(ii),H_z,freq,G_S);%����Ծ��Ӧʱ��
            %% **************************���� r �ų�**********************************************
            sumr = (exp(-lambda_Array*(z+h))-r_TE.*exp(lambda_Array*(z-h))).*lambda_Array.*besselj(1,lambda_Array.*r);
            sumr_zeros = exp(-lambda_Array*(z+h)).*lambda_Array.*besselj(1,lambda_Array.*r);%% ��Ӧֱ���ɷ֣���ʱ�ĺ˺����е�r_TE=0
            H_r = 0.5*a * Fast_Hankel(a,sumr,J_1);%����Ծ��ӦƵ��
            H_r_zeros = 0.5*a *  Fast_Hankel(a,sumr_zeros,J_1);
            h_r1_t(ii) = -GS_Trans2(t(ii),H_r,G_S);%GS_Trans2(t(ii),H_vertical_zeros); %��������Ӧʱ��
            h_r2_t(ii) = +GS_Trans(t(ii),H_r_zeros,freq,G_S)-GS_Trans(t(ii),H_r,freq,G_S);%����Ծ��Ӧʱ��
            h01_r1_t(ii) = GS_Trans2(t(ii),H_r,G_S);%GS_Trans2(t(ii),H_vertical_zeros); %��������Ӧʱ��
            h01_r2_t(ii) = GS_Trans(t(ii),H_r,freq,G_S);%����Ծ��Ӧʱ��
%             %% **************************phi����糡**********************************************
            sump = (exp(-lambda_Array*(z+h))+r_TE.*exp(lambda_Array*(z-h))).*besselj(1,lambda_Array.*r);
            sump_zeros = exp(-lambda_Array*(z+h)).*besselj(1,lambda_Array.*r);%% ��Ӧֱ���ɷ֣���ʱ�ĺ˺����е�r_TE=0
            H_p = -1i*2*pi*mu_0*freq'.*0.5*a .* Fast_Hankel(a,sump,J_1);%����Ծ��ӦƵ��
            H_p_zeros = -1i*2*pi*mu_0*freq'.*0.5*a .*  Fast_Hankel(a,sump_zeros,J_1);
            h_p1_t(ii) = -GS_Trans2(t(ii),H_p,G_S);%GS_Trans2(t(ii),H_vertical_zeros); %��������Ӧʱ��
            h_p2_t(ii) = +GS_Trans(t(ii),H_p_zeros,freq,G_S)-GS_Trans(t(ii),H_p,freq,G_S);%����Ծ��Ӧʱ��
            h01_p1_t(ii) = GS_Trans2(t(ii),H_p,G_S);%GS_Trans2(t(ii),H_vertical_zeros); %��������Ӧʱ��
            h01_p2_t(ii) = GS_Trans(t(ii),H_p,freq,G_S);%����Ծ��Ӧʱ��
        end
    end
 %% ������������Ӧ������Ծ��Ӧ
hzimpulse01 = h01_z1_t;
hzstep01 = h01_z2_t;
hrimpulse01 = h01_r1_t;
hrstep01 = h01_r2_t;
epimpulse01 = h01_p1_t;
epstep01 = h01_p2_t;
 %% ���渺������Ӧ�͸���Ծ��Ӧ
hzimpulse10 = h_z1_t;
hzstep10 = h_z2_t;
hrimpulse10 = h_r1_t;
hrstep10 = h_r2_t;
epimpulse10 = h_p1_t;
epstep10 = h_p2_t;
%% ������������
if fig ==1 
    save(['atem_magnetic_dipole_response_h' num2str(h) 'z' num2str(z) '_r' num2str(r) '_rou' num2str(rou) '.mat'],...
        'hzimpulse10','hzstep10','hrimpulse10','hrstep10','epimpulse10','epstep10',...
        'hzimpulse01','hzstep01','hrimpulse01','hrstep01','epimpulse01','epstep01');
end
if fig ==0 
    %% h0z0:���淢��������--���Ļ���
    save(['atem_circle_loop_response_h' num2str(h) 'z' num2str(z) '_r' num2str(r) '_rou' num2str(rou) '.mat'],...
        'hzimpulse10','hzstep10','hrimpulse10','hrstep10','epimpulse10','epstep10',...
        'hzimpulse01','hzstep01','hrimpulse01','hrstep01','epimpulse01','epstep01');
end
%% --------------------------------------------------------------------------
% ������
I = 1;
rou = 100;
 [step_H01,step_H10,impluse_H10] = Central_Circle_Loop_Step_Response(mu_0,rou,I,a,t);
 
%% ������͸���Ծ
%��ͼ��
figure;
loglog(t.*10^3,abs(h_z1_t),'r','Linewidth',1);
hold on
loglog(t.*10^3,(h_z2_t),'b','Linewidth',1);
hold on 
loglog(t.*1000,abs(diff([0 step_H10])*fs),'r--','Linewidth',1);
hold on;
loglog(t.*1000,step_H10,'b--','Linewidth',1);
hold on;
grid on;
legend('��ֵ����������Ӧ','��ֵ������Ծ��Ӧ','��������������Ӧ','����������Ծ��Ӧ');
title('����ԾԲ�λ��ߵĴ�ֱ�ų�')
xlabel('ʱ�䣨ms��')
ylabel('Hz/(A/m)');
%% ����Ծ
%��ͼ��
figure;
loglog(t.*10^3,abs(h01_z1_t),'r','Linewidth',1);
hold on
loglog(t.*10^3,(h01_z2_t),'b','Linewidth',1);
hold on 
loglog(t(1:end-1).*1000,abs(diff(step_H01)*fs),'r--','Linewidth',1);
hold on;
loglog(t.*1000,step_H01,'b--','Linewidth',1);
hold on;
grid on;
legend('��ֵ����������Ӧ','��ֵ������Ծ��Ӧ','��������������Ӧ','����������Ծ��Ӧ');
title('����ԾԲ�λ��ߵĴ�ֱ�ų�');
xlabel('ʱ�䣨ms��');
ylabel('Hz/(A/m)');
%% --------------------------------------------------------------------------
%��ͼ��
figure;
loglog(t.*10^3,abs(h_r1_t),'r','Linewidth',2);
hold on;
loglog(t.*10^3,0.5*(abs(h_r1_t)-h_r1_t),'r--','Linewidth',2);
hold on;
loglog(t.*10^3,abs(h_r2_t),'b','Linewidth',2);
hold on ;
loglog(t.*10^3,0.5*(abs(h_r2_t)-h_r2_t),'b--','Linewidth',2);
hold on ;
grid on;
title('˲��Բ�λ��ߵľ���ų�');
xlabel('ʱ�䣨ms��');
ylabel('Hr/(A/m)');
%% --------------------------------------------------------------------------
%��ͼ��
figure;
loglog(t.*10^3,abs(h_p1_t),'r','Linewidth',2);
hold on;
loglog(t.*10^3,0.5*(abs(h_p1_t)-h_p1_t),'r--','Linewidth',2);
hold on;
loglog(t.*10^3,abs(h_p2_t),'b','Linewidth',2);
hold on ;
loglog(t.*10^3,0.5*(abs(h_p2_t)-h_p2_t),'b--','Linewidth',2);
hold on ;
grid on;
title('˲��Բ�λ��ߵ�phi�糡');
xlabel('ʱ�䣨ms��');
ylabel('Ep/(A/m)');
close all;clear all;
%  end

