% %% ����ˮƽ�ӵس�����Դ˲������Ӧ
% ���ø�˹���õ»��ֵķ���������ɢ���ӵس�����Դ
%
%sigma1 :��һ��ĵ絼��
% �����Ӻ���
%{
[ temp_e1,temp_e2,temp_h1,temp_h2] = calculate_temp(r1,r2,z,t,G_S,m2,J_1,delta_1)
%}

function [hz_01,hz_10,hz_impulse,hx_01,hx_10,hx_impulse,hy_01,hy_10,hy_impulse,ex_01,ex_10,ex_impulse,ey_01,ey_10,ey_impulse,ez_01,ez_impulse] = ...
    Calculate_Horizontal_Finite_Electrical_Source_GuaLeg_out(I,L,h,x,y,z,r_array,Ak,t,n)
% ================================
% I �������ǿ�ȣ�
% L �ӵس����ߵĳ��ȣ�
% (x,y,z)�۲������ꣻ
% t�۲�ʱ��
% h ��ʾ�ӵس����߾����ľ��룬����h = 0;
% ================================
u0 = 4*pi*1e-7;
% r = (x^2+y^2).^(0.5);
%% ����ӵ���Ĳ���
R1 = ((x+L/2).^2+y.^2).^0.5; 
R2 = ((x-L/2).^2+y.^2).^0.5;
%% define the part two of the ex and  hy, in addition to  the hz  
% part 2 of the ex
ex_01_par2 =  zeros(1,length(t));
ex_impulse_par2 =  zeros(1,length(t));
% part 2 of the hy
hy_01_par2 =  zeros(1,length(t));
hy_10_par2 = zeros(1,length(t));
hy_impulse_par2 =  zeros(1,length(t));
% hz
hz_01 = zeros(1,length(t));% the + step response
hz_10 = zeros(1,length(t));% the - step response
hz_impulse = zeros(1,length(t));% the + impulse response
%% 
G_S=load ('G_S.txt')';% G_S������
m2 = 1:length(G_S);
%--------------------------------------------------------------------------
%��һ������ȡ�Ѿ��洢���˲���ϵ��,��ʾΪ��������
%%  ����e0,h00,h01
load J0_Gupt.txt;       
J_zero = J0_Gupt( :, 3)'; % ���ٺ��˶��任�˲�ϵ��
delta = J0_Gupt( :, 2)'; %  ������ĺ�����ƫ����
%% e1 h1 �������
 load J1_Gupt.txt;       
J_1 = J1_Gupt( :, 3)'; % ���ٺ��˶��任�˲�ϵ��
delta_1 = J1_Gupt( :, 2)'; %  ������ĺ�����ƫ����
% �����м�ֵ
% ��Ӧ�ĵ��е��м���C_e1~e11;
% ��Ӧ�ĵ��е��м���C_e2~e12;
% ��Ӧ�ĵ��е��м���C_h1~h11;
% ��Ӧ�ĵ��е��м���C_h2~h12;
% [ e11_impulse,e11_step,e12_impulse,e12_step,h11_impulse,h11_step,h12_impulse,h12_step,h11_10_step,h12_10_step] =...
%     calculate_temp(R1,R2,z,t,G_S,m2,J_1,delta_1);
[ e11_impulse,e11_step,e12_impulse,e12_step,h11_impulse,h11_step,h12_impulse,h12_step,h11_10_step,h12_10_step,ez_01_step,ez_1_impulse] =...
    calculate_temp(R1,R2,z,t,G_S,m2,J_1,delta_1,J_zero,delta);
% ����Ey-Hx
ey_01 =(I.*y./(4.*pi).*(1./R2.*e12_step-1./R1.*e11_step));
ey_impulse = (I.*y./(4.*pi).*(1./R2.*e12_impulse-1./R1.*e11_impulse));% ����Ծ��Ӧ
hx_01 = -I.*y./(4.*pi).*(1./R2.*h12_step-1./R1.*h11_step);% + step response,only the second response
% ����Hx��һ�γ�
HX0 = -I.*y.*(1./(sqrt(z.^2+R2.^2).*(sqrt(z.^2+R2.^2)-z)) - 1./(sqrt(z.^2+R1.^2).*(sqrt(z.^2+R1.^2)-z)));
% hx_10 = HX0-hx_01;
hx_10 = -I.*y./(4.*pi).*(1./R2.*h12_10_step-1./R1.*h11_10_step);
%
hx_impulse = -I.*y./(4.*pi).*(1./R2.*h12_impulse-1./R1.*h11_impulse);% + step response,only the second response 
% ����Ex-Hy����ӵ����йصĲ���
ex_01_par1 =  ( I./(4.*pi).*((x-L/2)./R2.*e12_step-(x+L/2)./R1.*e11_step));
ex_impulse_par1 =  (I./(4.*pi).*((x-L/2)./R2.*e12_impulse-(x+L/2)./R1.*e11_impulse));
% Ez
ez_01 = I./4./pi.*ez_01_step;
ez_impulse = I./4./pi.*ez_1_impulse;
%Hy
hy_01_par1 = I./(4.*pi).*((x-L./2)./R2.*h12_step-(x+L./2)./R1.*h11_step);
hy_10_par1 = -hy_01_par1;
hy_impulse_par1 = I./(4.*pi).*((x-L./2)./R2.*h12_impulse-(x+L./2)./R1.*h11_impulse);
%% ������Դ��ɢ��----�õ���ɢ����ĵ�ż��Դ��
% ���ø�˹���õ»��ַ��� �õ���ɢ�����
dx = L./2;
% obs_point = [x,y,z];
% tol = 1e-10;
% [Ak,xk,dxk,r_array] = GuaLeg_DiscreteSource(obs_point,L,n,tol);
% dx��λ�õ�ͬ��Ak��λ��
for kk = 1:length(r_array) % �۲�����ɢԴ��ˮƽ
    r = r_array(kk);
    coef = Ak(kk);
%% -----------------------------------------�����س����ߵĻ�����---------------------------------------------------------------
%����lambda������lambda��frequency��չ�ɶ�ά����-----��Դ�ֱ���򺺿˶��任
lambda_0 = (1./r) .*exp(delta); % ��μ���lambda���ɲ�����ĺ�����ƫ����ת��Ϊ���ֱ���lambda�����J0
lambda_1 = (1./r) .*exp(delta_1); % ��μ���lambda���ɲ�����ĺ�����ƫ����ת��Ϊ���ֱ���lambda,���J1���˶��任
     for ii=1:length(t)
            freq = (log(2)*1i/(t(ii)*2*pi))*m2;
            %--------------------------------------------------------------------------
            %���ݵ��ƹ�ʽ�������-����ķ���ϵ��
            [lambda0_Array,frequency0_Array] = Array_trans(lambda_0,freq); % ������J0���˶��任
            [lambda1_Array,frequency1_Array] = Array_trans(lambda_1,freq); % ������J1���˶��任   
            r_TE0=calculate_r_TE(lambda_0,freq); % ���J0
            r_TE1=calculate_r_TE(lambda_1,freq); % ���J1
    %         z1bar_1 = calculate_r_TM_Z1bar(lambda_1,freq);% ���J1
    %         z1bar_0 = calculate_r_TM_Z1bar(lambda_0,freq);% ���J0
            %h��Դλ��,zΪ�۲�λ��
            %% -----------------------------------�뺽��tem �糡���������� h=0----------------------------------------
    %         w0 = 2.*pi.*frequency0_Array;%����
    %         w1 = 2.*pi.*frequency1_Array;%����
             w =  2.*pi.*freq';
             z0_bar = 1i.*w.*u0;
            %% ���Ex-part 2
            f_e0 = (1+r_TE0).*exp(lambda0_Array.*(z));
            g_e0 = Fast_Hankel(r,f_e0,J_zero);
            ex_impulse_par2(ii) =ex_impulse_par2(ii) +coef.* (-I./(4.*pi)).*GS_Trans2(t(ii),z0_bar.*g_e0,G_S); %��������Ӧʱ��
            ex_01_par2(ii) =ex_01_par2(ii) + coef.*(-I./(4.*pi)).*GS_Trans(t(ii),z0_bar.*g_e0,freq,G_S);%����Ծ��Ӧʱ��
            
            %% -----------------------------------�뺽��tem �ų����������� h=0----------------------------------------
            %% ���Hz
                f_hz1 = (1+r_TE1).*exp(lambda1_Array.*(z)).*lambda1_Array;% r_TM1 = 1;
                primary_hz1 = exp(lambda1_Array.*(z)).*lambda1_Array;%% ��Ӧֱ���ɷ֣���ʱ�ĺ˺����е�r_TE=0  
                g_hz1 = Fast_Hankel(r,f_hz1,J_1);%����Ծ��ӦƵ��
                g_hz1_zeros = Fast_Hankel(r,primary_hz1,J_1);%����Ծ��Ӧ��Ƶ��Ӧ
    %         h1_0_impulse(ii) = -GS_Trans2(t(ii),g_h1,G_S);%��������Ӧʱ��
    %         h1_10(ii) = +GS_Trans(t(ii), g_h1_zeros,freq,G_S)-GS_Trans(t(ii),g_h1,freq,G_S);%����Ծ��Ӧʱ��
                hz_10(ii) =hz_10(ii)  + coef.*(I.*y./(4.*pi)).*(GS_Trans(t(ii), g_hz1_zeros,freq,G_S)-GS_Trans(t(ii),g_hz1,freq,G_S))./r;%����Ծ��Ӧʱ��
                hz_impulse(ii) =hz_impulse(ii) + coef.*(I.*y./(4.*pi)).*GS_Trans2(t(ii),g_hz1,G_S)./r; %��������Ӧʱ��
                hz_01(ii) =hz_01(ii) + coef.*(I.*y./(4.*pi)).*GS_Trans(t(ii),g_hz1,freq,G_S)./r;%����Ծ��Ӧʱ��

           %% ��� Hy part 2
%            y0 = 2.*exp(lambda0_Array.*(z)).*lambda0_Array;
%            primary_h00 =exp(lambda0_Array.*(z)).*lambda0_Array;%% ��Ӧֱ���ɷ֣���ʱ�ĺ˺����е�r_TE=0
           f_h00 = (1+r_TE0).*exp(lambda0_Array.*(z)).*lambda0_Array;     
           %{
           g_h00 = Fast_Hankel(r,f_h00,J_zero);%����Ծ��ӦƵ��
           g_h00_zeros = Fast_Hankel(r,primary_h00,J_zero);%����Ծ��Ӧ��Ƶ��Ӧ
           %}
           g_h00 = Fast_Hankel(r,f_h00,J_zero);%����Ծ��ӦƵ��
%            g_h00_zeros = Fast_Hankel(r,primary_h00,J_zero);%����Ծ��Ӧ��Ƶ��Ӧ
%            hy_10_par2(ii) = hy_10_par2(ii)+(-I./(4.*pi)).*(GS_Trans(t(ii), g_h00_zeros,freq,G_S)-GS_Trans(t(ii),g_h00,freq,G_S));%����Ծ��Ӧʱ��
           hy_impulse_par2(ii) =hy_impulse_par2(ii) + coef.*(I./(4.*pi)).*GS_Trans2(t(ii),g_h00,G_S); %��������Ӧʱ��
           hy_01_par2(ii) =hy_01_par2(ii) +coef.* (I./(4.*pi)).*GS_Trans(t(ii),g_h00,freq,G_S);%����Ծ��Ӧʱ��         
     end
end
    %% �糡����
    % ex
    ex_01 = ex_01_par1 + ex_01_par2.*dx;% ����Ծ
    ex_10 = ex_01(end) - ex_01;% ����Ծ
    ex_impulse = ex_impulse_par1 + ex_impulse_par2.*dx;   
    % ey
    ey_10 = ey_01(end) - ey_01;% ����Ծ
    %%  �ų�����  
    % ����Hy��һ�γ�
     HY0 = -I./4./pi.*z./(y.^2+z.^2).*((L/2-x)./sqrt(z.^2+R2.^2)+(L/2+x)./sqrt(z.^2+R1.^2)) +...
         I./4./pi.*((x-L/2)./sqrt(z.^2+R2.^2)./(sqrt(z.^2+R2.^2)-z) -(x+L/2)./sqrt(z.^2+R1.^2)./(sqrt(z.^2+R1.^2)-z));
    hy_01 = hy_01_par1 + hy_01_par2.*dx;% ����Ծ��Ӧ
%     hy_10 = hy_10_par1 + hy_10_par2.*dx;
%     hy_10 = hy_01(end) - hy_01 ;
        hy_10 = HY0 - hy_01 ;% ���㸺��Ծ��Ӧ
    hy_impulse = hy_impulse_par1 + hy_impulse_par2.*dx;
    %Hz
    hz_10 = hz_10.*dx;
    hz_01 = hz_01.*dx;
%     HZ0 = I.*y./4./pi./(y.^2+z.^2).*( (L./2-x)./sqrt(z.^2+R2.^2) + (L./2 +x)./sqrt(z.^2+R1.^2));
%     hz_10 = HZ0 - hz_01;
    hz_impulse = hz_impulse.*dx;
%         load parameters.txt;% ��ȡ����
%        save(['SemiAtem_Horizontal_Finite_Electrical_Source_L' num2str(L) '_h' num2str(h) '_z' num2str(z) '_x' num2str(x) '_y' num2str(y) '.mat'],...
%         'hz_01','hz_10','hz_impulse','hx_01','hx_10','hx_impulse','hy_01','hy_10','hy_impulse','ex_01','ex_10','ex_impulse','ey_01','ey_10','ey_impulse','parameters');
end