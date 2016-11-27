%--------------------------------------------------------------------------
%�������������洹ֱ��ż����ˮƽ�ų���Ƶ�ʱ仯���ɡ�
%��ż�����ڵ�����100 Ohm*m �ľ��ȴ�ر��棬�۲�����ż��100m��
%--------------------------------------------------------------------------
format long;
mu_0 = 4*pi*1e-7;
freq= logspace(-1,5,250);
r = 100;
%--------------------------------------------------------------------------
%��һ������ȡ�Ѿ��洢���˲���ϵ��,��ʾΪ��������
load FHT_J1_filter_Kong.txt;       
J_one = FHT_J1_filter_Kong( :, 3)';
delta = FHT_J1_filter_Kong( :, 2)';
%--------------------------------------------------------------------------
%����lambda������lambda��frequency��չ�ɶ�ά����
lambda=(1./r) .*exp(delta);
[lambda_Array,frequency_Array] = Array_trans(lambda,freq);
%--------------------------------------------------------------------------
% %���ݵ��ƹ�ʽ�����һ��ķ���ϵ��
r_TE = calculate_r_TE(lambda,freq);
%--------------------------------------------------------------------------
%����ų���ˮƽ���������ÿ��ٺ��˶��任����ֵ���ֵ�⡣  
sum = (1-r_TE).*lambda_Array.^2 ;
H_horizontal = 1/(4*pi) *  Fast_Hankel(r,sum,J_one);
%--------------------------------------------------------------------------
%��ͼ����ɫ��ʾʵ��������ɫ��ʾ�����������Ϊ��ֵ��ʵ��Ϊ��ֵ��
loglog(freq,0.5*(abs(real(H_horizontal))-real(H_horizontal)),'r--','Linewidth',2)
hold on
loglog(freq,(real(H_horizontal)),'r','Linewidth',2)
hold on
loglog(freq,0.5*(abs(imag(H_horizontal))-imag(H_horizontal)),'b--','Linewidth',2)
hold on 
loglog(freq,(imag(H_horizontal)),'b','Linewidth',2)
xlabel('Ƶ�ʣ�Hz��')
ylabel('Hz(A/m)')
title('��ֱ��ż����ˮƽ�ų���Ƶ�ʵı仯����')
axis([10^-1,10^5,10^-12,10^-6])