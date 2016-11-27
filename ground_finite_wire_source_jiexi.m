%% ���ӵص��߳�������ϵĴ�ֱ�ų�ʱ�䵼���Ľ�����
% ������֤ ��ֵ���󳤵���˲�䳡��Ӧ���ݵ���ȷ��
% �۲��@(0,y,0)
% ������Դ��Χ(-L/2,L/2)

function [ hz_impulse_jiexijie] = ground_finite_wire_source_jiexi(u0,rou,I,L,yr,Tob)
r = ((L/2).^2+yr.^2).^0.5;
% alph = (u0./rou).^0.5.*r;
sita = (u0./rou./Tob).^(0.5)./2;
temp_L = sita.*L./2;
temp_r = sita.*r;
temp = sita.*yr;
% ��������
erf_L = erf(temp_L);
erf_r = erf(temp_r);
hz_coeff = 2.*I.*rou./(pi.*u0.*yr.^3);
hz_impulse_jiexijie = hz_coeff.*((1+temp.^2).*exp(-temp.^2).*erf_L-L./(2.*r).*(1+(yr.^2)./2./r.^2).*erf_r + temp_L.*yr.^2./(pi).^0.5./r.^2.*exp(-temp_r.^2));
%% save data
save(['ground_finite_wire_source_jiexi_rou' num2str(rou) '_L' num2str(L) '_y' num2str(yr) '.mat'],'hz_impulse_jiexijie');
end