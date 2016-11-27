%% compute the survey area 
%  along the y-line (0,y,0) to get the max offset
%  far field early time the Hz and Uz formula
%% clear screen and all parameters
clc;close all;clear all;
%% define the parameters
L = 2500;
I = 40;
rho_1 = logspace(-3,3,70);
sigma_1 = 1./rho_1;
sigma_v = [0.1 1 10 ].*1e-9;
%% Formula
% $$ H^{E}_{Z}{(t)} = \frac{3tIdl } {2\pi{r^4} }{\rho_1}sin(\phi)$$ 
%%
% $$ U^{E}_{Z}{(t)} = \frac{3Idl } {2\pi{r^4} }{\rho_1}sin(\phi)$$
%%
% $$ r = (\frac{3\rho_1{Idl} }{2\pi{\sigma_v} })^{\frac{1}{4}}$$   
%%
max_r = (3.*I.*L./(2.*pi.*(sigma_1'*sigma_v)) ).^(1/4);
%% Display the survey area
figure;
loglog(rho_1,max_r,'linewidth' ,2);
grid on;
title(['P=' num2str(I*L) 'Am,max-offset vs NoiseLevel and Resistivity']);
legend(['Noise level' num2str(sigma_v(1)*1e9) 'nV/m^2'],['Noise level' num2str(sigma_v(2)*1e9) 'nV/m^2'],...
    ['Noise level' num2str(sigma_v(3)*1e9) 'nV/m^2']);
legend boxoff;
xlabel('apparent resistivity/(Ohm*m)');
ylabel('max offset along y axis / (m)');






