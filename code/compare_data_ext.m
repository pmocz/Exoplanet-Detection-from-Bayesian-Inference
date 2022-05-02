% Xinyi Guo, Marion Dierickx, Philip Mocz
% April 2013
% plot extended RV data
clear all
close all
clc
format long
rng(8)

%%
rv_data = load('data/HD159868_ext.txt');

t = rv_data( :, 1);
t = t - min(t);
tmin_axis = -400;0;
tt = linspace(tmin_axis,max(t)*1.1,1000);
rv_data( :, 1) = t;
sigma = rv_data( : ,3);

fig_hdl = figure(5);
set(fig_hdl,'Position',[100 100 1200 600]);
set(fig_hdl,'PaperPosition',[100 100 1200 600]);
errorbar(t(t<1700),rv_data((t<1700),2)',sigma(t<1700),'r.','linewidth',1,'markersize',16)
hold on
errorbar(t(t>1700),rv_data((t>1700),2)',sigma(t>1700),'m.','linewidth',1,'markersize',16)
%plot(tt,rv,'b','linewidth',2)
%hold off
axis( [tmin_axis max(t)*1.1 min(rv_data(:,2))-5 max(rv_data(:,2))+35] )
xlabel('time (d)','interpreter','latex','fontsize',12)
ylabel('radial velocity (m s$^{-1}$)','interpreter','latex','fontsize',12)
title('Extended Data and MAP model obtained through old data','interpreter','latex','fontsize',14)

N_planets_data = 2;
N_planets_model = 2;

if N_planets_model == 1
    K_MAP = [ 43.09 ];
    w_MAP = [ 4.849 ];
    e_MAP = [ 0.6985 ];
    P_MAP = [ 981.85 ];
    chi_MAP = [ 0.683 ];
    V_MAP =   0.851;
end

if N_planets_model == 2
    K_MAP = [  39.953262948547476  23.664635014227457];
    w_MAP = [   2.621801334796591   4.676768877623713];
    e_MAP = [   0.096568795218105   0.037592416464826];
    P_MAP = [   1169.371061923345   354.233854902685];
    chi_MAP = [   0.175125658297301   0.149718855846739];
    V_MAP =   6.283733996508242;
end

rv_MAP = rv_model( V_MAP, K_MAP, w_MAP, e_MAP, P_MAP, chi_MAP, tt );
plot(tt,rv_MAP,'b-','linewidth',2)
%line([1700 1700], [min(rv_data(:,2))-5 max(rv_data(:,2))+25], 'LineStyle', '--' )
legend('old data', 'extended Data', 'MAP model from old data')

hold off
