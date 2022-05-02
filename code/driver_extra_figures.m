% Xinyi Guo, Marion Dierickx, Philip Mocz
% April 2013
% generate extra post main run figures
clear all
close all
clc
format long
rng(32)

%% User inputs

filename = 'output/d2_';


%% plot rv data and all the best fit models

% color scheme
my_colors = lines(16);

% 1 planet fit
N_planets_data = 2;
N_planets_model = 1;
savefile = 'output/d2m1/results.mat';
load(savefile)

t = rv_data(:,1)';
sigma = rv_data(:,3)';

tt = linspace(0,max(t)*1.1,1000);

p = 1;
rv = rv_model( V, K(p), w(p), e(p), P(p), chi(p), tt);
for p = 2: N_planets_data
    rv = rv + rv_model( 0, K(p), w(p), e(p), P(p), chi(p), tt );
end
rv_MAP = rv_model( V, K_MAP, w_MAP, e_MAP, P_MAP, chi_MAP, tt);

% plot
figure(1);
errorbar(t,rv_data(:,2)',sigma,'r.','linewidth',1,'markersize',16)
hold on
plot(tt,rv,'k','linewidth',2)
plot(tt,rv_MAP,'--','color',my_colors(1,:),'linewidth',2)
hold off
axis( [0 max(t)*1.1 min(rv_data(:,2))-20 max(rv_data(:,2))+5+70] )
xlabel('time (d)','interpreter','latex','fontsize',12)
ylabel('radial velocity (m s$^{-1}$)','interpreter','latex','fontsize',12)
title('Synthetic Observations and Model Fit','interpreter','latex','fontsize',14)

% 2 planet fit
N_planets_model = 2;
savefile = 'output/d2m2/results.mat';
load(savefile)

tt = linspace(0,max(t)*1.1,1000);
rv_MAP = rv_model( V, K_MAP, w_MAP, e_MAP, P_MAP, chi_MAP, tt);
% plot
hold on
plot(tt,rv_MAP,'--','color',my_colors(2,:),'linewidth',2)
hold off


% 3 planet fit
N_planets_model = 3;
savefile = 'output/d2m3/results.mat';
load(savefile)

tt = linspace(0,max(t)*1.1,1000);
rv_MAP = rv_model( V, K_MAP, w_MAP, e_MAP, P_MAP, chi_MAP, tt);
% plot
hold on
plot(tt,rv_MAP,'--','color',my_colors(5,:),'linewidth',2)
hold off


% plot data
hold on
errorbar(t,rv_data(:,2)',sigma,'r.','linewidth',1,'markersize',16)
hold off

% legend
legend_hdl = legend('data','true model','1 planet fit','2 planet fit','3 planet fit');
set(legend_hdl,'interpreter','latex');
saveas(gcf,[filename 'data_and_model.eps'],'psc2');



