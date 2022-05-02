% Xinyi Guo, Marion Dierickx, Philip Mocz
% April 2013
% Bayesian analysis of RV curves using nested sampling mcmc
clear all
close all
clc
format long
rng(8)

%% User inputs
save_parameters = 0; 
save_figs = 0;
filename = 'output/HD159868/';

N_planets_model = 3;

rv_data = load('data/HD159868.txt');
t = rv_data( :, 1);
t = t - min(t);
rv_data( :, 1) = t; 
sigma = rv_data( : ,3); 

% plot
figure(1);
errorbar(t,rv_data(:,2)',sigma,'r.','linewidth',1,'markersize',16)
%hold on
%plot(tt,rv,'b','linewidth',2)
%hold off
axis( [0 max(t)*1.1 min(rv_data(:,2))-5 max(rv_data(:,2))+25] )
xlabel('time (d)','interpreter','latex','fontsize',12)
ylabel('radial velocity (m s$^{-1}$)','interpreter','latex','fontsize',12)
title('Real Data and MAP model','interpreter','latex','fontsize',14)

%% Nested Sampling

N_planets = N_planets_model;

V_bounds = [-100 100];   [-2000 2000];
K_bounds = [0 200];      [0 2000];
w_bounds = [0 2*pi];
e_bounds = [0 0.8]; [0 1];
P_bounds = [0.2 2000];    [0.2 15*10^3];
chi_bounds = [0 1];
s_bounds = [0 2];        [0 2000];

prior_bounds = [V_bounds;
    K_bounds;
    w_bounds;
    e_bounds;
    P_bounds;
    chi_bounds;
    s_bounds];

N_alive = 250;
N_mcmc = 20;
N_posterior_sample = 2*10^5;
log_tolerance = -Inf;
max_iter = 20000;
live_plot = 1;

[ posterior, posterior2, logZ, L, W ] = nested_sample( rv_data, N_planets, ...
    prior_bounds, N_alive, N_mcmc, N_posterior_sample, log_tolerance, max_iter, live_plot, save_figs, filename );

%% read out the posteriors

V_posterior2 = posterior2{1}(1, :);
s_posterior2 = posterior2{1}(7, :);
for p = 1:N_planets;
    K_posterior2{p} = posterior2{p}(2, :);
    w_posterior2{p} = posterior2{p}(3, :);
    e_posterior2{p} = posterior2{p}(4, :);
    P_posterior2{p} = posterior2{p}(5, :);
    chi_posterior2{p} = posterior2{p}(6, :);
end


%% plot posteriors

N_bins = 20;
fig_hdl = figure(3);
set(fig_hdl,'position',[100,100,1600,800]);

% plot the marginalized posterior for each planet
acc = 1;
for p = 1:N_planets;
    subplot(N_planets + 1, 5, acc);     acc = acc + 1;
    [K_MAP(p) K_err(p)] = plot_posterior( K_posterior2{p}, N_bins, ['K_' int2str(p)], 0 )
    
    subplot(N_planets + 1, 5, acc);     acc = acc + 1;
    [w_MAP(p) w_err(p)] = plot_posterior( w_posterior2{p}, N_bins, ['w_' int2str(p)], 2*pi )
    
    subplot(N_planets + 1, 5, acc);     acc = acc + 1;
    [e_MAP(p) e_err(p)] = plot_posterior( e_posterior2{p}, N_bins, ['e_' int2str(p)], 0 )
    
    subplot(N_planets + 1, 5, acc);     acc = acc + 1;
    [P_MAP(p) P_err(p)] = plot_posterior( P_posterior2{p}, N_bins, ['P_' int2str(p)], 0 )
    
    subplot(N_planets + 1, 5, acc);     acc = acc + 1;
    [chi_MAP(p) chi_err(p)] = plot_posterior( chi_posterior2{p}, N_bins, ['\chi_' int2str(p)], 1 )
end

subplot(N_planets + 1, 5, acc)
[V_MAP V_err] = plot_posterior( V_posterior2, N_bins, 'V', 0 )
acc = acc + 1;
subplot(N_planets + 1, 5, acc)
[s_MAP s_err] = plot_posterior( s_posterior2, N_bins, 's', 0 )

if save_figs
    saveas(gcf,[filename 'posteriors.eps'],'psc2');
end


%% plot MAP model

figure(1)
hold on
p = 1;
tt = linspace(0,max(t)*1.1,1000);
rv_MAP = rv_model( V_MAP, K_MAP, w_MAP, e_MAP, P_MAP, chi_MAP, tt );
plot(tt,rv_MAP,'b-','linewidth',2)
legend_hdl = legend('data','MAP model');
set(legend_hdl,'interpreter','latex');
hold off
if save_figs
    saveas(gcf,[filename 'data_and_model.eps'],'psc2');
end


%% diagnostics -- plot convergence of Z

figure(4);
logZ_partial = 0*L;
for i = 1:length(W)
    L_partial = L(1:i);
    W_partial = W(1:i);
    logZ_partial(i) = log(exp(L_partial-max(L_partial))*W_partial) + max(L_partial);
end
semilogy(logZ_partial,'b','linewidth',2)
xlabel('\# nested sampling steps','interpreter','latex','fontsize',12)
ylabel('$\log Z$','interpreter','latex','fontsize',12)
title('Convergence of Bayesian Evidence','interpreter','latex','fontsize',14)
if save_figs
    saveas(gcf,[filename 'Z_convergence.eps'],'psc2');
end

%% save the results of the simulation
if save_parameters
    savefile = [ filename 'p' int2str(N_planets_model) '_results.mat' ];
    save( savefile, 'posterior2', 'rv_data', 'V_MAP', 'K_MAP', 'w_MAP', 'e_MAP', 'P_MAP', 'chi_MAP', 's_MAP', 'V_err', 'K_err', 'w_err', 'e_err', 'P_err', 'chi_err', 's_err', 'logZ');
end
