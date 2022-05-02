% Xinyi Guo, Marion Dierickx, Philip Mocz
% April 2013
% Bayesian analysis of planet mass
clear all
close all
clc
format long
rng(47)

%% User inputs

filename = 'output/HD159868/p2_results.mat';
M_star = 1.087; % solar masses


%% Physical constants

G = 6.67384*10^(-11);
M_sun = 1.9891*10^30;
day_to_s = 86400;
parsec_to_m = 3.08567758*10^16;
au_to_m = 149597870700;
mass_of_jupiter = 1.898*10^27;


%% load and reformat posteriors

load(filename);

N_planets = length(posterior2);
V_posterior2 = posterior2{1}(1, :);
s_posterior2 = posterior2{1}(7, :);
for p = 1:N_planets;
    K_posterior2{p} = posterior2{p}(2, :);
    w_posterior2{p} = posterior2{p}(3, :);
    e_posterior2{p} = posterior2{p}(4, :);
    P_posterior2{p} = posterior2{p}(5, :);
    chi_posterior2{p} = posterior2{p}(6, :);
end


%% calculate planet masses and semi-major-axes (Feroz 2011)

for p = 1:N_planets;
    astar_sin_i{p} = K_posterior2{p} .* P_posterior2{p}*day_to_s .* ...
        sqrt(1-e_posterior2{p}.^2) / (2*pi) / au_to_m;
    m_sin_i{p} = K_posterior2{p} * (M_star*M_sun)^(2/3) .* ...
        (P_posterior2{p}*day_to_s).^(1/3) .* ...
        sqrt(1-e_posterior2{p}.^2) / (2*pi*G)^(1/3) / mass_of_jupiter;
    a{p} = (M_star*M_sun) * astar_sin_i{p} ./ (m_sin_i{p}*mass_of_jupiter);
end


%% plot posteriors

N_bins = 20;
fig_hdl = figure(3);
set(fig_hdl,'position',[100,100,1600,800]);
set(fig_hdl,'PaperPosition',[100,100,1600,800]);

acc = 1;
for p = 1:N_planets;
    subplot(N_planets, 3, acc);     acc = acc + 1;
    [astar_sin_i_MAP(p) astar_sin_i_err(p)] = plot_posterior( astar_sin_i{p}, N_bins, [ 'a_*\sin(i_' int2str(p) ')/{\rm au}' ], 0 )
    
    subplot(N_planets, 3, acc);     acc = acc + 1;
    [m_sin_i_MAP(p) m_sin_i_err(p)] = plot_posterior( m_sin_i{p}, N_bins, [ 'm_' int2str(p) '\sin(i_' int2str(p) ')/{\rm M_J}' ], 2*pi )
    
    subplot(N_planets, 3, acc);    acc = acc + 1;
    [a_MAP(p) a_err(p)] = plot_posterior( a{p}, N_bins, ['a_' int2str(p) '/{\rm au}' ], 0 )
end
