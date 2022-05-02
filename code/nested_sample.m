function [ posterior, posterior2, logZ, L, W ] = nested_sample( rv_data, N_planets, ...
    prior_bounds, N_alive, N_mcmc, N_posterior_sample, log_tolerance, max_iter, live_plot, save_figs, filename)
%NESTED_SAMPLE nested sampling algorithm for RV data with multiple planets

% rv_data                radial velocity data
% N_planets              number of planets to fit
% prior_bounds           bounds on priors
% N_alive                number of alive particles in nested sampling
% N_mcmc                 number of MCMC steps to take when sampling prior
% N_posterior_sample     size of posterior sample after resampling
% log_tolerance          tolerance limit for delta_logZi (best to set to -Inf)
% max_iter               maximum iterations to take
% live_plot              plot the nested sampling in real time
% save_figs              save figures (and movie) flag
% filename               beginning part of filename

% posterior              unweighted posterior distributions
% posterior2             resampled posterior with weights
% logZ                   Bayesian evidence
% L                      Likelihoods of unweighted posteriors
% W                      weights of unweighted posteriors

t = rv_data(:, 1)';
sigma = rv_data(:, 3)';
L = [];

% (1) draw sample from prior

% prior domain ranges
V_bounds = prior_bounds(1, :);
K_bounds = prior_bounds(2, :);
w_bounds = prior_bounds(3, :);
e_bounds = prior_bounds(4, :);
P_bounds = prior_bounds(5, :);
chi_bounds = prior_bounds(6, :);
s_bounds = prior_bounds(7, :);

k0 = 1;
s0 = 1;

e_mu = 0; 
e_sig = 0.3;
s_mu = 2.3; 
s_sig = 0.1; 

% functions
prior_V = @(V)     1 / (V_bounds(2) - V_bounds(1));
prior_K = @(K)     1 / (K_bounds(2) - K_bounds(1)); %(K + k0).^(-1) / log((k0 + K_bounds(2)) / k0);
prior_w = @(w)     1 / (w_bounds(2) - w_bounds(1));
prior_e = @(e)     1 / (e_bounds(2) - e_bounds(1));
prior_P = @(P)     1 / (P_bounds(2) - P_bounds(1)); %P.^(-1) / log(P_bounds(2) / P_bounds(1));
prior_chi = @(chi) 1 / (chi_bounds(2) - chi_bounds(1));
prior_s = @(s)     (s + s0).^(-1) /  log(( s0 + s_bounds(2) ) / s0);

prior = @(V, K, w, e, P, chi, s) prior_V(V).*prod(prior_K(K).*prior_w(w).*prior_e(e).*prior_P(P).*prior_chi(chi)).*prior_s(s);
log_likelihood = @(rv_data, rv_predicted, s, sigma) sum( log(1./sqrt(2*pi*(sigma.^2+s^2))) + (-(rv_predicted-rv_data).^2./(2*(sigma.^2+s^2))) );

% draw alive sample from priors
V_alive = sample_prior(1, N_alive, V_bounds, 'uniform');
K_alive = sample_prior(N_planets, N_alive, K_bounds, 'uniform'); % modjeffreys
w_alive = sample_prior(N_planets, N_alive, w_bounds, 'uniform');
e_alive = sample_prior(N_planets, N_alive, e_bounds, 'uniform');
% e_alive = sample_prior(N_planets, N_alive, e_bounds, 'halfnormal');
P_alive = sample_prior(N_planets, N_alive, P_bounds, 'uniform'); % jeffreys
chi_alive = sample_prior(N_planets, N_alive, chi_bounds, 'uniform');
s_alive = sample_prior(1, N_alive, s_bounds, 'modjeffreys');
% s_alive = random('Normal', s_mu, s_sig, [1, N_alive]); 

% calculate likelihoods of alive sample
L_alive = zeros(1,N_alive);
for j = 1:N_alive
    % initialize the rv_model with 1 planet fit
    rv_i = rv_model( V_alive(j), K_alive(:, j), w_alive(:, j), e_alive(:, j), P_alive(:, j), chi_alive(:, j), t );
    L_alive(j) = log_likelihood( rv_data(:,2)', rv_i, s_alive(j), sigma );
end

% proposal distribution Gaussian widths
sigma_V = 0.01*(V_bounds(2) - V_bounds(1));
sigma_K = 0.01*(K_bounds(2) - K_bounds(1));
sigma_w = 0.01*(w_bounds(2) - w_bounds(1));
sigma_e = 0.01*(e_bounds(2) - e_bounds(1));
sigma_P = 0.01*(P_bounds(2) - P_bounds(1));
sigma_chi = 0.01*(chi_bounds(2) - chi_bounds(1));
sigma_s = 0.01*(s_bounds(2) - s_bounds(1));

% initialize movie file
% movie
if save_figs
    writerObj = VideoWriter([filename 'nested_sampling.avi'],'Motion JPEG AVI');
    open(writerObj);
end

% main loop
i = 1;
delta_logZi = 1;
last_logZ = -Inf;
while (delta_logZi > log_tolerance) && (i <= max_iter)
    % (2) sort likelihoods, find the smallest one
    [L_alive, indices] = sort(L_alive);
    V_alive = V_alive(indices);
    K_alive = K_alive(:, indices);
    w_alive = w_alive(:, indices);
    e_alive = e_alive(:, indices);
    P_alive = P_alive(:, indices);
    chi_alive = chi_alive(:, indices);
    s_alive = s_alive(indices);
    
    % store the smallest likelihood point
    L(i) = L_alive(1);
    V_posterior(i) = V_alive(1);
    K_posterior(:, i) = K_alive(:, 1);
    w_posterior(:, i) = w_alive(:, 1);
    e_posterior(:, i) = e_alive(:, 1);
    P_posterior(:, i) = P_alive(:, 1);
    chi_posterior(:, i) = chi_alive(:, 1);
    s_posterior(i) = s_alive(1);
    delta_logZi = L_alive(N_alive) - i/N_alive;
    
    % (3) replace the point with a higher likelihood sampled point, using Metropolis-Hastings MCMC
    rnd_index = randi(N_alive-1)+1;
    V_new = V_alive(rnd_index);
    K_new = K_alive(:, rnd_index);
    w_new = w_alive(:, rnd_index);
    e_new = e_alive(:, rnd_index);
    P_new = P_alive(:, rnd_index);
    chi_new = chi_alive(:, rnd_index);
    s_new = s_alive(rnd_index);
    L_new = L_alive(rnd_index);
    N_accept = 0;
    N_reject = 0;
    
    while N_accept < N_mcmc
        V_prop = V_bounds(1);
        while V_prop <= V_bounds(1) || V_prop >= V_bounds(2)
            V_prop = random('normal',V_new,sigma_V);
        end
        K_prop = K_bounds(1)*ones(size(K_new));
        
        while sum( +(K_prop <= K_bounds(1)) ) > 0 || sum( +(K_prop >= K_bounds(2)) ) > 0
            K_prop = random('normal',K_new,sigma_K);
        end
        w_prop = w_bounds(1)*ones(size(w_new));
        while sum( +(w_prop <= w_bounds(1)) ) > 0 || sum( +(w_prop >= w_bounds(2)) ) > 0
            w_prop = mod(random('normal',w_new,sigma_w),2*pi);
        end
        e_prop = e_bounds(1)*ones(size(e_new));
        while sum( +(e_prop <= e_bounds(1)) ) > 0 || sum( +(e_prop >= e_bounds(2)) ) > 0
            e_prop = random('normal',e_new,sigma_e);
        end
        P_prop = P_bounds(1)*ones(size(P_new));
        while sum( +(P_prop <= P_bounds(1)) ) > 0 || sum( +(P_prop >= P_bounds(2)) ) > 0
            P_prop = random('normal',P_new,sigma_P);
        end
        chi_prop = chi_bounds(1)*ones(size(chi_new));
        while sum( +(chi_prop <= chi_bounds(1)) ) > 0 || sum( +(chi_prop >= chi_bounds(2)) ) > 0
            chi_prop = mod(random('normal',chi_new,sigma_chi),1);
        end
        s_prop = s_bounds(1);
        while s_prop <= s_bounds(1) || s_prop >= s_bounds(2)
            s_prop = random('normal',s_new,sigma_s);
        end
        
        rv_prop = rv_model( V_prop, K_prop, w_prop, e_prop, P_prop, chi_prop, t );
        L_prop = log_likelihood( rv_data(:,2)', rv_prop, s_prop, sigma );
        
        prior_new = prior(V_new, K_new, w_new, e_new, P_new, chi_new, s_new);
        prior_prop = prior(V_prop, K_prop, w_prop, e_prop, P_prop, chi_prop, s_prop);
        A = prior_prop/prior_new;
        U = rand();
        
        N_reject = N_reject + 1;
        if L_prop > L(i)
            if U < A
                V_new = V_prop;
                K_new = K_prop;
                w_new = w_prop;
                e_new = e_prop;
                P_new = P_prop;
                chi_new = chi_prop;
                s_new = s_prop;
                L_new = L_prop;
                N_accept = N_accept + 1;
                N_reject = N_reject - 1;
            end
        end
        
    end
    
    
    V_alive(1) = V_new;
    K_alive(:, 1) = K_new;
    w_alive(:, 1) = w_new;
    e_alive(:, 1) = e_new;
    P_alive(:, 1) = P_new;
    chi_alive(:, 1) = chi_new;
    s_alive(1) = s_new;
    L_alive(1) = L_new;
    
    
    % update the proposal Gaussian widths (adaptive mcmc step sizes)
    if N_accept > N_reject
        if sigma_s / (s_bounds(2) - s_bounds(1)) < 0.1 % cap maximum growth
            sigma_V = sigma_V*exp(1/N_accept);
            sigma_K = sigma_K*exp(1/N_accept);
            sigma_w = sigma_w*exp(1/N_accept);
            sigma_e = sigma_e*exp(1/N_accept);
            sigma_P = sigma_P*exp(1/N_accept);
            sigma_chi = sigma_chi*exp(1/N_accept);
            sigma_s = sigma_s*exp(1/N_accept);
        end
    else
        if sigma_s / (s_bounds(2) - s_bounds(1)) > 10^-4 % cap minimum shrinkage
            sigma_V = sigma_V*exp(-1/N_reject);
            sigma_K = sigma_K*exp(-1/N_reject);
            sigma_w = sigma_w*exp(-1/N_reject);
            sigma_e = sigma_e*exp(-1/N_reject);
            sigma_P = sigma_P*exp(-1/N_reject);
            sigma_chi = sigma_chi*exp(-1/N_reject);
            sigma_s = sigma_s*exp(-1/N_reject);
        end
    end
    
    if live_plot
        % plot the sampling in real time
        tt = linspace(0,max(t)*1.1,1000);
        rv = rv_model( V_posterior(i), K_posterior(:, i), w_posterior(:, i), e_posterior(:, i), P_posterior(:, i), chi_posterior(:, i), tt );
        rv_new = rv_model( V_new, K_new, w_new, e_new, P_new, chi_new, tt );
        
        % plot
        fig_hdl = figure(2);
        clf
        errorbar(t,rv_data(:,2)',sigma,'r.','linewidth',1,'markersize',16)
        hold on
        for q = max(1,i-10):i
            rvq = rv_model( V_posterior(q), K_posterior(:, q), w_posterior(:, q), e_posterior(:, q), P_posterior(:, q), chi_posterior(:, q), tt );
            plot(tt,rvq,'g','linewidth',1);
        end
        plot(tt,rv,'b--','linewidth',2)
        plot(tt,rv_new,'b','linewidth',2)
        hold off
        axis( [0 max(t)*1.1 min(rv_data(:,2))-5 max(rv_data(:,2))+5] )
        xlabel('time (d)','interpreter','latex','fontsize',12)
        ylabel('radial velocity (m s$^{-1}$)','interpreter','latex','fontsize',12)
        title(['Nested Sampling: step $' num2str(i) '$'],'interpreter','latex','fontsize',14)
        drawnow
        if save_figs
            % add frame to movie
            frame = getframe(fig_hdl);
            writeVideo(writerObj,frame);
            % same every 100th frame as eps
            if mod(i,100) == 0
                saveas(gcf,[filename 'nested' num2str(i) '.eps'],'psc2');
            end
        end
    end
    
    % calculate evidence weights
    N_L = i;
    X = exp(-(1:N_L)/N_alive)';
    W = 0.5*( circshift(X,[0 1]) + circshift(X,[0 -1]) );
    W(1) = X(1);
    W(end) = X(end);
    % compute Bayesian evidence Z (trapezoid rule integration)
    logZ = log(exp(L-max(L))*W) + max(L);
    delta_logZi = logZ - last_logZ; 
    last_logZ = logZ; 

    % print info about step
    fprintf(['step: ' num2str(i) ', N_reject = ' num2str(N_reject) ', prop. size = ' num2str(sigma_s / (s_bounds(2) - s_bounds(1))) ', current_logZ = ' num2str(logZ) ', delta_logZ = ' num2str(delta_logZi) '\n']);
       
    % increment counter
    i = i+1;
end

if save_figs
    close(writerObj);
end

%% Post-process results of nested sampling to find Z and posteriors

% calculate evidence weights
N_L = length(L);
X = exp(-(1:N_L)/N_alive)';
W = 0.5*( circshift(X,[0 1]) + circshift(X,[0 -1]) );
W(1) = X(1);
W(end) = X(end);
% compute Bayesian evidence Z (trapezoid rule integration)
logZ = log(exp(L-max(L))*W) + max(L);
delta_logZi = logZ - last_logZ; 
last_logZ = logZ; 
% print info about step
fprintf(['step: ' num2str(i) ', current_logZ = ' num2str(logZ) ', delta_logZi = ' num2str(delta_logZi) '\n']);


% calculate posterior sample weights and cumulative distribution
W_sample = exp(L+log(W')-logZ);
W_sample_cumulative = W_sample;
for i = 1:N_L
    W_sample_cumulative(i) = sum(W_sample(1:i));
end
W_sample_cumulative(N_L) = 1.0;

% resample the posterior with the weights
for i = 1:N_posterior_sample
    ss = rand;
    for j = 1:N_L
        if (j == 1 && (ss <= W_sample_cumulative(j))) || ...
                (j > 1 && ss <= W_sample_cumulative(j) && ss > W_sample_cumulative(j-1))
            V_posterior2(i) = V_posterior(j);
            K_posterior2(:, i) = K_posterior(:, j);
            w_posterior2(:, i) = w_posterior(:, j);
            e_posterior2(:, i) = e_posterior(:, j);
            P_posterior2(:, i) = P_posterior(:, j);
            chi_posterior2(:, i) = chi_posterior(:, j);
            s_posterior2(i) = s_posterior(j);
        end
    end
end


%% concatenate the posteriors for the output

posterior = cell(N_planets);
posterior2 = cell(N_planets);

for p = 1: N_planets;
    posterior{p} = [V_posterior; ...
        K_posterior(p, :); ...
        w_posterior(p, :); ...
        e_posterior(p, :); ...
        P_posterior(p, :); ...
        chi_posterior(p, :); ...
        s_posterior];
    
    posterior2{p} = [V_posterior2; ...
        K_posterior2(p, :); ...
        w_posterior2(p, :); ...
        e_posterior2(p, :); ...
        P_posterior2(p, :); ...
        chi_posterior2(p, :); ...
        s_posterior2];
end

end

