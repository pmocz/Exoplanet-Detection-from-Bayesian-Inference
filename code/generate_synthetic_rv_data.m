function rv = generate_synthetic_rv_data( V, K, w, e, P, chi, t, sigma, s )
%GENERATE_SYNTHETIC_RV_DATA generates synthetic observations of radial
%velocity curve

% V       systemic velocity
% K       velocity semiamplitude
% w       longitude of periastron
% e       eccentricity
% P       orbital period
% chi     fraction of an orbit, prior to t=0, at which periastron occured
% t       array of times to output radail velocities
% sigma   array of Gaussian measurement errors at each time
% s       stellar jitter Gaussian error

% units: V and K should have same units. P and t should have same units.

% returns: times in the first column, rv_data in the second column, sigma
% in the third column

rv = zeros(numel(t),3);

rv(:,1) = t;
rv(:,2) = rv_model( V, K, w, e, P, chi, t ) + normrnd(0,sqrt(sigma.^2 + s^2));
rv(:,3) = sigma;

end

