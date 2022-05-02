function sample = sample_prior(N_planets, N, bounds, flag)
% SAMPLE_PRIOR samples the prior distribution as determined by the flag
% parameter
% Flag options: 'uniform', 'jeffreys', 'modjeffreys', 'halfnormal'

% N_planets     number of planets
% N             sample size
% bounds        lower and upper range of probability distribution: [l u]
% flag          distribution to sample from 

% break in 'modjeffreys' distribution function is at k0 = 1

k0 = 1;

lowerbound = bounds(1);
upperbound = bounds(2);

switch flag 
    case 'uniform' 
        sample = random('uniform',lowerbound, upperbound, [N_planets,N]);
    case 'jeffreys'
        rmin = log(lowerbound);
        rmax = log(upperbound);
        unifsample = random('uniform',rmin, rmax, [N_planets,N]);
        sample = exp(unifsample);
    case 'modjeffreys'
        rmin = log(lowerbound + k0);
        rmax = log(upperbound + k0);
        unifsample = random('uniform',rmin, rmax, [N_planets,N]);
        sample = exp(unifsample) - k0;    
    case 'halfnormal'
        rmin = 1/2*erf(2.35702*lowerbound);
        rmax = 1/2*erf(2.35702*upperbound);
        unifsample = random('uniform', rmin, rmax, [N_planets, N]);
        sample = erfinv(2.*unifsample)./2.35702;
end
