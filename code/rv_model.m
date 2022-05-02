function rv = rv_model( V, K, w, e, P, chi, t )
%RV_MODEL generates radial velocity curve for given orbital parameters at
%the array of times specified for a single planet

% V       systemic velocity
% K       velocity semiamplitude
% w       longitude of periastron
% e       eccentricity
% P       orbital period
% chi     fraction of an orbit, prior to t=0, at which periastron occured
% t       array of times to output radail velocities

% units: V and K should have same units. P and t should have same units.

%% initialize

N = length(t);
N_planets = length(K);

%initialize RV curve
rv = 0*t + V;

for p = 1:N_planets
    %% calculate the true anomalies f at each of the times
    % mean anomalies M in [0,2pi]
    M = ( 2*pi/P(p) )*( t + chi(p)*P(p) );
    M = mod(M, 2*pi);
    % solve Kepler's equation with Newton-Raphson iterator to find eccentric
    % anomalies E
    E = 0*M;
    for i = 1:N
        tolerance = 10^-8;
        max_iter = 100;
        Ek = pi;
        % f = @(E) E - e*sin(E) - M(i);
        % f_prime = @(E) = 1 - e*cos(E);
        iter = 0;
        delta_E = 1;
        while iter <= max_iter && abs(delta_E) > tolerance
            delta_E = -(Ek - e(p)*sin(Ek) - M(i))/(1 - e(p)*cos(Ek));   % -f(Ek)/f_prime(Ek);
            Ek = Ek + delta_E;
            iter = iter + 1;
        end
        if abs(delta_E) > tolerance
            fprintf('Error: Newton-Raphson Iterative failed!');
        end
        E(i) = Ek;
    end
    % true anomalies (http://en.wikipedia.org/wiki/True_anomaly)
    f = 2*atan2( sqrt(1+e(p))*sin(E/2), sqrt(1-e(p))*cos(E/2) );
    cos_f = (cos(E) - e(p))./(1 - e(p)*cos(E));
    sin_f = (sqrt(1 - e(p)^2)*sin(E))./(1 - e(p)*cos(E));
    
    %% add effect of planet to the RV curve
    %rv = rv - K*(sin(f+w) + e*sin(w));
    rv = rv - K(p)*(cos(f+w(p)) + e(p)*cos(w(p)));
    %rv = rv - K*(sin_f*cos(w) - cos_f*sin(w) + e*sin(w));
    %rv = rv - K*(cos(w)*(e+cos_f) - sin(w)*sin_f);
end

end

