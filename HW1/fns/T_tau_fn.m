function [ T_time, T_tau ] = T_tau_fn( oe, k, mu, alpha )

    % get useful orbital elements 
    a = oe(1) ; 
    e = oe(2) ; 

    % period in time (s) 
    T_time = 2*pi*sqrt( a^3 / mu ) ; 
    
    if alpha == 0 
        T_tau   = 1/k * T_time ; 
    elseif alpha == 1 
        T_tau   = T_time / (k*a) ; 
    elseif alpha == 2 
        h       = sqrt( mu * a * (1-e^2) ) ; 
        n       = 2*pi / T_time ; 
        T_tau   = T_time * n / ( h * k ) ; 
    elseif alpha == 3/2 
        Q       = ellipke( 2*e / (1+e), 1e-18 ) ; 
        P       = 4*Q / ( k * sqrt( mu * (1+e) ) ) ; 
        T_tau   = 1/k * P ; 
    else
        disp("alpha =/= 0, 1, 2, or 3/2. Exiting"); 
        return 
    end 

end 