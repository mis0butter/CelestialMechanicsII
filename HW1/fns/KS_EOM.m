function KS_t_eps_prime = KS_EOM( tau, KS_t_eps, a, mu ) 
% Integrate equations of motion for u, uprime, t, and energy 

    % get u and uprime (tau velocity)  
    u        = KS_t_eps(1:4) ; 
    u_prime  = KS_t_eps(5:8) ; 
    
    % a is in VTN frame. Transform into u basis ? 
        
    % get uprimeprime (tau acceleration) 
    r               = dot( u, u ) ; 
    eps             = 2 * dot( u_prime, u_prime ) / r - mu / r ;     
    u_prime_prime   = eps/2 * u + 1/2 * r * Lu_fn(u)' * a ; 
    
    % time and energy derivatives wrt tau 
    t_prime   = r ;     
    eps_prime = 2 * dot(Lu_fn(u)' * a, u_prime ) ; 

    % output 
    KS_t_eps_prime = [ u_prime ; u_prime_prime ; t_prime ; eps_prime ] ; 

end 