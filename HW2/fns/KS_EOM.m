function KS_t_eps_prime = KS_EOM( tau, KS_t_eps, a_VTN, mu ) 
% Integrate equations of motion for u, uprime, t, and energy 

    % get u and uprime (tau velocity)  
    u        = KS_t_eps(1:4) ; 
    u_prime  = KS_t_eps(5:8) ; 
    
    % a is in VTN frame. Transform into inertial, then u basis ? 
    % first get r, v, and a stuff ... 
    rv  = KS2rv( [ u; u_prime ] ) ; 
    r   = rv(1:3) ; v = rv(4:6) ; 
    
    % get vectors you need 
    r_hat = r / norm(r) ; v_hat = v / norm(v) ; 
    h     = cross(r, v) ; h_hat = h / norm(h) ; 
    T     = cross(h, v) ; T_hat = T / norm(T) ; 
    
    % a_VTN --> a_Cartesian 
    a_V   = a_VTN(1) ; a_T = a_VTN(2) ; a_N = a_VTN(3) ; 
    a_XYZ = a_V * v_hat + a_T * ( cross(h_hat, v_hat) ) + a_N * h_hat ; 
    a_XYZ = [ a_XYZ ; 0 ] ; 
    
    % get uprimeprime (tau acceleration) 
    r               = dot( u, u ) ; 
    eps             = 2 * dot( u_prime, u_prime ) / r - mu / r ;     
    u_prime_prime   = eps/2 * u + 1/2 * r * Lu_fn(u)' * a_XYZ ; 
    
    % time and energy derivatives wrt tau 
    t_prime   = r ;     
    eps_prime = 2 * dot(Lu_fn(u)' * a_XYZ, u_prime ) ; 

    % output 
    KS_t_eps_prime = [ u_prime ; u_prime_prime ; t_prime ; eps_prime ] ; 

end 