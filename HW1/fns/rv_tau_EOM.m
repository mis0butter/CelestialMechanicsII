function d_rv_t = rv_tau_EOM( rv_tau, a_VTN, mu, k, alpha ) 
% Integrate equations of motion for Cartesian state and tau velocity 

    r       = rv_tau(1:3) ; 
    v_tau   = rv_tau(4:6) ; 
    t       = rv_tau(end) ; 
    
    d_rv_t  = zeros(7,1) ; 
    d_rv_t(1:3) = v_tau ; 
    d_rv_t(4:6) = alpha * norm(v_tau) / norm(r) * v_tau ... 
                  - k^2 * mu / ( r^(3-2*alpha) ) * r ... 
                  + k^2 * norm(r)^(2*alpha) * a ; 
    d_rv_t(end) = k * r^(alpha) ; 

end 

