function d_rvt_dtau = rv_tau_EOM( tau, r_v_tau_t, a, mu, k, alpha ) 
% Integrate equations of motion for Cartesian state and tau velocity 

    r       = r_v_tau_t(1:3) ; 
    v_tau   = r_v_tau_t(4:6) ; 
    t       = r_v_tau_t(end) ; 
    
    % r, v, t derivatives wrt tau 
    dr_dtau = v_tau ; 
    dv_dtau = alpha * dot( r, v_tau ) / norm(r)^2 * v_tau ... 
                  - k^2 * mu / ( norm(r)^(3-2*alpha) ) * r ... 
                  + k^2 * norm(r)^(2*alpha) * a ; 
    dt_dtau = k * norm(r)^alpha ; 
    
    % output 
    d_rvt_dtau = [ dr_dtau ; dv_dtau ; dt_dtau ] ; 

end 

