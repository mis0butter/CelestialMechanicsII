function d_rvt_dtau = rv_t_EOM( tau, rvt, a, mu, k, alpha ) 
% Integrate equations of motion for Cartesian state and time velocity 

    r = rvt(1:3) ; 
    v = rvt(4:6) ; 
    t = rvt(end) ; 
    
    % r, v, t derivatives wrt time 
    dr = v ; 
    dv = -mu/norm(r)^3 * r + a ; 
    dt = 1 ; 
    
    % Sundman transformation 
    d_rvt_dtau = k * norm(r)^alpha * [ dr ; dv ; dt ] ; 
    
end 

