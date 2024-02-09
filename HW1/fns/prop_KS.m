function [state, rv_hist] = prop_KS( state0, a, mu, T_tau, N, tol ) 

    % Solve ODE 
    [tau_hist, state] = ode78rpr(@(tau, KS_t_eps) KS_EOM( tau, KS_t_eps, a, mu), 0, T_tau, T_tau / N, state0, tol) ; 

    % extract 
    u_hist      = state(:, 1:4) ; 
    uprime_hist = state(:, 5:8) ; 
    t_hist      = state(:, 9) ; 
    eps_hist    = state(:, 10) ; 

    % get Cartesian hist 
    rv_hist = [] ; 
    for i = 1 : length(u_hist) 
        KS      = [ u_hist(i,:)' ; uprime_hist(i,:)' ] ; 
        rv      = KS2rv( KS )' ;    
        rv_hist = [ rv_hist ; rv ] ; 
    end 

end 