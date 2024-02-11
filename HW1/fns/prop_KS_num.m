function [state, rv_hist] = prop_KS_num( state0, a, mu, T_tau, N, tol ) 
% propagate fixed-step KS EOM. state contains: 
%   u, uprime, t, eps (energy) 

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
    
    % add to state 
    state(:,11) = tau_hist ; 
    
    % sometimes last state is propagated basically 0 tau, check for that 
    dtau_hist = diff(tau_hist) ; 
    if any( dtau_hist < 1e-10 ) 
        state = state(1:end-1,:) ; 
    end 

end 