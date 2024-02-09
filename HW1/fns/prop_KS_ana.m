function [  ] = prop_KS_ana( state0,  ) 

tau_hist   = [ 0 : T_tau / N : T_tau ] ; 
state_hist = [] ; 
for tau = tau_hist  

    u_uprime_t = KS_time_fn( u0, uprime0, tau, eps ) ; 
    state_hist = [ state_hist ; u_uprime_t' ] ; 
    
end 

% get Cartesian hist 
rv_hist = [] ; 
for i = 1 : length(u_hist) 
    KS      = [ u_hist(i,:)' ; uprime_hist(i,:)' ] ; 
    rv      = KS2rv( KS )' ;    
    rv_hist = [ rv_hist ; rv ] ; 
end 
    
    
    