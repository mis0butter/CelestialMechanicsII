function [ K_A_hist, K_B_hist ] = K_check( state, mu )
% Compute K at each step. State must contain u and uprime in columns 1-8 

    % extract 
    u_hist      = state(:, 1:4) ; 
    uprime_hist = state(:, 5:8) ; 
    
    % energy 
    K_A_hist = [] ; 
    K_B_hist = [] ; 
    for i = 1 : length(u_hist) 

        u           = u_hist(i,:) ; 
        uprime      = uprime_hist(i,:) ; 

        % energy 
        eps = 2 * dot( uprime,uprime ) / dot( u,u ) - mu / dot( u,u ) ; 

        K_A         = - eps/2 * dot(u,u) + dot(uprime, uprime) - mu/2 ;     
        K_A_hist    = [ K_A_hist ; K_A ] ;     

        K2          = eps - 2 * dot( uprime,uprime ) / dot( u,u ) + mu / dot( u,u ) ; 
        K_B         = - dot( u,u ) / 2 * K2 ; 
        K_B_hist    = [ K_B_hist ; K_B ] ; 

    end 

end 


