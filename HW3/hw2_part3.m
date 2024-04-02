
syms x e M 

n = 10 ; 
truncate = 10 ; 

% Jn_x( n, e, truncate ) 

E_Jn_x_out = E_Jn_x( M, n, e, truncate ) 


E_Jn_x_out    = E_Jn_x( M, n, e, truncate ) 
E_besselj_out = E_besselj( M, n, e ) 

M = 1 ; 
e = 0.5 ; 
n = 3 ; 
truncate = n ; 

Jn_x( n, e, truncate ) - besselj( n, e ) 

E_Jn_x_out    = E_Jn_x( M, n, e, truncate ) 
E_besselj_out = E_besselj( M, n, e ) 

E_Jn_x_out - E_besselj_out 


%% part 3.b 

M   = 1 ; 
tol = 1e-15 ; 

% loop for all eccentricities 
n_hist         = [] ; 
E_Jn_x_hist    = [] ; 
E_besselj_hist = [] ; 

e_hist = [ 0 : 0.05 : 0.75 ] ; 
for e = e_hist 
    
    n = 1 ; 
    
    E_Jn_x_old    = 0 ; 
    E_Jn_x_out    = E_Jn_x( M, n, e, n ) ; 
    E_besselj_out = E_besselj( M, n, e ) ; 

%     while abs( E_Jn_x_out - E_besselj_out ) > tol 
    while abs( E_Jn_x_out - E_Jn_x_old ) > tol 
        
        n = n + 1 ; 
        
        E_Jn_x_old    = E_Jn_x_out ; 
        E_Jn_x_out    = E_Jn_x( M, n, e, n ) ; 
%         E_besselj_out = E_besselj( M, n, e ) ; 
        
    end 
    
    n_hist         = [ n_hist ; n ] ; 
    E_Jn_x_hist    = [ E_Jn_x_hist ; E_Jn_x_out ] ; 
%     E_besselj_hist = [ E_besselj_hist ; E_besselj_out ] ; 
    
end 


%% 3.b plot 

figure(1) ; clf ; hold on ; grid on ; 
    scatter( e_hist, n_hist, 'linewidth', 2 ) 
    plot( e_hist, n_hist ) 
    xlabel( 'Eccentricity' ) 
    ylabel( '# terms' ) 
    title( 'Eccentricity vs. Terms needed for 1e-15 tol' ) 


%% subfunctions 

function J_out = Jn_x( n, x, order ) 

    J_out = 0 ;  
    j     = 0 ; 
    
    while true 
        
        % term at each j 
        num  = (-1)^j ; 
        den  = factorial( j ) * factorial( j + n ) ; 
        term = num / den * ( x / 2 )^( 2 * j + n ) ; 
        
        % sum 
        J_out = J_out + term ; 
        
        % increment j 
        j = j + 1 ; 
        
        % check if exponent for next iter exceeds desired order  
        if (2 * j + n) > order 
            break 
        end 
        
    end 

end 


function E = E_Jn_x( M, n, x, truncate ) 

    E = M ; 
    
    for n_i = 1 : n 
        
        bessel_term = Jn_x( n_i, n_i * x, truncate ) ; 
        term        = 2 * bessel_term / n_i * sin( n_i * M ) ; 
        
        E = E + term ; 
        
    end 

end 