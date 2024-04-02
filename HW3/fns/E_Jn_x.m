
function E = E_Jn_x( M, n, x, truncate ) 

    E = M ; 
    
    for n_i = 1 : n 
        
        bessel_term = Jn_x( n_i, n_i * x, truncate ) ; 
        term        = 2 * bessel_term / n_i * sin( n_i * M ) ; 
        
        E = E + term ; 
        
    end 

end 


