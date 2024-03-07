
function E = E_besselj( M, n, x ) 

    E = M ; 
    
    for n_i = 1 : n 
        
        bessel_term = besselj( n_i, n_i * x ) ; 
        term        = 2 * bessel_term / n_i * sin( n_i * M ) ; 
        
        E = E + term ; 
        
    end 

end 


