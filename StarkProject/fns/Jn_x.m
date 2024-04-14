


function J_out = Jn_x( n, x, truncate ) 

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
        
        % check if exponent for next iter exceeds truncate power 
        if (2 * j + n) > truncate 
            break 
        end 
        
    end 

end 