function rv = KS2rv( KS )

    if ~iscolumn( KS ) 
        KS = KS' ; 
    end 

    % KS in u basis 
    u       = KS(1:4) ;
    uprime  = KS(5:8) ; 
    
    % recover position 
    r = Lu_fn(u) * u ; 
    r = r(1:3) ; 

    % recover velocity 
    v = 2/norm(r) * Lu_fn(u) * uprime ; 
    v = v(1:3) ; 
    
    % Cartesian state 
    rv = [ r ; v ] ; 

end 