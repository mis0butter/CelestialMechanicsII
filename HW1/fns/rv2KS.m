function KS = rv2KS( rv ) 

    if ~iscolumn( rv ) 
        rv = rv' ; 
    end 

    % get position and velocity 
    r = rv(1:3) ;   x = r(1) ; y = r(2) ; z = r(3) ; 
    v = rv(4:6) ; 
    
    % deal with redundant variable 
    if r(1) >= 0 
        u1 = sqrt( 1/2 * ( x + norm(r) ) ) ; 
        u2 = y / ( 2 * u1 ) ; 
        u3 = z / ( 2 * u1 ) ; 
        u4 = 0 ;         
    else 
        u2 = sqrt( 1/2 * ( norm(r) - x ) ) ; 
        u1 = y / ( 2 * u2 ) ; 
        u3 = 0 ; 
        u4 = z / ( 2 * u2 ) ; 
    end 
    
    % KS parameters 
    u       = [ u1 ; u2 ; u3 ; u4 ] ; 
    uprime  = 1/2 * Lu_fn(u)' * [ v ; 0 ] ; 
    KS      = [ u ; uprime ] ; 

end 