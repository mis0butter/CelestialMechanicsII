function u_uprime_t = KS_time_fn( u0, uprime0, tau, eps ) 
% Take the KS variables, energy and tau as input and
% converts them to the KS variables and time at the specified tau

    % frequency? 
    w = sqrt( - eps / 2 ) ;  
    
    I4   = eye(4) ;  
    wtau = w * tau ; 
    C    = [    cos(wtau) * I4      sin(wtau) / w * I4 ; 
                -w*sin(wtau) * I4   cos(wtau) * I4 ] ; 
            
    % get u and uprime 
    u_uprime = C * [ u0 ; uprime0 ] ; 
    u        = u_uprime(1:4) ; 
    uprime   = u_uprime(5:8) ; 
    
    % time wtf ... ok simplifying this 
    u01 = u0(1) ;       u02 = u0(2) ;       u03 = u0(3) ;       u04 = u0(4) ; 
    up01 = uprime0(1) ; up02 = uprime0(2) ; up03 = uprime0(3) ; up04 = uprime0(4) ; 
    
    % color coded the notes 
    red    = u01 * up01 + u02 * up02 + u03 * up03 + u04 * up04 ; 
    blue   = -u01^2 - u02^2 - u03^2 - u04^2 ; 
    yellow = u01^2 + u02^2 + u03^2 + u04^2 ; 
    black1 = blue * w^2 + up04^2 + up01^2 + up02^2 + up03^2 ; 
    black2 = yellow * w^2 + up04^2 + up01^2 + up02^2 + up03^2 ; 
    
    % time 
    t = 1/( 2*w^3 ) * ( ... 
        -2*w * red * cos(wtau)^2 - sin(wtau) * black1 * cos(wtau) + w * black2 * tau ... 
        ) ;  
    
    % output 
    u_uprime_t = [ u; uprime; t ] ; 

end 



