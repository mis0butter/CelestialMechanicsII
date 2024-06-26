function [ tau_p, alphaSun ] = taup_caseSun( caseSun, a, e, mu, cSun ) 

    % tau period 
    n = sqrt( mu / a^3 ) ; 

    % 0: alpha = 0, tau = time 
    if caseSun == 0 
        tau_p = 2 * pi / ( n * cSun ) ; 
        alphaSun = 0 ; 

    % 1: alpha = 1, tau = eccentric anomaly 
    elseif caseSun == 1 
        tau_p = 2 * pi / ( n * cSun * a ) ;  
        alphaSun = 1 ; 

    % 2: alpha = 2, tau = true anomaly 
    elseif caseSun == 2 
        tau_p = 2 * pi / ( n * cSun * sqrt( a * ( 1 - e^2 ) ) ) ; 
        alphaSun = 2 ; 

    % 3: alpha = 3/2, tau = intermediate anomaly 
    elseif caseSun == 3 
        M     = 2 * e / ( 1 + e ) ; 
        tau_p = 4 * ellipke( M ) / ( cSun * sqrt( mu * ( 1 + e ) ) ) ; 
        alphaSun = 3/2 ;     

    % 4: alpha = alphaSun, generic Sundman transformation: user provided alpha  
    else 
        tau_p = 2 * pi ; 

    end 

end 