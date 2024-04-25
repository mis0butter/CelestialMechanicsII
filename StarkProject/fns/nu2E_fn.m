% "Convert true anomaly to eccentric anomaly"
function E = nu2E_fn( nu, e ) 

    % elliptic 
    if e < 1.0 
        E = 2 * atan( sqrt( (1-e)/(1+e) ) * tan(nu/2) ) ; 

    % hyperbolic 
    else 
        E = acosh( (e+cos(nu)) / (1+e*cos(nu)) ) ; 
        if nu < 0.0 || nu > pi 
            E = -E ; 
        end 
    end 

    if E < 0.0 
        E = 2*pi + E ; 
    end

end 