function FGHT = starkCoeffs_ecc(X, acc, order, cSun, FGHT) 

    r  = sqrt( X(1)*X(1) + X(2)*X(2) + X(3)*X(3) ) ; 
    v2 = X(4)*X(4) + X(5)*X(5) + X(6)*X(6) ; 
    RdotV = dot(X(1:3),X(4:6)) ; 
    PdotR = dot(X(1:3),acc) ; 
    PdotV = dot(X(4:6),acc) ; 
    p2 = acc(1)*acc(1) + acc(2)*acc(2) + acc(3)*acc(3) ; 
    k  = cSun ; 
    
    switch (order) 
        
    end 
    
    

end 