function FGHT = starkCoeffs_ecc(X, acc, order, cSun, FGHT) 

    r  = sqrt( X(1)*X(1) + X(2)*X(2) + X(3)*X(3) ) ; 
    v2 = X(4)*X(4) + X(5)*X(5) + X(6)*X(6) ; 
    RdotV = dot(X(1:3),X(4:6)) ; 
    PdotR = dot(X(1:3),acc) ; 
    PdotV = dot(X(4:6),acc) ; 
    p2 = acc(1)*acc(1) + acc(2)*acc(2) + acc(3)*acc(3) ; 
    k  = cSun ; 
    
    switch (order) 

    case (0)

    case (1) 

        stark_ecc_order_1 

    case (2) 

        stark_ecc_order_2 

    case (3) 

        stark_ecc_order_3 

    case (4) 

        stark_ecc_order_4 

    case (5) 

        stark_ecc_order_5 

    case (6) 

        stark_ecc_order_6 

    case (7) 

        stark_ecc_order_7 

    case (8) 

        stark_ecc_order_8 

    case (9) 

        stark_ecc_order_9 

    case (10) 

        stark_ecc_order_10 

    case (11) 

        stark_ecc_order_11 

    case (12) 

        stark_ecc_order_12 

    case (13) 

        stark_ecc_order_13 

    case (14) 

        stark_ecc_order_14 

    case (15) 

        stark_ecc_order_15 

    case (16) 

        stark_ecc_order_16 

    case (17) 

        stark_ecc_order_17 

    case (18) 

        stark_ecc_order_18 
        
    end 

end 