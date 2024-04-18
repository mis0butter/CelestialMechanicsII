function FGHT = starkCoeffs_time(X, acc, order, cSun, FGHT) 

%         use maxOrder
%         implicit double precision (t)
% 
%         integer, intent(in) :: order
%         real*8, intent(in)  :: X(7), acc(3), cSun
% 
%         real*8  :: FGHT(MAXORD*4)
%         real*8  :: r, RdotV, PdotR, PdotV, v2, p2, k

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
            
            run starkCoeffs/time/order_1 
            
        case (2) 
            
            run starkCoeffs/time/order_2 
            
        case (3) 
            
            run starkCoeffs/time/order_3 
            
        case (4) 
            
            run starkCoeffs/time/order_4 

        case (5) 
            
            run starkCoeffs/time/order_5 
            
        case (6) 
            
            run starkCoeffs/time/order_6 
            
        case (7) 
            
            run starkCoeffs/time/order_7 

        case (8) 
            
            run starkCoeffs/time/order_8 

        case (9) 
            
            run starkCoeffs/time/order_9 

        case (10) 
            
            run starkCoeffs/time/order_10 

        case (11) 
            
            run starkCoeffs/time/order_11 

        case (12) 
            
            run starkCoeffs/time/order_12 
            
%         case (13) 
%             
%             run starkCoeffs/time/order_13 
%             
%         case (14) 
%             
%             run starkCoeffs/time/order_14 
%             
%         case (15) 
%             
%             run starkCoeffs/time/order_15 
%             
%         case (16) 
%             
%             run starkCoeffs/time/order_16 
%             
%         case (17) 
%             
%             run starkCoeffs/time/order_17 
%             
%         case (18) 
%             
%             run starkCoeffs/time/order_18 
            
        otherwise 
            
            disp('STOPPING: invalid order in coeffs file') ; 
            disp('max order available is 18') ; 
          	disp('requested order is ', order, stop) ; 
      
    end 

end 