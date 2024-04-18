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
            
            stark_time_order_1 
            
        case (2) 
            
            stark_time_order_2 
            
        case (3) 
            
            stark_time_order_3 
            
        case (4) 
            
            stark_time_order_4 

        case (5) 
            
            stark_time_order_5 
            
        case (6) 
            
            stark_time_order_6 
            
        case (7) 
            
            stark_time_order_7 

        case (8) 
            
            stark_time_order_8 

        case (9) 
            
            stark_time_order_9 

        case (10) 
            
            stark_time_order_10 

        case (11) 
            
            stark_time_order_11 

        case (12) 
            
            stark_time_order_12 
            
%         case (13) 
%             
%             stark_time_order_13 
%             
%         case (14) 
%             
%             stark_time_order_14 
%             
%         case (15) 
%             
%             stark_time_order_15 
%             
%         case (16) 
%             
%             stark_time_order_16 
%             
%         case (17) 
%             
%             stark_time_order_17 
%             
%         case (18) 
%             
%             stark_time_order_18 
            
        otherwise 
            
            disp('STOPPING: invalid order in coeffs file') ; 
            disp('max order available is 18') ; 
          	disp('requested order is ', order, stop) ; 
      
    end 

end 