

%         case(2)
      t1 = k * r ; 
      t2 = k ^ 2 ; 
      t7 = 0.5000000000000000D0 * t2 * RdotV ; 
      t8 = r ^ 2 ; 
      FGHT(1) = 0.0D0 ; 
      FGHT(2) = t1 ; 
      FGHT(3) = 0.0D0 ; 
      FGHT(4) = t1 ; 
      FGHT(5) = -0.5000000000000000D0 * t2 / r ; 
      FGHT(6) = t7 ; 
      FGHT(7) = 0.5000000000000000D0 * t2 * t8 ; 
      FGHT(8) = t7 ; 
      