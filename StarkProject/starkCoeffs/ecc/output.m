%         case(7) ;
      t1 = k * r ;
      t2 = k ^ 2 ;
      t3 = 0.1D1 / r ;
      t7 = 0.5000000000000000D0 * t2 * RdotV ;
      t8 = r ^ 2 ;
      t11 = t2 * k ;
      t12 = v2 + PdotR ;
      t13 = t12 * r ;
      t14 = -0.2D1 + t13 ;
      t23 = t2 ^ 2 ;
      t24 = 0.1D1 / t8 ;
      t28 = t8 * r ;
      t29 = t28 * PdotV ;
      t37 = 0.1250000000000000D0 * t23 * (t29 + 0.3333333333333333D0 * t12 * RdotV * r - 0.6666666666666667D0 * RdotV) * t3 ;
      t38 = RdotV ^ 2 ;
      t40 = v2 * t8 ;
      t43 = PdotR * t8 ;
      t48 = t23 * k ;
      t53 = t8 * PdotV * RdotV ;
      t54 = 0.9D1 * t53 ;
      t55 = t8 ^ 2 ;
      t56 = t55 * p2 ;
      t58 = PdotR * r ;
      t60 = v2 ^ 2 ;
      t61 = t8 * t60 ;
      t62 = v2 * r ;
      t63 = 0.4D1 * t62 ;
      t64 = t40 * PdotR ;
      t65 = 0.2D1 * t64 ;
      t66 = PdotR ^ 2 ;
      t67 = t8 * t66 ;
      t76 = PdotR * v2 ;
      t78 = RdotV * PdotV ;
      t92 = t23 * t2 ;
      t100 = p2 * RdotV ;
      t102 = PdotV * t12 ;
      t103 = 0.1200000000000000D1 * t100 + t102 ;
      t104 = t103 * t55 ;
      t107 = t12 ^ 2 ;
      t111 = 0.6666666666666667D-1 * (0.9D1 * t78 + t107) * RdotV * t8 ;
      t116 = 0.2666666666666667D0 * (v2 + 0.3250000000000000D1 * PdotR) * RdotV * r ;
      t117 = 0.2666666666666667D0 * RdotV ;
      t126 = t38 * r ;
      t139 = t28 * t66 ;
      t141 = t55 * r ;
      t142 = t141 * p2 ;
      t144 = -0.30D2 * t38 + 0.99D2 * t29 * RdotV + 0.15D2 * t126 * v2 + 0.15D2 * t126 * PdotR - 0.40D2 * t40 + 0.25D2 * r - 0.58D2 * t43 + 0.16D2 * t28 * t60 + 0.32D2 * t28 * v2 * PdotR + 0.16D2 * t139 + 0.18D2 * t142 ;
      t152 = t23 * t11 ;
      t171 = PdotV ^ 2 ;
      t174 = t66 * PdotR ;
      t176 = t60 * v2 ;
      t178 = t28 * RdotV ;
      t84 = t178 * PdotV ;
      t193 = -0.8D1 - 0.108D3 * t53 + 0.12D2 * t62 + 0.60D2 * t58 - 0.57D2 * t56 - 0.6D1 * t61 - 0.30D2 * t67 - 0.36D2 * t64 + 0.33D2 * t142 * v2 + 0.33D2 * t142 * PdotR + 0.45D2 * t141 * t171 + t28 * t174 + t28 * t176 + 0.54D2 * t84 * v2 + 0.54D2 * t84 * PdotR + 0.3D1 * t28 * PdotR * t60 + 0.3D1 * t139 * v2 + 0.45D2 * t28 * t38 * p2 ;
      FGHT(1) = 0.0D0 ;
      FGHT(2) = t1 ;
      FGHT(3) = 0.0D0 ;
      FGHT(4) = t1 ;
      FGHT(5) = -0.5000000000000000D0 * t2 * t3 ;
      FGHT(6) = t7 ;
      FGHT(7) = 0.5000000000000000D0 * t2 * t8 ;
      FGHT(8) = t7 ;
      FGHT(9) = 0.0D0 ;
      FGHT(10) = 0.1666666666666667D0 * t11 * t14 ;
      FGHT(11) = 0.5000000000000000D0 * t11 * r * RdotV ;
      FGHT(12) = 0.1666666666666667D0 * t11 * (-0.1D1 + t13) ;
      FGHT(13) = -0.4166666666666667D-1 * t23 * t24 * t14 ;
      FGHT(14) = t37 ;
      FGHT(15) = 0.4166666666666667D-1 * t23 * (0.3D1 * t38 + 0.4D1 * t40 - 0.5D1 * r + 0.4D1 * t43) ;
      FGHT(16) = t37 ;
      FGHT(17) = -0.5000000000000000D-1 * t48 * PdotV ;
      FGHT(18) = 0.8333333333333333D-2 * t48 * t3 * (t54 + 0.3D1 * t56 - 0.7D1 * t58 + t61 - t63 + t65 + t67 + 0.4D1) ;
      FGHT(19) = 0.1250000000000000D0 * t48 * ((-0.1400000000000000D1 + t13) * RdotV + t29) ;
      FGHT(20) = 0.2500000000000000D-1 * t48 * (0.6666666666666667D0 + t56 + (0.6666666666666667D0 * t76 + 0.3D1 * t78 + 0.3333333333333333D0 * t66 + 0.3333333333333333D0 * t60) * t8 + (-0.1D1 * v2 - 0.2D1 * PdotR) * r) * t3 ;
      FGHT(21) = -0.1388888888888889D-2 * t92 * (0.9D1 * t56 - 0.13D2 * t58 + t54 + t61 - t63 + t65 + t67 + 0.4D1) / t28 ;
      FGHT(22) = 0.2083333333333333D-1 * t92 * (t104 - 0.2D1 * t29 + t111 - t116 + t117) * t24 ;
      FGHT(23) = 0.1388888888888889D-2 * t92 * t3 * t144 ;
      FGHT(24) = 0.2083333333333333D-1 * t92 * (t104 - 0.1400000000000000D1 * t29 + t111 - t116 + t117) * t24 ;
      FGHT(25) = -0.5952380952380952D-2 * t152 * (t103 * r - 0.2D1 * PdotV) * t3 ;
      FGHT(26) = 0.1984126984126984D-3 * t152 * t24 * t193 ;
      FGHT(27) = 0.4166666666666667D-1 * t152 * ((0.9000000000000000D0 * t100 + t102) * t55 - 0.1271428571428571D1 * t29 + 0.3000000000000000D0 * (0.4D1 * t78 + t107) * RdotV * t8 - 0.9428571428571429D0 * (v2 + 0.1636363636363636D1 * PdotR) * RdotV * r + 0.6857142857142857D0 * RdotV) * t3 ;
      FGHT(28) = 0.1984126984126984D-3 * (t152 * (0.33D2 * p2 * v2 + 0.33D2 * p2 * PdotR + 0.45D2 * t171) * t141 - 0.39D2 * t152 * p2 * t55 + t152 * (t174 + 0.3D1 * t66 * v2 + (0.54D2 * t78 + 0.3D1 * t60) * PdotR + 0.45D2 * t38 * p2 + t176 + 0.54D2 * t78 * v2) * t28 + t152 * (-0.90D2 * t78 - 0.5D1 * t60 - 0.29D2 * t66 - 0.34D2 * t76) * t8 + t152 * (0.8D1 * v2 + 0.38D2 * PdotR) * r - 0.4D1 * t152) * t24 ;
      