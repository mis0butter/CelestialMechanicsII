% case(9) ;
t1 = k ^ 2 ;
t2 = r ^ 2 ;
t3 = t2 * r ;
t4 = 0.1D1 / t3 ;
t8 = t1 * k ;
t9 = t2 ^ 2 ;
t10 = t9 * r ;
t11 = 0.1D1 / t10 ;
t17 = t1 ^ 2 ;
t18 = v2 + PdotR ;
t19 = t18 * t2 ;
t21 = RdotV ^ 2 ;
t22 = 0.5D1 * t21 ;
t25 = t9 * t3 ;
t26 = 0.1D1 / t25 ;
t34 = PdotV * t9 ;
t37 = 0.5D1 * RdotV * t18 * t2 ;
t38 = RdotV * r ;
t40 = t21 * RdotV ;
t41 = 0.1166666666666667D2 * t40 ;
t43 = t17 * k ;
t45 = t9 ^ 2 ;
t7 = 0.1D1 / t45 ;
t12 = 0.1D1 / r ;
t47 = t7 * t12 ;
t58 = t9 * t2 ;
t59 = t58 * p2 ;
t60 = RdotV * PdotV ;
t62 = t18 ^ 2 ;
t65 = (-0.20D2 * t60 - 0.5D1 * t62) * t9 ;
t73 = (0.70D2 * v2 + 0.70D2 * PdotR) * t21 ;
t76 = t21 * r ;
t78 = t21 ^ 2 ;
t79 = 0.105D3 * t78 ;
t81 = t17 * t1 ;
t84 = t7 * t4 ;
t97 = p2 * RdotV ;
t99 = PdotV * t18 ;
t101 = (0.5000000000000000D0 * t97 + t99) * t58 ;
t102 = PdotV * t10 ;
t108 = 0.3500000000000000D1 * RdotV * (0.2D1 * t60 + t62) * t9 ;
t114 = t18 * t21 ;
t119 = r * t40 ;
t121 = t78 * RdotV ;
t122 = 0.2310000000000000D2 * t121 ;
t124 = t17 * t8 ;
t127 = t7 * t11 ;
t146 = PdotV ^ 2 ;
t148 = p2 * v2 ;
t149 = p2 * PdotR ;
t151 = (0.2D1 * t146 + t148 + t149) * t45 ;
t152 = t25 * p2 ;
t158 = t62 * t18 ;
t161 = (-0.7D1 * t21 * p2 - 0.28D2 * t99 * RdotV - 0.2333333333333333D1 * t158) * t58 ;
t169 = t40 * PdotV ;
t170 = 0.84D2 * t169 ;
t171 = t62 * t21 ;
t172 = 0.63D2 * t171 ;
t187 = r * t78 ;
t190 = 0.2002000000000000D3 * t78 * t21 ;
t192 = t17 ^ 2 ;
t195 = t7 * t26 ;
t207 = (0.21D2 * PdotR + 0.21D2 * v2) * t40 ;
t227 = t192 * k ;
t287 = t45 ^ 2 ;
FGHT(1) = 0.0D0 ;
FGHT(2) = k ;
FGHT(3) = 0.0D0 ;
FGHT(4) = k ;
FGHT(5) = -0.5000000000000000D0 * t1 * t4 ;
FGHT(6) = 0.0D0 ;
FGHT(7) = 0.5000000000000000D0 * t1 ;
FGHT(8) = 0.0D0 ;
FGHT(9) = 0.5000000000000000D0 * t8 * t11 * RdotV ;
FGHT(10) = -0.1666666666666667D0 * t8 * t4 ;
FGHT(11) = 0.0D0 ;
FGHT(12) = 0.0D0 ;
FGHT(13) = 0.1250000000000000D0 * t17 * (t19 - 0.6666666666666667D0 * r - t22) * t26 ;
FGHT(14) = 0.2500000000000000D0 * t17 * t11 * RdotV ;
FGHT(15) = -0.4166666666666667D-1 * t17 * t4 ;
FGHT(16) = 0.0D0 ;
FGHT(17) = 0.7500000000000000D-1 * (t34 - t37 + 0.3333333333333333D1 * t38 + t41) * t43 * t47 ;
FGHT(18) = 0.7500000000000000D-1 * t43 * (t19 - 0.8888888888888889D0 * r - t22) * t26 ;
FGHT(19) = 0.7500000000000000D-1 * t43 * t11 * RdotV ;
FGHT(20) = 0.0D0 ;
FGHT(21) = 0.1250000000000000D-1 * (t59 + t65 + (0.6333333333333333D1 * PdotR + 0.7333333333333333D1 * v2) * t3 + (-0.2444444444444444D1 + t73) * t2 - 0.4666666666666667D2 * t76 - t79) * t81 * t84 ;
FGHT(22) = 0.5000000000000000D-1 * t81 * (t34 - t37 + 0.4166666666666667D1 * t38 + t41) * t47 ;
FGHT(23) = 0.2500000000000000D-1 * (t19 - 0.9444444444444444D0* r - t22) * t81 * t26 ;
FGHT(24) = 0.0D0 ;
FGHT(25) = -0.8928571428571429D-1 * (t101 - 0.7400000000000000D0 * t102 - t108 + 0.5040000000000000D1 * (0.8888888888888889D0 * PdotR + v2) * RdotV * t3 + 0.21D2 * (-0.8000000000000000D-1 + t114) * RdotV * t2 - 0.14D2 * t119 - t122) * t124 * t127 ;
FGHT(26) = 0.8928571428571429D-2 * (t59 + t65 + (0.8800000000000000D1 * v2 + 0.7800000000000000D1 * PdotR) * t3 + (-0.3822222222222222D1 + t73) * t2 - 0.56D2 * t76 - t79) * t124 * t84 ;
FGHT(27) = 0.1785714285714286D-1 * (t34 - t37 + 0.4400000000000000D1 * t38 + t41) * t124 * t47 ;
FGHT(28) = 0.0D0 ;
FGHT(29) = -0.1674107142857143D-1 * (t151 - 0.7600000000000000D0 * t152 + t161 + (0.2032000000000000D2 * t60 + 0.5360000000000000D1 * (v2 + 0.8059701492537313D0 * PdotR) * t18) * t10 + (t170 + t172 - 0.3893333333333333D1 * v2 - 0.3093333333333333D1 * PdotR) * t9 + (0.8651851851851852D0 + (-0.8960000000000000D2 * v2 - 0.8120000000000000D2 * PdotR) * t21) * t3 - 0.231D3 * (-0.1292929292929293D0 + t114) * t21 * t2 + 0.154D3 * t187 + t190) * t192 * t195 ;
FGHT(30) = -0.6696428571428571D-1 * t192 * (t101 - 0.8800000000000000D0 * t102 - t108 + 0.5880000000000000D1 * (0.9047619047619048D0 * PdotR + v2) * RdotV * t3 + (t207 - 0.2426666666666667D1 * RdotV) * t2 - 0.1633333333333333D2 * t119 - t122) * t127 ;
FGHT(31) = 0.3348214285714286D-2 * t192 * (t59 + t65 + (0.9200000000000000D1 * v2 + 0.8200000000000000D1 * PdotR) * t3 + (-0.4207407407407407D1 + t73) * t2 - 0.5880000000000000D2 * t76 - t79) * t84 ;
FGHT(32) = 0.0D0 ;
FGHT(33) = -0.1302083333333333D-1 * t227 * (t45 * t2 * p2 * PdotV + ((-0.7D1 * t149 - 0.7D1 * t148 - 0.14D2 * t146) * RdotV - 0.7D1 * PdotV * t62) * t45 + (0.5200000000000000D1 * t97 + 0.1077714285714286D2 * PdotV * (0.9056203605514316D0 * PdotR + v2)) * t25 + (0.21D2 * t40 * p2 + 0.126D3 * t99 * t21 + 0.21D2 * t158 * RdotV - 0.3954285714285714D1 * PdotV) * t58 - 0.4725714285714286D2 * (0.1909310761789601D1 * t60 + (v2 + 0.8331318016928658D0 * PdotR) * t18) * RdotV * t10 - 0.231D3 * RdotV * (t169 + t171 - 0.1464440321583179D0 * v2 - 0.1207173778602350D0 * PdotR) * t9 + 0.3256000000000000D3 * RdotV * (-0.2308802308802309D-1 + (v2 + 0.9189189189189189D0 * PdotR) * t21) * t3 + 0.6006000000000000D3 * t40 * (-0.1807081807081807D0 + t114) * t2 - 0.4004000000000000D3 * r * t121 - 0.429D3 * t78 * t40) / t287 * t12 ;
FGHT(34) = -0.1302083333333333D-1 * (t151 - 0.8971428571428571D0 * t152 + t161 + (0.2354285714285714D2 * t60 + 0.6125714285714286D1 * (v2 + 0.8302238805970149D0 * PdotR) * t18) * t10 + (t170 + t172 - 0.5302857142857143D1 * v2 - 0.4365714285714286D1 * PdotR) * t9 + (0.1510264550264550D1 + (-0.1024000000000000D3 * v2 - 0.94D2 * PdotR) * t21) * t3 - 0.231D3 * (-0.1766233766233766D0 + t114) * t21 * t2 + 0.176D3 * t187 + t190) * t227 * t195 ;
FGHT(35) = -0.2604166666666667D-1 * (t101 - 0.9171428571428571D0 * t102 - t108 + 0.6125714285714286D1 * (0.9085820895522388D0 * PdotR + v2) * RdotV * t3 + (t207 - 0.2651428571428571D1 * RdotV) * t2 - 0.1706666666666667D2 * t119 - t122) * t227 * t127 ;
