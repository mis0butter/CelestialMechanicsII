
% case 3             

t1 = k^2;
t2 = r^2;
t5 = 0.1e1 / r;
t4 = 0.1e1 / t2 * t5;
t8 = t1 * k;
t9 = t2^2;
FGHT(1) = 0.0e0;
FGHT(2) = k;
FGHT(3) = 0.0e0;
FGHT(4) = k;
FGHT(5) = -0.5000000000000000e0 * t1 * t4;
FGHT(6) = 0.0e0;
FGHT(7) = 0.5000000000000000e0 * t1;
FGHT(8) = 0.0e0;
FGHT(9) = 0.5000000000000000e0 * t8 / t9 * t5 * RdotV;
FGHT(10) = -0.1666666666666667e0 * t8 * t4;
FGHT(11) = 0.0e0;
FGHT(12) = 0.0e0;