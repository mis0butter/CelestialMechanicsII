% case 12 

            
            t1 = k^2;
            t2 = r^2;
            t3 = t2 * r;
            t4 = 0.1 / t3;
            t8 = t1 * k;
            t9 = t2^2;
            t10 = t9 * r;
            t11 = 0.1 / t10;
            t17 = t1^2;
            t18 = v2 + PdotR;
            t19 = t18 * t2;
            t21 = RdotV^2;
            t22 = 0.5 * t21;
            t25 = t9 * t3;
            t26 = 0.1 / t25;
            t34 = PdotV * t9;
            t35 = RdotV * t18;
            t37 = 0.5 * t35 * t2;
            t38 = RdotV * r;
            t40 = t21 * RdotV;
            t41 = 0.1166666666666667 * t40;
            t43 = t17 * k;
            t45 = t9^2;
            t46 = t45 * r;
            t47 = 0.1 / t46;
            t58 = t9 * t2;
            t59 = t58 * p2;
            t60 = RdotV * PdotV;
            t62 = t18^2;
            t65 = (-0.20 * t60 - 0.5 * t62) * t9;
            t73 = (0.70 * v2 + 0.70 * PdotR) * t21;
            t76 = t21 * r;
            t78 = t21^2;
            t79 = 0.105 * t78;
            t81 = t17 * t1;
            t83 = t45 * t3;
            t84 = 0.1 / t83;
            t97 = p2 * RdotV;
            t99 = PdotV * t18;
            t101 = (0.5000000000000000 * t97 + t99) * t58;
            t102 = PdotV * t10;
            t108 = 0.3500000000000000 * RdotV * (0.2 * t60 + t62) * t9;
            t114 = t18 * t21;
            t119 = r * t40;
            t121 = t78 * RdotV;
            t122 = 0.2310000000000000 * t121;
            t124 = t17 * t8;
            t126 = t45 * t10;
            t127 = 0.1 / t126;
            t146 = PdotV^2;
            t148 = p2 * v2;
            t149 = p2 * PdotR;
            t151 = (0.2 * t146 + t148 + t149) * t45;
            t152 = t25 * p2;
            t154 = t21 * p2;
            t158 = t62 * t18;
            t161 = (-0.7 * t154 - 0.28 * t99 * RdotV - 0.2333333333333333 * t158) * t58;
            t169 = t40 * PdotV;
            t170 = 0.84 * t169;
            t171 = t62 * t21;
            t172 = 0.63 * t171;
            t187 = r * t78;
            t189 = t78 * t21;
            t190 = 0.2002000000000000 * t189;
            t192 = t17^2;
            t195 = 0.1 / t45 * t26;
            t207 = (0.21 * PdotR + 0.21 * v2) * t40;
            t227 = t192 * k;
            t228 = t45 * t2;
            t230 = t228 * p2 * PdotV;
            t236 = PdotV * t62;
            t239 = ((-0.7 * t149 - 0.7 * t148 - 0.14 * t146) * RdotV - 0.7 * t236) * t45;
            t247 = t40 * p2;
            t248 = 0.21 * t247;
            t250 = 0.126 * t99 * t21;
            t252 = 0.21 * t158 * RdotV;
            t281 = r * t121;
            t283 = t78 * t40;
            t284 = 0.429 * t283;
            t287 = t45^2;
            t64 = 0.1 / t287;
            t289 = t64 / r;
            t333 = t192 * t1;
            t334 = t45 * t9;
            t335 = p2^2;
            t336 = t334 * t335;
            t338 = p2 * PdotV * RdotV;
            t341 = t148 + t149 + 0.4 * t146;
            t342 = t18 * t341;
            t345 = (-0.56 * t338 - 0.14 * t342) * t228;
            t355 = (0.252 * t148 + 0.504 * t146 + 0.252 * t149) * t21;
            t357 = 0.504 * t236 * RdotV;
            t358 = PdotR^2;
            t359 = t358 * PdotR;
            t360 = v2 * t359;
            t361 = 0.84 * t360;
            t363 = v2^2;
            t364 = t363^2;
            t365 = 0.21 * t364;
            t366 = t358^2;
            t367 = 0.21 * t366;
            t368 = t363 * v2;
            t369 = t368 * PdotR;
            t370 = 0.84 * t369;
            t371 = t363 * t358;
            t372 = 0.126 * t371;
            t387 = t78 * p2;
            t388 = 0.462 * t387;
            t389 = t99 * t40;
            t390 = 0.3696 * t389;
            t391 = t158 * t21;
            t392 = 0.924 * t391;
            t396 = v2 * PdotR;
            t410 = t121 * PdotV;
            t411 = 0.4804800000000000 * t410;
            t412 = t62 * t78;
            t413 = 0.6006 * t412;
            t431 = r * t189;
            t433 = t78^2;
            t434 = 0.7293 * t433;
            t435 = t336 + t345 + (0.4331428571428571 * t146 + 0.2197714285714286 * t148 + 0.1996000000000000 * t149) * t46 + (t355 + t357 + t361 - 0.8257142857142857 * p2 + t365 + t367 + t370 + t372) * t45 + (-0.1842285714285714 * t154 - 0.7598857142857143 * (v2 + 0.9186343811099413 * PdotR) * PdotV * RdotV - 0.6592380952380952 * t62 * (v2 + 0.7741982086102283 * PdotR)) * t25 + (-t388 - t390 - t392 + 0.2744000000000000 * t60 + 0.7496000000000000 * t363 + 0.5241142857142857 * t358 + 0.1263542857142857 * t396) * t58 + (0.2621142857142857 * t169 + 0.2047885714285714 * t18 * (v2 + 0.8535911602209945 * PdotR) * t21 - 0.2708317460317460 * PdotR - 0.3604317460317460 * v2) * t10 + (0.6007195767195767 + t411 + t413 + (-0.1228228571428571 * PdotR - 0.1450742857142857 * v2) * t21) * t9 - 0.8408400000000000 * (-0.3834109956558936 + (v2 + 0.9285714285714286 * PdotR) * t21) * t21 * t3 - 0.12012 * (-0.2333333333333333 + t114) * t78 * t2 + 0.8008 * t431 + t434;
            t438 = t64 * t4;
            t512 = (0.6666666666666667 * t146 + t148 + t149) * PdotV;
            t514 = (0.2500000000000000 * t335 * RdotV + t512) * t334;
            t516 = t83 * p2 * PdotV;
            t522 = PdotV * t158;
            t525 = (-0.9 * t154 * PdotV - 0.4500000000000000 * t342 * RdotV - 0.3 * t522) * t228;
            t542 = (0.33 * t149 + 0.66 * t146 + 0.33 * t148) * t40;
            t544 = 0.99 * t236 * t21;
            t546 = 0.33 * t369;
            t547 = 0.8250000000000000 * t364;
            t548 = 0.33 * t360;
            t549 = 0.4950000000000000 * t371;
            t550 = 0.8250000000000000 * t366;
            t573 = 0.3000000000000000 * t387;
            t574 = 0.3 * t389;
            t594 = 0.6666666666666667 * t410;
            t614 = r * t283;
            t617 = 0.5498690476190476 * t433 * RdotV;
            t618 = t514 - 0.7878684807256236 * t516 + t525 + ((0.6357596371882086 * t149 + 0.6905215419501134 * t148 + 0.1363083900226757 * t146) * RdotV + 0.7088548752834467 * PdotV * (v2 + 0.8526894963292334 * PdotR) * t18) * t46 + (t542 + t544 + (-0.2546870748299320 * p2 + t546 + t547 + t548 + t549 + t550) * RdotV - 0.5406530612244898 * PdotV * (v2 + 0.8466119751035952 * PdotR)) * t45 + (-0.2384081632653061 * t247 - 0.1469585034013605 * (v2 + 0.9284586791588166 * PdotR) * PdotV * t21 - 0.2539229024943311 * t62 * (v2 + 0.7996070726915521 * PdotR) * RdotV + 0.1313892668178382 * PdotV) * t25 - 0.143 * (t573 + t574 + t391 - 0.3668585382871097 * t60 - 0.3416989359846503 * t396 - 0.1984371184371184 * t363 - 0.1453444967730682 * t358) * RdotV * t58 + 0.3132380952380952 * (0.9641304347826087 * t169 + (v2 + 0.8695652173913043 * PdotR) * t18 * t21 - 0.3361894144502840 * PdotR - 0.4303445346923608 * v2) * RdotV * t10 + 0.6435000000000000 * (0.3491331322019153 * 0.01 + t594 + t412 + (-0.2955908289241623 * PdotR - 0.3421516754850088 * v2) * t21) * RdotV * t9 - 0.8961333333333333 * (-0.5459867161994822 * 0.01 + (v2 + 0.9361702127659574 * PdotR) * t21) * t40 * t3 - 0.1041857142857143 * (-0.2867102396514161 + t114) * t78 * t2 + 0.6945714285714286 * t614 + t617;
            t619 = t192 * t8;
            t622 = t64 * t11;
            t679 = t336 + t345 + (0.4903111111111111 * t146 + 0.2488126984126984 * t148 + 0.2286412698412698 * t149) * t46 + (t355 + t357 + t361 - 0.1096063492063492 * p2 + t365 + t367 + t370 + t372) * t45 + (-0.2084317460317460 * t154 - 0.8513269841269841 * PdotV * (v2 + 0.9273738850272963 * PdotR) * RdotV - 0.7324867724867725 * (v2 + 0.7967783877492054 * PdotR) * t62) * t25 + (-t388 - t390 - t392 + 0.3548952380952381 * t60 + 0.9492571428571429 * t363 + 0.6936507936507936 * t358 + 0.1632736507936508 * t396) * t58 + (0.2935847619047619 * t169 + 0.2275428571428571 * (v2 + 0.8682320441988950 * PdotR) * t18 * t21 - 0.4231026455026455 * PdotR - 0.5408169312169312 * v2) * t10 + (0.1140467960023516 + t411 + t413 + (-0.1588749206349206 * PdotR - 0.1839758730158730 * v2) * t21) * t9 - 0.9342666666666667 * (-0.5210934190526027 * 0.01 + (v2 + 0.9357142857142857 * PdotR) * t21) * t21 * t3 - 0.12012 * (-0.3174603174603175 + t114) * t78 * t2 + 0.8897777777777778 * t431 + t434;
            t723 = t192 * t17;
            t774 = 0.264 * t247 * PdotV + 0.198 * t342 * t21 + 0.264 * t522 * RdotV + 0.6600000000000000 * t366 * PdotR + 0.33 * v2 * t366 + 0.66 * t363 * t359 + 0.66 * t368 * t358 + (0.33 * t364 - 0.9501804988662132 * p2) * PdotR + 0.6600000000000000 * t364 * v2 - 0.1117600000000000 * t148 - 0.2181601814058957 * t146;
            t890 = p2 * t341 * t45 * t58 - 0.8017233560090703 * t126 * t335 + (-0.9 * t21 * t335 - 0.72 * t512 * RdotV - 0.6 * t62 * (t148 + 0.6 * t146 + t149)) * t334 + (0.5544671201814059 * t338 + 0.1232154195011338 * p2 * t358 + (0.2671655328798186 * t148 + 0.5278820861678005 * t146) * PdotR + 0.5689433106575964 * v2 * t146 + 0.1439501133786848 * p2 * t363) * t83 + t774 * t228 + ((-0.2778068027210884 * t149 - 0.2986748299319728 * t148 - 0.5903455782312925 * t146) * t21 - 0.6115374149659864 * (v2 + 0.8690664768176560 * PdotR) * PdotV * t35 - 0.8584308390022676 * t360 - 0.1984308390022676 * t366 - 0.9878458049886621 * t369 - 0.1384707482993197 * t371 + 0.2775449735449735 * p2 - 0.2631383219954649 * t364) * t46 + ((-0.858 * t149 - 0.1716 * t146 - 0.858 * t148) * t78 - 0.3432 * t236 * t40 + (0.1086704943310658 * p2 - 0.429 * t366 - 0.429 * t364 - 0.2574 * t371 - 0.1716 * t360 - 0.1716 * t369) * t21 + 0.4581869569160998 * PdotV * (v2 + 0.8640769518634345 * PdotR) * RdotV + 0.4092221919879063 * (0.6443325109730134 * t358 + 0.1618096576703521 * t396 + t363) * t18) * t45 + (0.6142190476190476 * t387 + 0.5033600000000000 * (v2 + 0.9361471861471861 * PdotR) * PdotV * t40 + 0.1300165079365079 * (v2 + 0.8198324022346369 * PdotR) * t62 * t21 - 0.1098551534391534 * t60 - 0.3077306727135299 * t363 - 0.4930878306878307 * t396 - 0.1940835071806500 * t358) * t25 + (0.858 * t189 * p2 + 0.10296 * t99 * t121 + 0.4290 * t158 * t78 - 0.1781371428571429 * t169 + (-0.1088434285714286 * t358 - 0.1433813333333333 * t363 - 0.2510262857142857 * t396) * t21 + 0.7961945074326027 * PdotR + 0.1103617031997984 * v2) * t58 + (-0.1471489375997313 - 0.7206110476190476 * t410 - 0.9309980952380952 * (v2 + 0.8823873610298420 * PdotR) * t18 * t78 + (0.6750205291005291 * v2 + 0.5431881481481481 * PdotR) * t21) * t10 - 0.14586 * (0.7713109935332158 * 0.01 + 0.5714285714285714 * t410 + t412 + (-0.4459383753501401 * v2 - 0.3914098972922502 * PdotR) * t21) * t21 * t9 + 0.2022592000000000 * (-0.7146448322918911 * 0.01 + (0.9423076923076923 * PdotR + v2) * t21) * t78 * t3 + 0.1979528571428571 * (-0.3405847953216374 + t114) * t189 * t2 - 0.1319685714285714 * r * t433 - 0.9237800000000000 * t433 * t21;
            t970 = t514 - 0.8883900226757370 * t516 + t525 + ((0.7160770975056689 * t149 + 0.7708390022675737 * t148 + 0.1521473922902494 * t146) * RdotV + 0.7850566893424036 * t99 * (v2 + 0.8669882441292857 * PdotR)) * t46 + (t542 + t544 + (-0.3263841269841270 * p2 + t546 + t547 + t548 + t549 + t550) * RdotV - 0.6787791383219955 * (v2 + 0.8627735002418641 * PdotR) * PdotV) * t45 + (-0.2660204081632653 * t247 - 0.1627251700680272 * (v2 + 0.9353904166283454 * PdotR) * PdotV * t21 - 0.2793151927437642 * t62 * (v2 + 0.8178246115377746 * PdotR) * RdotV + 0.1936278155706727 * PdotV) * t25 - 0.143 * (t573 + t574 + t391 - 0.4609314495028781 * t60 - 0.4290465724751439 * t396 - 0.2452412349555207 * t363 - 0.1858880167451596 * t358) * RdotV * t58 + 0.3445619047619048 * RdotV * (0.9703557312252964 * t169 + (v2 + 0.8814229249011858 * PdotR) * t18 * t21 - 0.4516025037764168 * PdotR - 0.5605544214239866 * v2) * t10 + 0.6435000000000000 * (0.6102198017541933 * 0.01 + t594 + t412 + (-0.3714285714285714 * PdotR - 0.4232804232804233 * v2) * t21) * RdotV * t9 - 0.9857466666666667 * (-0.7143368845496505 * 0.01 + (v2 + 0.9419729206963250 * PdotR) * t21) * t40 * t3 - 0.1041857142857143 * (-0.3550326797385621 + t114) * t121 * t2 + 0.7640285714285714 * t614 + t617;
            t1028 = t336 + t345 + (0.5059174603174603 * t146 + 0.2567873015873016 * t148 + 0.2366158730158730 * t149) * t46 + (t355 + t357 + t361 - 0.1171707936507937 * p2 + t365 + t367 + t370 + t372) * t45 + (-0.2155873015873016 * t154 - 0.8781968253968254 * (v2 + 0.9295959975707800 * PdotR) * PdotV * RdotV - 0.7539894179894180 * t62 * (v2 + 0.8025739628361309 * PdotR)) * t25 + (-t388 - t390 - t392 + 0.3791032380952381 * t60 + 0.1009046349206349 * t363 + 0.7451898412698413 * t358 + 0.1744064761904762 * t396) * t58 + (0.3033904761904762 * t169 + 0.2346247619047619 * (v2 + 0.8722093230934095 * PdotR) * t18 * t21 - 0.4704753439153439 * PdotR - 0.5960296296296296 * v2) * t10 + (0.1309727454438566 + t411 + t413 + (-0.1704547428571429 * PdotR - 0.1963891809523810 * v2) * t21) * t9 - 0.9647733333333333 * t21 * (-0.5605544214239866 + (v2 + 0.9377470355731225 * PdotR) * t21) * t3 - 0.12012 * (-0.3174603174603175 + t114) * t78 * t2 + 0.9200302222222222 * t431 + t434;
            t723 = t192 * t17;
            t774 = 0.264 * t247 * PdotV + 0.198 * t342 * t21 + 0.264 * t522 * RdotV + 0.66 * t366 * PdotR + 0.33 * v2 * t366 + 0.66 * t363 * t359 + 0.66 * t368 * t358 + (0.33 * t364 - 0.9501804988662132 * p2) * PdotR + 0.66 * t364 * v2 - 0.11176 * t148 - 0.2181601814058957 * t146; 
            t890 = p2 * t341 * t45 * t58 - 0.8017233560090703 * t126 * t335 + (-0.9 * t21 * t335 - 0.72 * t512 * RdotV - 0.6 * t62 * (t148 + 0.6 * t146 + t149)) * t334 + (0.5544671201814059 * t338 + 0.1232154195011338 * p2 * t358 + (0.2671655328798186 * t148 + 0.5278820861678005 * t146) * PdotR + 0.5689433106575964 * v2 * t146 + 0.1439501133786848 * p2 * t363) * t83 + t774 * t228 + ( (-0.2778068027210884 * t149 - 0.2986748299319728 * t148 - 0.5903455782312925 * t146) * t21 - 0.6115374149659864 * (v2 + 0.8690664768176560 * PdotR) * PdotV * t35 - 0.8584308390022676 * t360 - 0.1984308390022676 * t366 - 0.9878458049886621 * t369 - 0.1384707482993197 * t371 + 0.2775449735449735 * p2 - 0.2631383219954649 * t364) * t46 + ((-0.858 * t149 - 0.1716 * t146 - 0.858 * t148) * t78 - 0.3432 * t236 * t40 + (0.1086704943310658 * p2 - 0.429 * t366 - 0.429 * t364 - 0.2574 * t371 - 0.1716 * t360 - 0.1716 * t369) * t21 + 0.4581869569160998 * PdotV * (v2 + 0.8640769518634345 * PdotR) * RdotV + 0.4092221919879063 * (0.6443325109730134 * t358 + 0.1618096576703521 * t396 + t363) * t18) * t45 + (0.6142190476190476 * t387 + 0.5033600000000000 * (v2 + 0.9361471861471861 * PdotR) * PdotV * t40 + 0.1300165079365079 * (v2 + 0.8198324022346369 * PdotR) * t62 * t21 - 0.1098551534391534 * t60 - 0.3077306727135299 * t363 - 0.4930878306878307 * t396 - 0.1940835071806500 * t358) * t25 + (0.858 * t189 * p2 + 0.10296 * t99 * t121 + 0.4290 * t158 * t78 - 0.1781371428571429 * t169 + (-0.1088434285714286 * t358 - 0.1433813333333333 * t363 - 0.2510262857142857 * t396) * t21 + 0.7961945074326027 * PdotR + 0.1103617031997984 * v2) * t58 + (-0.1471489375997313 - 0.7206110476190476 * t410 - 0.9309980952380952 * (v2 + 0.8823873610298420 * PdotR) * t18 * t78 + (0.6750205291005291 * v2 + 0.5431881481481481 * PdotR) * t21) * t10 - 0.14586 * (0.7713109935332158 * 0.01 + 0.5714285714285714 * t410 + t412 + (-0.4459383753501401 * v2 - 0.3914098972922502 * PdotR) * t21) * t21 * t9 + 0.2022592000000000 * (-0.7146448322918911 * 0.1 + (0.9423076923076923 * PdotR + v2) * t21) * t78 * t3 + 0.1979528571428571 * (-0.3405847953216374 + t114) * t189 * t2 - 0.1319685714285714 * r * t433 - 0.9237800000000000 * t433 * t21;
            t970 = t514 - 0.8883900226757370 * t516 + t525 + ((0.7160770975056689 * t149 + 0.7708390022675737 * t148 + 0.1521473922902494 * t146) * RdotV + 0.7850566893424036 * t99 * (v2 + 0.8669882441292857 * PdotR)) * t46 + (t542 + t544 + (-0.3263841269841270 * p2 + t546 + t547 + t548 + t549 + t550) * RdotV - 0.6787791383219955 * (v2 + 0.8627735002418641 * PdotR) * PdotV) * t45 + (-0.2660204081632653 * t247 - 0.1627251700680272 * (v2 + 0.9353904166283454 * PdotR) * PdotV * t21 - 0.2793151927437642 * t62 * (v2 + 0.8178246115377746 * PdotR) * RdotV + 0.1936278155706727 * PdotV) * t25 - 0.143 * (t573 + t574 + t391 - 0.4609314495028781 * t60 - 0.4290465724751439 * t396 - 0.2452412349555207 * t363 - 0.1858880167451596 * t358) * RdotV * t58 + 0.3445619047619048 * RdotV * (0.9703557312252964 * t169 + (v2 + 0.8814229249011858 * PdotR) * t18 * t21 - 0.4516025037764168 * PdotR - 0.5605544214239866 * v2) * t10 + 0.6435000000000000 * (0.6102198017541933 * 0.001 + t594 + t412 + (-0.3714285714285714 * PdotR - 0.4232804232804233 * v2) * t21) * RdotV * t9 - 0.9857466666666667 * (-0.7143368845496505 * 0.1 + (v2 + 0.9419729206963250 * PdotR) * t21) * t40 * t3 - 0.1041857142857143 * (-0.3550326797385621 + t114) * t121 * t2 + 0.7640285714285714 * t614 + t617;
            t1028 = t336 + t345 + (0.5059174603174603 * t146 + 0.2567873015873016 * t148 + 0.2366158730158730 * t149) * t46 + (t355 + t357 + t361 - 0.1171707936507937 * p2 + t365 + t367 + t370 + t372) * t45 + (-0.2155873015873016 * t154 - 0.8781968253968254 * (v2 + 0.9295959975707800 * PdotR) * PdotV * RdotV - 0.7539894179894180 * t62 * (v2 + 0.8025739628361309 * PdotR)) * t25 + (-t388 - t390 - t392 + 0.3791032380952381 * t60 + 0.1009046349206349 * t363 + 0.7451898412698413 * t358 + 0.1744064761904762 * t396) * t58 + (0.3033904761904762 * t169 + 0.2346247619047619 * (v2 + 0.8722093230934095 * PdotR) * t18 * t21 - 0.4704753439153439 * PdotR - 0.5960296296296296 * v2) * t10 + (0.1309727454438566 + t411 + t413 + (-0.1704547428571429 * 0.1 * PdotR - 0.1963891809523810 * 0.1 * v2) * t21) * t9 - 0.9647733333333333 * t21 * (-0.5605544214239866 * 0.1 + (v2 + 0.9377470355731225 * PdotR) * t21) * t3 - 0.12012 * (-0.3174603174603175 + t114) * t78 * t2 + 0.9200302222222222 * t431 + t434;

            FGHT(1) = 0.0;
            FGHT(2) = k;
            FGHT(3) = 0.0;
            FGHT(4) = k;
            FGHT(5) = -0.5000000000000000 * t1 * t4;
            FGHT(6) = 0.0;
            FGHT(7) = 0.5000000000000000 * t1;
            FGHT(8) = 0.0;
            FGHT(9) = 0.5000000000000000 * t8 * t11 * RdotV;
            FGHT(10) = -0.1666666666666667 * t8 * t4;
            FGHT(11) = 0.0;
            FGHT(12) = 0.0;
            FGHT(13) = 0.1250000000000000 * t17 * (t19 - 0.6666666666666667 * r - t22) * t26;
            FGHT(14) = 0.2500000000000000 * t17 * t11 * RdotV;
            FGHT(15) = -0.4166666666666667 * t17 * t4;
            FGHT(16) = 0.0;
            FGHT(17) = 0.7500000000000000 * (t34 - t37 + 0.3333333333333333 * t38 + t41) * t43 * t47;
            FGHT(18) = 0.7500000000000000 * t43 * (t19 - 0.8888888888888889 * r - t22) * t26;
            FGHT(19) = 0.7500000000000000 * t43 * t11 * RdotV;
            FGHT(20) = 0.0;
            FGHT(21) = 0.1250000000000000 * (t59 + t65 + (0.6333333333333333 * PdotR + 0.7333333333333333 * v2) * t3 + (-0.2444444444444444 + t73) * t2 - 0.4666666666666667 * t76 - t79) * t81 * t84;
            FGHT(22) = 0.5000000000000000 * t81 * (t34 - t37 + 0.4166666666666667 * t38 + t41) * t47;
            FGHT(23) = 0.2500000000000000 * (t19 - 0.9444444444444444 * r - t22) * t81 * t26;
            FGHT(24) = 0.0;

            FGHT(25) = -0.8928571428571429E-1 * (t101 - 0.74 * t102 - t108 + 0.504E1 * (0.8888888888888889 * PdotR + v2) * RdotV * t3 + 0.21E2 * (-0.08 + t114) * RdotV * t2 - 0.14E2 * t119 - t122) * t124 * t127;
            FGHT(26) = 0.8928571428571429E-2 * (t59 + t65 + (0.88 * v2 + 0.78 * PdotR) * t3 + (-0.3822222222222222 + t73) * t2 - 0.56E2 * t76 - t79) * t124 * t84;
            FGHT(27) = 0.1785714285714286E-1 * (t34 - t37 + 0.44 * t38 + t41) * t124 * t47;
            FGHT(28) = 0.0;
            FGHT(29) = -0.1674107142857143E-1 * (t151 - 0.76 * t152 + t161 + (0.2032E2 * t60 + 0.536E1 * (v2 + 0.8059701492537313 * PdotR) * t18) * t10 + (t170 + t172 - 0.3893333333333333 * v2 - 0.3093333333333333 * PdotR) * t9 + (0.8651851851851852 + (-0.896E2 * v2 - 0.812E2 * PdotR) * t21) * t3 - 0.231E3 * (-0.1292929292929293 + t114) * t21 * t2 + 0.154E3 * t187 + t190) * t192 * t195;
            FGHT(30) = -0.6696428571428571E-1 * t192 * (t101 - 0.88 * t102 - t108 + 0.588E1 * (0.9047619047619048 * PdotR + v2) * RdotV * t3 + (t207 - 0.2426666666666667 * RdotV) * t2 - 0.1633333333333333E2 * t119 - t122) * t127;
            FGHT(31) = 0.3348214285714286E-2 * t192 * (t59 + t65 + (0.92 * v2 + 0.82 * PdotR) * t3 + (-0.4207407407407407 + t73) * t2 - 0.588E2 * t76 - t79) * t84;
            FGHT(32) = 0.0;
            FGHT(33) = -0.1302083333333333E-1 * t227 * (t230 + t239 + (0.52 * t97 + 0.1077714285714286E2 * PdotV * (0.9056203605514316 * PdotR + v2)) * t25 + (t248 + t250 + t252 - 0.3954285714285714 * PdotV) * t58 - 0.4725714285714286E2 * (0.1909310761789601 * t60 + (v2 + 0.8331318016928658 * PdotR) * t18) * RdotV * t10 - 0.231E3 * RdotV * (t169 + t171 - 0.1464440321583179 * v2 - 0.1207173778602350 * PdotR) * t9 + 0.3256E3 * RdotV * (-0.2308802308802309E-1 + (v2 + 0.9189189189189189 * PdotR) * t21) * t3 + 0.6006E3 * (-0.1807081807081807 + t114) * t40 * t2 - 0.4004E3 * t281 - t284) * t289;
            FGHT(34) = -0.1302083333333333E-1 * (t151 - 0.8971428571428571 * t152 + t161 + (0.2354285714285714E2 * t60 + 0.6125714285714286E1 * (v2 + 0.8302238805970149 * PdotR) * t18) * t10 + (t170 + t172 - 0.5302857142857143 * v2 - 0.4365714285714286 * PdotR) * t9 + (0.1510264550264550E1 + (-0.1024E3 * v2 - 0.94E2 * PdotR) * t21) * t3 - 0.231E3 * (-0.1766233766233766 + t114) * t21 * t2 + 0.176E3 * t187 + t190) * t227 * t195;
            FGHT(35) = -0.2604166666666667E-1 * (t101 - 0.9171428571428571 * t102 - t108 + 0.6125714285714286E1 * (0.9085820895522388 * PdotR + v2) * RdotV * t3 + (t207 - 0.2651428571428571 * RdotV) * t2 - 0.1706666666666667E2 * t119 - t122) * t227 * t127;
            FGHT(36) = 0.0;
            FGHT(37) = -0.1302083333333333E-2 * t333 * t435 * t438;
            FGHT(38) = -0.1041666666666667E-1 * t333 * (t230 + t239 + (0.599 * t97 + 0.1225428571428571E2 * PdotV * (0.9169969689904407 * PdotR + v2)) * t25 + (t248 + t250 + t252 - 0.5308571428571429 * PdotV) * t58 - 0.5316428571428571E2 * (0.1929060862555421 * t60 + t18 * (v2 + 0.8516727126158807 * PdotR)) * RdotV * t10 - 0.231E3 * (t169 + t171 - 0.1915893630179344 * v2 - 0.1620902906617192 * PdotR) * RdotV * t9 + ((0.3663E3 * v2 + 0.3399E3 * PdotR) * t40 - 0.1207936507936508E2 * RdotV) * t3 + 0.6006E3 * (-0.2368742368742369 + t114) * t40 * t2 - 0.45045E3 * t281 - t284) * t289;
            FGHT(39) = -0.5208333333333333;
            FGHT(40) = 0.0;
            FGHT(41) = 0.2982954545454545E-1 * t618 * t619 * t622;
            FGHT(42) = -0.1065340909090909E-2 * t679 * t619 * t438;
            FGHT(43) = -0.4261363636363636E-2 * t619 * (t230 + t239 + (0.6211746031746032 * t97 + 0.1266634920634921E2 * (0.9196972355196872 * PdotR + v2) * PdotV) * t25 + (t248 + t250 + t252 - 0.5694444444444444 * PdotV) * t58 - 0.5493650793650794E2 * (0.1934209765963594 * t60 + (v2 + 0.8564576711932967 * PdotR) * t18) * RdotV * t10 - 0.231E3 * (t169 + t171 - 0.2054669140383426 * v2 - 0.1749171992029135 * PdotR) * RdotV * t9 + 0.3792380952380952E3 * (-0.3565154305485797E-1 + (v2 + 0.9303867403314917 * PdotR) * t21) * RdotV * t3 + 0.6006E3 * t40 * (-0.2552667790763029 * t114) * t2 - 0.4671333333333333E3 * t281 - t284) * t289;
            FGHT(44) = 0.0;
            FGHT(45) = 0.3107244318181818E-2 * t723 * t890 * t64 * t26;
            FGHT(46) = 0.2485795454545455E-1 * t970 * t723 * t622;
            FGHT(47) = -0.4438920454545455E-3 * t723 * t1028 * t438;
            FGHT(48) = 0.0;