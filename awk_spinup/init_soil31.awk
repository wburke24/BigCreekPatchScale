{a = 0;}
($2 == "litter_cs.litr1c") {printf("%f %s\n",$1=0.07280740,$2); a=1;}
($2 == "litter_ns.litr1n") {printf("%f %s\n",$1=0.00290842,$2); a=1;}
($2 == "litter_cs.litr2c") {printf("%f %s\n",$1=0.32813543,$2); a=1;}
($2 == "litter_cs.litr3c") {printf("%f %s\n",$1=0.17280141,$2); a=1;}
($2 == "litter_cs.litr4c") {printf("%f %s\n",$1=0.57037561,$2); a=1;}
($2 == "soil_cs.soil1c") {printf("%f %s\n",$1=0.03057119,$2); a=1;}
($2 == "soil_ns.sminn") {printf("%f %s\n",$1=0.00000000,$2); a=1;}
($2 == "soil_ns.nitrate") {printf("%f %s\n",$1=0.00000000,$2); a=1;}
($2 == "soil_cs.soil2c") {printf("%f %s\n",$1=0.19131987,$2); a=1;}
($2 == "soil_cs.soil3c") {printf("%f %s\n",$1=1.73252175,$2); a=1;}
($2 == "soil_cs.soil4c") {printf("%f %s\n",$1=10.64227175,$2); a=1;}

(a == 0) {printf("%s	%s\n",$1,$2);}
