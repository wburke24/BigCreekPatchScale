{a = 0;}
($2 == "litter_cs.litr1c") {printf("%f %s\n",$1=0.00000720,$2); a=1;}
($2 == "litter_ns.litr1n") {printf("%f %s\n",$1=0.00000296,$2); a=1;}
($2 == "litter_cs.litr2c") {printf("%f %s\n",$1=0.05865499,$2); a=1;}
($2 == "litter_cs.litr3c") {printf("%f %s\n",$1=0.10694621,$2); a=1;}
($2 == "litter_cs.litr4c") {printf("%f %s\n",$1=0.20805751,$2); a=1;}
($2 == "soil_cs.soil1c") {printf("%f %s\n",$1=0.00894134,$2); a=1;}
($2 == "soil_ns.sminn") {printf("%f %s\n",$1=0.00000000,$2); a=1;}
($2 == "soil_ns.nitrate") {printf("%f %s\n",$1=0.00000000,$2); a=1;}
($2 == "soil_cs.soil2c") {printf("%f %s\n",$1=0.15463644,$2); a=1;}
($2 == "soil_cs.soil3c") {printf("%f %s\n",$1=1.56780744,$2); a=1;}
($2 == "soil_cs.soil4c") {printf("%f %s\n",$1=9.94377106,$2); a=1;}

(a == 0) {printf("%s	%s\n",$1,$2);}
