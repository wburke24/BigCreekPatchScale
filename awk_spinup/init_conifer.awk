BEGIN{h=0;}
{a = 0;}
($2 == "veg_parm_ID")&&($1==7){h=1;}
($2 == "rootzone.depth") && (h==1) {printf("%f %s\n",$1=4.46856801,$2); a=1;}
($2 == "cs.cpool") && (h==1) {printf("%f %s\n",$1=-0.89499289,$2); a=1;}
($2 == "cs.leafc") && (h==1) {printf("%f %s\n",$1=1.21796096,$2); a=1;}
($2 == "cs.dead_leafc") && (h==1) {printf("%f %s\n",$1=0.00000000,$2); a=1;}
($2 == "cs.leafc_store") && (h==1) {printf("%f %s\n",$1=0.18051718,$2); a=1;}
($2 == "cs.leafc_transfer") && (h==1) {printf("%f %s\n",$1=0.00000000,$2); a=1;}
($2 == "cs.live_stemc") && (h==1) {printf("%f %s\n",$1=0.08049065,$2); a=1;}
($2 == "cs.livestemc_store") && (h==1) {printf("%f %s\n",$1=0.04044029,$2); a=1;}
($2 == "cs.livestemc_transfer") && (h==1) {printf("%f %s\n",$1=0.00000000,$2); a=1;}
($2 == "cs.dead_stemc") && (h==1) {printf("%f %s\n",$1=11.95978001,$2); a=1;}
($2 == "cs.deadstemc_store") && (h==1) {printf("%f %s\n",$1=0.02395436,$2); a=1;}
($2 == "cs.deadstemc_transfer") && (h==1) {printf("%f %s\n",$1=0.00000000,$2); a=1;}
($2 == "cs.live_crootc") && (h==1) {printf("%f %s\n",$1=0.04020466,$2); a=1;}
($2 == "cs.livecrootc_store") && (h==1) {printf("%f %s\n",$1=0.01986724,$2); a=1;}
($2 == "cs.livecrootc_transfer") && (h==1) {printf("%f %s\n",$1=0.00000000,$2); a=1;}
($2 == "cs.dead_crootc") && (h==1) {printf("%f %s\n",$1=8.05220923,$2); a=1;}
($2 == "cs.deadcrootc_store") && (h==1) {printf("%f %s\n",$1=0.01176814,$2); a=1;}
($2 == "cs.deadcrootc_transfer") && (h==1) {printf("%f %s\n",$1=0.00000000,$2); a=1;}
($2 == "cs.frootc") && (h==1) {printf("%f %s\n",$1=0.85683872,$2); a=1;}
($2 == "cs.frootc_store") && (h==1) {printf("%f %s\n",$1=0.15188185,$2); a=1;}
($2 == "cs.frootc_transfer") && (h==1) {printf("%f %s\n",$1=0.00000000,$2); a=1;}
($2 == "cs.cwdc") && (h==1) {printf("%f %s\n",$1=0.32395774,$2); a=1;}
($2 == "epv.prev_leafcalloc") && (h==1) {printf("%f %s\n",$1=1.45000858,$2); a=1;}
($2 == "ns.npool") && (h==1) {printf("%f %s\n",$1=0.00274899,$2); a=1;}
($2 == "ns.leafn") && (h==1) {printf("%f %s\n",$1=0.03044902,$2); a=1;}
($2 == "ns.dead_leafn") && (h==1) {printf("%f %s\n",$1=0.00000000,$2); a=1;}
($2 == "ns.leafn_store") && (h==1) {printf("%f %s\n",$1=0.00451293,$2); a=1;}
($2 == "ns.leafn_transfer") && (h==1) {printf("%f %s\n",$1=0.00000000,$2); a=1;}
($2 == "ns.live_stemn") && (h==1) {printf("%f %s\n",$1=0.00134207,$2); a=1;}
($2 == "ns.livestemn_store") && (h==1) {printf("%f %s\n",$1=0.00067400,$2); a=1;}
($2 == "ns.livestemn_transfer") && (h==1) {printf("%f %s\n",$1=0.00000000,$2); a=1;}
($2 == "ns.dead_stemn") && (h==1) {printf("%f %s\n",$1=0.03189250,$2); a=1;}
($2 == "ns.deadstemn_store") && (h==1) {printf("%f %s\n",$1=0.00006388,$2); a=1;}
($2 == "ns.deadstemn_transfer") && (h==1) {printf("%f %s\n",$1=0.00000000,$2); a=1;}
($2 == "ns.live_crootn") && (h==1) {printf("%f %s\n",$1=0.00067015,$2); a=1;}
($2 == "ns.livecrootn_store") && (h==1) {printf("%f %s\n",$1=0.00033112,$2); a=1;}
($2 == "ns.livecrootn_transfer") && (h==1) {printf("%f %s\n",$1=0.00000000,$2); a=1;}
($2 == "ns.dead_crootn") && (h==1) {printf("%f %s\n",$1=0.02147201,$2); a=1;}
($2 == "ns.deadcrootn_store") && (h==1) {printf("%f %s\n",$1=0.00003138,$2); a=1;}
($2 == "ns.deadcrootn_transfer") && (h==1) {printf("%f %s\n",$1=0.00000000,$2); a=1;}
($2 == "ns.frootn") && (h==1) {printf("%f %s\n",$1=0.01428065,$2); a=1;}
($2 == "ns.frootn_store") && (h==1) {printf("%f %s\n",$1=0.00253136,$2); a=1;}
($2 == "ns.frootn_transfer") && (h==1) {printf("%f %s\n",$1=0.00000000,$2); a=1;}
($2 == "ns.cwdn") && (h==1) {printf("%f %s\n",$1=-0.00220383,$2); a=1;}
($2 == "ns.retransn") && (h==1) {printf("%f %s\n",$1=0.00359515,$2); a=1;h=0;}

(a == 0) {printf("%s	%s\n",$1,$2);}
