BEGIN{h=0;}
{a = 0;}
($2 == "veg_parm_ID")&&($1==51){h=1;}
($2 == "rootzone.depth") && (h==1) {printf("%f %s\n",$1=1.07925291,$2); a=1;}
($2 == "cs.cpool") && (h==1) {printf("%f %s\n",$1=-0.17221906,$2); a=1;}
($2 == "cs.leafc") && (h==1) {printf("%f %s\n",$1=0.15188539,$2); a=1;}
($2 == "cs.dead_leafc") && (h==1) {printf("%f %s\n",$1=0.00000000,$2); a=1;}
($2 == "cs.leafc_store") && (h==1) {printf("%f %s\n",$1=0.05606024,$2); a=1;}
($2 == "cs.leafc_transfer") && (h==1) {printf("%f %s\n",$1=0.00000000,$2); a=1;}
($2 == "cs.live_stemc") && (h==1) {printf("%f %s\n",$1=0.12043962,$2); a=1;}
($2 == "cs.livestemc_store") && (h==1) {printf("%f %s\n",$1=0.00773109,$2); a=1;}
($2 == "cs.livestemc_transfer") && (h==1) {printf("%f %s\n",$1=0.00000000,$2); a=1;}
($2 == "cs.dead_stemc") && (h==1) {printf("%f %s\n",$1=0.56667709,$2); a=1;}
($2 == "cs.deadstemc_store") && (h==1) {printf("%f %s\n",$1=0.00165070,$2); a=1;}
($2 == "cs.deadstemc_transfer") && (h==1) {printf("%f %s\n",$1=0.00000000,$2); a=1;}
($2 == "cs.live_crootc") && (h==1) {printf("%f %s\n",$1=0.08622412,$2); a=1;}
($2 == "cs.livecrootc_store") && (h==1) {printf("%f %s\n",$1=0.00525789,$2); a=1;}
($2 == "cs.livecrootc_transfer") && (h==1) {printf("%f %s\n",$1=0.00000000,$2); a=1;}
($2 == "cs.dead_crootc") && (h==1) {printf("%f %s\n",$1=0.47847014,$2); a=1;}
($2 == "cs.deadcrootc_store") && (h==1) {printf("%f %s\n",$1=0.00112263,$2); a=1;}
($2 == "cs.deadcrootc_transfer") && (h==1) {printf("%f %s\n",$1=0.00000000,$2); a=1;}
($2 == "cs.frootc") && (h==1) {printf("%f %s\n",$1=0.13139143,$2); a=1;}
($2 == "cs.frootc_store") && (h==1) {printf("%f %s\n",$1=0.02394502,$2); a=1;}
($2 == "cs.frootc_transfer") && (h==1) {printf("%f %s\n",$1=0.00000000,$2); a=1;}
($2 == "cs.cwdc") && (h==1) {printf("%f %s\n",$1=0.07775238,$2); a=1;}
($2 == "epv.prev_leafcalloc") && (h==1) {printf("%f %s\n",$1=0.20424644,$2); a=1;}
($2 == "ns.npool") && (h==1) {printf("%f %s\n",$1=0.00688272,$2); a=1;}
($2 == "ns.leafn") && (h==1) {printf("%f %s\n",$1=0.00379713,$2); a=1;}
($2 == "ns.dead_leafn") && (h==1) {printf("%f %s\n",$1=0.00000000,$2); a=1;}
($2 == "ns.leafn_store") && (h==1) {printf("%f %s\n",$1=0.00140151,$2); a=1;}
($2 == "ns.leafn_transfer") && (h==1) {printf("%f %s\n",$1=0.00000000,$2); a=1;}
($2 == "ns.live_stemn") && (h==1) {printf("%f %s\n",$1=0.00200711,$2); a=1;}
($2 == "ns.livestemn_store") && (h==1) {printf("%f %s\n",$1=0.00012885,$2); a=1;}
($2 == "ns.livestemn_transfer") && (h==1) {printf("%f %s\n",$1=0.00000000,$2); a=1;}
($2 == "ns.dead_stemn") && (h==1) {printf("%f %s\n",$1=0.00151116,$2); a=1;}
($2 == "ns.deadstemn_store") && (h==1) {printf("%f %s\n",$1=0.00000440,$2); a=1;}
($2 == "ns.deadstemn_transfer") && (h==1) {printf("%f %s\n",$1=0.00000000,$2); a=1;}
($2 == "ns.live_crootn") && (h==1) {printf("%f %s\n",$1=0.00143732,$2); a=1;}
($2 == "ns.livecrootn_store") && (h==1) {printf("%f %s\n",$1=0.00008763,$2); a=1;}
($2 == "ns.livecrootn_transfer") && (h==1) {printf("%f %s\n",$1=0.00000000,$2); a=1;}
($2 == "ns.dead_crootn") && (h==1) {printf("%f %s\n",$1=0.00127569,$2); a=1;}
($2 == "ns.deadcrootn_store") && (h==1) {printf("%f %s\n",$1=0.00000299,$2); a=1;}
($2 == "ns.deadcrootn_transfer") && (h==1) {printf("%f %s\n",$1=0.00000000,$2); a=1;}
($2 == "ns.frootn") && (h==1) {printf("%f %s\n",$1=0.00218986,$2); a=1;}
($2 == "ns.frootn_store") && (h==1) {printf("%f %s\n",$1=0.00039908,$2); a=1;}
($2 == "ns.frootn_transfer") && (h==1) {printf("%f %s\n",$1=0.00000000,$2); a=1;}
($2 == "ns.cwdn") && (h==1) {printf("%f %s\n",$1=-0.00077090,$2); a=1;}
($2 == "ns.retransn") && (h==1) {printf("%f %s\n",$1=0.00080553,$2); a=1;h=0;}

(a == 0) {printf("%s	%s\n",$1,$2);}
