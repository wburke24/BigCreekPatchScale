BEGIN{h=0;}
{a = 0;}
($2 == "veg_parm_ID")&&($1==50){h=1;}
($2 == "rootzone.depth") && (h==1) {printf("%f %s\n",$1=3.09238901,$2); a=1;}
($2 == "cs.cpool") && (h==1) {printf("%f %s\n",$1=-0.31992622,$2); a=1;}
($2 == "cs.leafc") && (h==1) {printf("%f %s\n",$1=0.59809751,$2); a=1;}
($2 == "cs.dead_leafc") && (h==1) {printf("%f %s\n",$1=0.00000000,$2); a=1;}
($2 == "cs.leafc_store") && (h==1) {printf("%f %s\n",$1=0.00000000,$2); a=1;}
($2 == "cs.leafc_transfer") && (h==1) {printf("%f %s\n",$1=0.29944004,$2); a=1;}
($2 == "cs.live_stemc") && (h==1) {printf("%f %s\n",$1=0.41808279,$2); a=1;}
($2 == "cs.livestemc_store") && (h==1) {printf("%f %s\n",$1=0.00000000,$2); a=1;}
($2 == "cs.livestemc_transfer") && (h==1) {printf("%f %s\n",$1=0.04344035,$2); a=1;}
($2 == "cs.dead_stemc") && (h==1) {printf("%f %s\n",$1=1.11089716,$2); a=1;}
($2 == "cs.deadstemc_store") && (h==1) {printf("%f %s\n",$1=0.00000000,$2); a=1;}
($2 == "cs.deadstemc_transfer") && (h==1) {printf("%f %s\n",$1=0.00492769,$2); a=1;}
($2 == "cs.live_crootc") && (h==1) {printf("%f %s\n",$1=0.34381953,$2); a=1;}
($2 == "cs.livecrootc_store") && (h==1) {printf("%f %s\n",$1=0.00000000,$2); a=1;}
($2 == "cs.livecrootc_transfer") && (h==1) {printf("%f %s\n",$1=0.03400810,$2); a=1;}
($2 == "cs.dead_crootc") && (h==1) {printf("%f %s\n",$1=1.03143472,$2); a=1;}
($2 == "cs.deadcrootc_store") && (h==1) {printf("%f %s\n",$1=0.00000000,$2); a=1;}
($2 == "cs.deadcrootc_transfer") && (h==1) {printf("%f %s\n",$1=0.00385773,$2); a=1;}
($2 == "cs.frootc") && (h==1) {printf("%f %s\n",$1=0.43722315,$2); a=1;}
($2 == "cs.frootc_store") && (h==1) {printf("%f %s\n",$1=0.00000000,$2); a=1;}
($2 == "cs.frootc_transfer") && (h==1) {printf("%f %s\n",$1=0.12231667,$2); a=1;}
($2 == "cs.cwdc") && (h==1) {printf("%f %s\n",$1=0.10941957,$2); a=1;}
($2 == "epv.prev_leafcalloc") && (h==1) {printf("%f %s\n",$1=0.89828786,$2); a=1;}
($2 == "ns.npool") && (h==1) {printf("%f %s\n",$1=0.00135407,$2); a=1;}
($2 == "ns.leafn") && (h==1) {printf("%f %s\n",$1=0.01495244,$2); a=1;}
($2 == "ns.dead_leafn") && (h==1) {printf("%f %s\n",$1=0.00000000,$2); a=1;}
($2 == "ns.leafn_store") && (h==1) {printf("%f %s\n",$1=0.00000000,$2); a=1;}
($2 == "ns.leafn_transfer") && (h==1) {printf("%f %s\n",$1=0.00748600,$2); a=1;}
($2 == "ns.live_stemn") && (h==1) {printf("%f %s\n",$1=0.00696805,$2); a=1;}
($2 == "ns.livestemn_store") && (h==1) {printf("%f %s\n",$1=0.00000000,$2); a=1;}
($2 == "ns.livestemn_transfer") && (h==1) {printf("%f %s\n",$1=0.00072401,$2); a=1;}
($2 == "ns.dead_stemn") && (h==1) {printf("%f %s\n",$1=0.00296239,$2); a=1;}
($2 == "ns.deadstemn_store") && (h==1) {printf("%f %s\n",$1=0.00000000,$2); a=1;}
($2 == "ns.deadstemn_transfer") && (h==1) {printf("%f %s\n",$1=0.00001314,$2); a=1;}
($2 == "ns.live_crootn") && (h==1) {printf("%f %s\n",$1=0.00573033,$2); a=1;}
($2 == "ns.livecrootn_store") && (h==1) {printf("%f %s\n",$1=0.00000000,$2); a=1;}
($2 == "ns.livecrootn_transfer") && (h==1) {printf("%f %s\n",$1=0.00056680,$2); a=1;}
($2 == "ns.dead_crootn") && (h==1) {printf("%f %s\n",$1=0.00275049,$2); a=1;}
($2 == "ns.deadcrootn_store") && (h==1) {printf("%f %s\n",$1=0.00000000,$2); a=1;}
($2 == "ns.deadcrootn_transfer") && (h==1) {printf("%f %s\n",$1=0.00001029,$2); a=1;}
($2 == "ns.frootn") && (h==1) {printf("%f %s\n",$1=0.00728705,$2); a=1;}
($2 == "ns.frootn_store") && (h==1) {printf("%f %s\n",$1=0.00000000,$2); a=1;}
($2 == "ns.frootn_transfer") && (h==1) {printf("%f %s\n",$1=0.00203861,$2); a=1;}
($2 == "ns.cwdn") && (h==1) {printf("%f %s\n",$1=-0.00294779,$2); a=1;}
($2 == "ns.retransn") && (h==1) {printf("%f %s\n",$1=0.00480572,$2); a=1;h=0;}

(a == 0) {printf("%s	%s\n",$1,$2);}
