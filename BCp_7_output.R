# -=================================================-
# -============= Output_BigCreek_Patch ==============
# ==================================================-

setwdhere()

library(data.table)
library(tidyverse)
library(plotly)
library(chron)
library(RHESSysIOinR)
library(randomForest)
library(randomForestExplainer)
library(gridExtra)
library(scales)
library(grid)

load("scenarios.rdata")
load("scen.rdata")

scen = plyr::rename(scen, c("start_period" = "Climatic Period",
                            "treatments" = "Treatment Method & Intensity",
                            "interval" = "Treatment Interval",
                            "soils" = "Soil Water Capacity",
                            "sharing" = "Root Sharing Coefficient",
                            "treat_type" = "Treatment Type",
                            "treat_intensity" = "Treatment Intensity",
                            "presc_fire" = "Prescribed Fire",
                            "veg" = "Vegetation Type",
                            "aspect" = "Aspect",
                            "clim_chg" = "Climate Change"))

# -=================================================-
# -=============== Helper Functions =================
# ==================================================-

read_robj = function(path) {
  DT_all = NULL
  load(path)
  if ("familyID" %in% colnames(DT_all)) {
    DT_all$familyID = NULL
  }
  days = length(DT_all$run[DT_all$run == 1])
  DT_all[, sim_day :=  rep(c(1:days), nrow(DT_all)/days), ]
  DT_all[, mn := lubridate::month(date), ]
  DT_all[, yr := lubridate::year(date), ]
  DT_all[, wy := data.table::fifelse(mn >= 10, yr + 1, yr), ]
  DT_all[, wym := data.table::fifelse(mn >= 10, mn - 9, mn + 3), ]
  DT_all[, sim_yr := wy - min(DT_all$yr), ]
  DT_all[, mn := NULL]
  DT_all[, yr := NULL]
  return(DT_all)
}

# chg from treatment
treat_change = function(data) {
  dt_merge = merge.data.table(data[data$"Treatment Method & Intensity" != "no treatment",],
                              data[data$"Treatment Method & Intensity" == "no treatment", 
                                   c("Climatic Period", 
                                     "Soil Depth",
                                     "Root Sharing Coefficient",
                                     "Vegetation Type", 
                                     "Aspect", 
                                     "Climate Change",
                                     "V1")],
                              by = c("Climatic Period", 
                                     "Soil Depth",
                                     "Root Sharing Coefficient",
                                     "Vegetation Type", 
                                     "Aspect", 
                                     "Climate Change"))
  
  dt_merge$treat_chg = dt_merge$V1.x - dt_merge$V1.y
  colnames(dt_merge)[colnames(dt_merge) == "V1.x"] = "value"
  colnames(dt_merge)[colnames(dt_merge) == "V1.y"] = "notreat_value"
  
  return(dt_merge)
}

treat_change2 = function(data) {
  
  trt = data[data$"Treatment Method & Intensity" != "no treatment", ]
  
  no_trt = data[data$"Treatment Method & Intensity" == "no treatment", !c("Treatment Method & Intensity", 
                                                                         "Treatment Intensity", 
                                                                         "Treatment Interval", 
                                                                         "Prescribed Fire",
                                                                         "Treatment Type") ]
  
  bycols = colnames(data)[!colnames(data) %in% c("V1", 
                                                 "value", 
                                                 "run", 
                                                 "Treatment Method & Intensity", 
                                                 "Treatment Intensity", 
                                                 "Treatment Interval", 
                                                 "Prescribed Fire",
                                                 "Treatment Type")]
  
  dt_merge = merge.data.table(trt, no_trt, by = bycols)
  
  if ("V1.x" %in% colnames(dt_merge)) {
    dt_merge$treat_chg = dt_merge$V1.x - dt_merge$V1.y
    colnames(dt_merge)[colnames(dt_merge) == "V1.x"] = "value"
    colnames(dt_merge)[colnames(dt_merge) == "V1.y"] = "notreat_value"
  } else if ("value.x" %in% colnames(dt_merge)) {
    dt_merge$treat_chg = dt_merge$value.x - dt_merge$value.y
  } else {
    cat("no V1 or value")
  }
  
  return(dt_merge)
}

# autogen subset
get_fixed_vars = function(DT) {
  isf = which(unname(sapply(DT, is.factor)))
  un = sapply(DT[ , ..isf ], FUN = unique )
  isun = sapply(un, FUN = function(X) length(X) == 1 )
  res = sapply(un[isun], as.character) 
  out = paste(names(res), res, collapse = " | ",sep = ":")
  return(out)
}
# group things
mean2sim = function(DT, scenarios) {
  DT_mn = DT[, .(mean(value)), by = c("run")]
  DT_mn = merge.data.table(DT_mn, scenarios, by = "run", all = FALSE)
  return(DT_mn)
}
mean2simCS = function(DT, scenarios) {
  DT_mn = DT[, .(mean(value)), by = c("run", "stratumID")]
  DT_mn = merge.data.table(DT_mn, scenarios, by = "run", all = FALSE)
  return(DT_mn)
}
# MONTHLY
mean2mn = function(DT, scenarios) {
  DT_mn = DT[, .(mean(value)), by = c("run", "sim_yr", "wym")]
  DT_mn = merge.data.table(DT_mn, scenarios, by = "run")
  DT_mn$yr_mn = (DT_mn$wym + (DT_mn$sim_yr - min(DT_mn$sim_yr)) * 12)
  return(DT_mn)
}
mean2mnCS = function(DT, scenarios) {
  DT_mn = DT[, .(mean(value)), by = c("run", "sim_yr", "wym", "stratumID")]
  DT_mn = merge.data.table(DT_mn, scenarios, by = "run")
  DT_mn$yr_mn = (DT_mn$wym + (DT_mn$sim_yr - min(DT_mn$sim_yr)) * 12)
  return(DT_mn)
}
# ANNUAL
mean2yr = function(DT, scenarios) {
  DT_yr = DT[, .(mean(value)), by = c("run", "sim_yr")]
  DT_yr = merge.data.table(DT_yr, scenarios, by = "run")
  return(DT_yr)
}
max2yr = function(DT, scenarios) {
  DT_yr = DT[, .(max(value)), by = c("run", "sim_yr")]
  DT_yr = merge.data.table(DT_yr, scenarios, by = "run", all = FALSE)
  return(DT_yr)
}
# CENTER OF MASS
CoM = function(in_data) {
  CoM =  match( TRUE, cumsum(in_data/sum(in_data)) > .5 ) - 1
  return(CoM)
}
CoM2yr = function(DT, scenarios) {
  DT_yr = DT[, .(CoM(value)), by = c("run", "sim_yr")]
  DT_yr = merge.data.table(DT_yr, scenarios, by = "run")
  return(DT_yr)
}

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}


# -=================================================-
# -================= Calc Fire Risk =================
# ==================================================-
# generate the fire risk output
# ijwf 2017 paper, assume max slope and wind
# relative_def = 1 - fire_et / pet
# litter loading litterC
# relaitve def w respect to litter loading

make_reldef = F
if (make_reldef) {
  path = "output/R_obj/fire_et.rdata"
  DT = read_robj(path)
  
  path = "output/R_obj/pet.rdata"
  DT2 = read_robj(path)
  
  fire = merge.data.table(x = DT, y = DT2[,c("run", "date", "value")], by = c("run", "date"))
  colnames(fire)[colnames(fire) == "value.x"] = "fire_et"
  colnames(fire)[colnames(fire) == "value.y"] = "pet"
  rm(DT, DT2)
  gc()
  
  DT_all = fire[, value := fifelse(pet > 0, (1 - fire_et/pet), 0)]
  DT_all = DT_all[, .(run, date, value)]

  rm(fire)
  gc()
  
  save(DT_all, file = "output/R_obj/reldef.rdata")
}

make_pspread = F
if (make_pspread) {
  
  path = "output/R_obj/litrc_sum.rdata"
  litrc = read_robj(path)[,c("run", "date", "value")]
  
  path = "output/R_obj/reldef.rdata"
  load(path)
  reldef = read_robj(path)[,c("run", "date", "value")]
  
  fire = merge.data.table(x = reldef, y = litrc, by = c("run", "date"))
  colnames(fire)[colnames(fire) == "value.x"] = "rel_def"
  colnames(fire)[colnames(fire) == "value.y"] = "litrc"
  
  rm(litrc, reldef)
  
  # pspread = product of constituent probs (litter load, relative deficit, topographic slope, wind dir)
  # litter load pl = 1 / 1 + e^ (- k1_l (1 - k2_l))
  k1_l = 3.9
  k2_l = 0.07
  # deficit pd = 1 / 1 + e^ (- k1_d (1 - k2_d))
  k1_d = 3.8
  k2_d = 0.27
  
  fire[, pdef := (1 / (1 + exp(1)^(-k1_d * (rel_def - k2_d)))  )]
  fire[, pload := (1 / (1 + exp(1)^(-k1_l * (litrc - k2_l)))  )]
  fire[, pspread := (pdef * pload)]
  
  DT_all = fire[, .(date, run, pspread)]
  colnames(DT_all)[colnames(DT_all) == "pspread"] = "value"
  
  rm(fire)
  
  save(DT_all, file = "output/R_obj/pspread.rdata")  
}

# -=================================================-
# -============== Random Forest gen =================
# -=================================================-
# output vars: fire et, height, lai, leafc, litrc_sum, litter.S, pet, plantc, psn, Qout, rz_storage
# random forests - loop through all vars
vars = c("plantc", "psn", "fire_et", "pspread")

genrf = T
if (genrf) {
  paths = file.path("output/R_obj/", paste0(vars, ".rdata"))
  pct_chg = list()
  
  for (i in seq_along(vars)) {
    cat(vars[i], ": reading data\n")
    read_dt = read_robj(paths[i])
    cat("NAs:",sum(is.na(read_dt$value)), "\n")
    if (vars[i] == "pspread") {
      read_dt = max2yr(read_dt, scen)
      colnames(read_dt)[colnames(read_dt) == "V1"] = "value"
    }
    
    cat(vars[i], ": avg data over sim\n")
    sim_avg = mean2sim(read_dt, scen)
    cat(vars[i], ": getting treatment change over sim\n")
    sim_chg = treat_change(sim_avg)
    cat(vars[i], ": getting pct change over sim\n")
    sim_chg[, pct_chg := treat_chg/value]
    
    # because it's already a percent change, dont need a %/%, thats just confusing
    if (vars[i] == "pspread") {
      sim_chg[, pct_chg := treat_chg]
    }
    pct_chg[[i]] = sim_chg
    cat(vars[i], ": calc random forest\n")
    rf = randomForest(x = sim_chg[,!c("treat_chg", "run", "value", "notreat_value", "Treatment Type", 
                                      "Treatment Intensity", "Prescribed Fire", "pct_chg")], y = sim_chg$treat_chg, localImp = TRUE)
    save(rf, file = paste0("output/random_forests/rf_", vars[i], ".rdata"))
  }
  
  save(pct_chg, file = "output/R_obj/pct_chg.rdata")
  
  rm(read_dt)
  gc()
  
  # fire severity, shrub and conifer seperate
  path = "output/R_obj/height.rdata"
  DT = read_robj(path)
  
  h_sim = mean2simCS(DT, scen)
  rm(DT); gc()
  
  h_shrub_chg = treat_change(h_sim[`Vegetation Type` == "shrub", ])
  h_shrub_chg[, pct_chg := treat_chg/value]
  h_sim_conifer_diff = merge.data.table(h_sim[stratumID == 11 & `Vegetation Type` == "conifer", ], h_sim[stratumID == 12 & `Vegetation Type` == "conifer", ], 
                                        by = colnames(h_sim)[!colnames(h_sim) %in% c("V1", "stratumID")] )
  h_sim_conifer_diff$V1 = h_sim_conifer_diff$V1.x - h_sim_conifer_diff$V1.y
  h_sim_conifer_diff[,c("V1.x", "V1.y","stratumID.y")] = NULL
  names(h_sim_conifer_diff)[names(h_sim_conifer_diff) == "stratumID.x"] = "stratumID"
  h_conifer_diffchg = treat_change(h_sim_conifer_diff)
  h_conifer_diffchg[, pct_chg := treat_chg/value]
  
  save(h_shrub_chg,file = "h_shrub_chg")
  save(h_conifer_diffchg, file = "h_conifer_diffchg")
  
  
  rf_h_shrub = randomForest(x = h_shrub_chg[,!c("treat_chg", "run", "value", "notreat_value", "Treatment Type", "Vegetation Type", 
                                    "Treatment Intensity", "Prescribed Fire", "pct_chg")], y = h_shrub_chg$treat_chg, localImp = TRUE)
  save(rf_h_shrub, file = "output/random_forests/rf_h_shrub.rdata")
  rm(rf_h_shrub) ; gc()
  
  rf_h_conifer_diff = randomForest(x = h_conifer_diffchg[,!c("treat_chg", "run", "value", "notreat_value", "Treatment Type", "Vegetation Type",
                                            "Treatment Intensity", "Prescribed Fire", "pct_chg")], y = h_conifer_diffchg$treat_chg, localImp = TRUE)
  save(rf_h_conifer_diff, file = "output/random_forests/rf_h_conifer_diff.rdata")
  rm(rf_h_conifer_diff) ; gc()
  
  
  beepr::beep(sound = 2)
}


# -=================================================-
# -================= Effect Sizes ===================
# ==================================================-

# plot effect sizes
plot_effect = T

if (plot_effect) {

  vartitle = c("Stand Carbon", "Net Primary Productivity", "Evapotranspiration", "Fire Spread Probability")
  vars = c("plantc", "psn", "fire_et", "pspread")
  load("output/R_obj/pct_chg.rdata")
  
  p_effect = list()
  for (i in 1:4) {
    p_effect[[i]] =  ggplot(pct_chg[[i]]) +
      aes(x = pct_chg, fill = `Treatment Type`) +
      geom_histogram(bins = 40L) +
      #labs(x = "Percent Change", y = "Scenario Count", title = vartitle[i]) +
      labs(title = vartitle[i]) +
      scale_x_continuous(labels = scales::percent) +
      scale_fill_hue(labels = c("Prescribed Fire", "Overstory Thinning", "Understory Thinning")) +
      theme_light() +
      theme(axis.title.x=element_blank(), axis.title.y=element_blank(),)
    
    if (vars[i] == "psn") {
      # quantile(pct_chg[[i]]$pct_chg, c(0.03, 0.97))
      p_effect[[i]] = p_effect[[i]] + scale_x_continuous(labels = scales::percent, limits = c(-1.5, .75) )
      mylegend<-g_legend(p_effect[[i]])
    }
    
  }
  
  p_ef = grid.arrange(arrangeGrob(p_effect[[1]] + theme(legend.position="none"), 
                                  p_effect[[2]] + theme(legend.position="none"),
                                  p_effect[[3]] + theme(legend.position="none"),
                                  p_effect[[4]] + theme(legend.position="none"), ncol = 2,
                                  bottom=textGrob("Percent Change", gp=gpar(fontsize=11)),
                                  left=textGrob("Scenario Count", gp=gpar(fontsize=11), rot=90)),
                      mylegend, ncol = 2, widths = c(8,2))
  
  
  ggsave(filename = "plots/PaperFigs/hist_effect.svg", p_ef, width = 9, height = 6)
  ggsave(filename = "plots/PaperFigs/hist_effect.tiff", p_ef, width = 9, height = 6, dpi = 300)
  beepr::beep(sound = 2) 
}

# -=================================================-
# -======= Fire severity - height chg/height =========
# ==================================================-

calcheights = F
if (calcheights) {
  path = "output/R_obj/height.rdata"
  DT = read_robj(path)
  
  # whole sim
  h_sim = mean2simCS(DT, scen)
  rm(DT); gc()
  
  h_shrub_chg = treat_change(h_sim[`Vegetation Type` == "shrub", ])
  h_shrub_chg[, pct_chg := treat_chg/value]
  
  h_sim_conifer_diff = merge.data.table(h_sim[stratumID == 11 & `Vegetation Type` == "conifer", ], h_sim[stratumID == 12 & `Vegetation Type` == "conifer", ], 
                                        by = colnames(h_sim)[!colnames(h_sim) %in% c("V1", "stratumID")] )
  h_sim_conifer_diff$V1 = h_sim_conifer_diff$V1.x - h_sim_conifer_diff$V1.y
  h_sim_conifer_diff[,c("V1.x", "V1.y","stratumID.y")] = NULL
  names(h_sim_conifer_diff)[names(h_sim_conifer_diff) == "stratumID.x"] = "stratumID"
  
  h_conifer_diffchg = treat_change(h_sim_conifer_diff)
  h_conifer_diffchg[, pct_chg := treat_chg/value]
  
  save(h_shrub_chg,file = "output/R_obj/h_shrub_chg.rdata")
  save(h_conifer_diffchg, file = "output/R_obj/h_conifer_diffchg.rdata")
}

load("output/R_obj/h_shrub_chg.rdata")
load("output/R_obj/h_conifer_diffchg.rdata")


# plot the changes
p_FE_conifer_type = ggplot(h_conifer_diffchg) +
  aes(x = pct_chg, fill = `Treatment Type`) +
  geom_histogram(bins = 40L) +
  labs(x = "Percent Change", subtitle = "Conifer-Shrub Mean Height Difference Change") +
  scale_x_continuous(labels = scales::percent) +
  scale_fill_hue(labels = c("Prescribed Fire", "Overstory Thinning", "Understory Thinning")) +
  theme_light() +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),)

g <- ggplot_build(p_FE_conifer_type)
fill_cols = unique(g$data[[1]]["fill"])[c(1,3),1]


p_FE_shrub_type = ggplot(h_shrub_chg) +
  aes(x = pct_chg, fill = `Treatment Type`) +
  geom_histogram(bins = 40L) +
  labs(x = "Percent Change", subtitle = "Shrub Mean Height Change") +
  scale_x_continuous(labels = scales::percent) +
  #scale_fill_hue(labels = c("Prescribed Fire", "Understory Thinning")) +
  scale_fill_manual(values = c(none = fill_cols[1], under = fill_cols[2]), labels = c("Prescribed Fire", "Understory Thinning"))+
  theme_light() +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),)

mylegend<-g_legend(p_FE_conifer_type)

p_ef2 = grid.arrange(arrangeGrob(p_FE_shrub_type + theme(legend.position="none"), 
                                 p_FE_conifer_type + theme(legend.position="none"), ncol = 2,
                                 top = textGrob("Fire Severity", gp=gpar(fontsize=16), x = 0.035, just = "left"),
                                 bottom=textGrob("Percent Change", gp=gpar(fontsize=11)),
                                 left=textGrob("Scenario Count", gp=gpar(fontsize=11), rot=90)),
                     mylegend, ncol = 2, widths = c(8,2))


ggsave(filename = "plots/PaperFigs/Fire_effects_type.svg", plot = p_ef2 , width = 12, height = 4)
ggsave(filename = "plots/PaperFigs/Fire_effects_type.tiff", plot = p_ef2, width = 12, height = 4, dpi = 300)


summary(h_conifer_diffchg$pct_chg)
summary(h_shrub_chg$pct_chg)

# -=================================================-
# -============ Effect Size Stats ===================
# ==================================================-

effect_stats = T
if (effect_stats) {
  
  options(scipen=999)
  
  pctgt0 = function(X) sum(X > 0)/length(X)
  pctlt0 = function(X) sum(X < 0)/length(X)
  
  vartitle = c("Stand Carbon", "Net Primary Productivity", "Evapotranspiration", "Fire Spread Probability")
  load("output/R_obj/pct_chg.rdata")
  
  
  stats = data.frame(vartitle)
  for (i in 1:4) {
    stats[i,2:7] = summary(pct_chg[[i]]$pct_chg)
    stats$pctgt0[i] = pctgt0(pct_chg[[i]]$pct_chg)
    stats$pctlt0[i] = pctlt0(pct_chg[[i]]$pct_chg)
  }
  colnames(stats)[2:7] = names(summary(pct_chg[[1]]$pct_chg))
  
  
    
}

# -=================================================-
# -============ Effect Size Scatter Plots ===================
# ==================================================-

effect_scatter = T
if (effect_scatter) {
  
  options(scipen=999)
  
  vartitle = c("Stand Carbon", "Net Primary Productivity", "Evapotranspiration", "Fire Spread Probability")
  load("output/R_obj/pct_chg.rdata")
  
  pct_chg_DT = data.table(pct_chg[[1]], pct_chg[[2]][,14:16], pct_chg[[3]][,14:16], pct_chg[[4]][,14:16])
  names(pct_chg_DT)[14:16] = paste0(names(pct_chg_DT)[14:16], "_plantc")
  names(pct_chg_DT)[17:19] = paste0(names(pct_chg_DT)[17:19], "_npp")
  names(pct_chg_DT)[20:22] = paste0(names(pct_chg_DT)[20:22], "_et")
  names(pct_chg_DT)[23:25] = paste0(names(pct_chg_DT)[23:25], "_pspread")
  
  
  library(dplyr)
  library(ggplot2)
  
  plantc_x_et = pct_chg_DT %>%
    ggplot() +
    aes(x = pct_chg_plantc, y = pct_chg_et, colour = `Treatment Type`) +
    geom_point(size = 0.5) +
    scale_x_continuous(labels = scales::percent) +
    scale_y_continuous(labels = scales::percent) +
    labs(x = "Percent Change Stand Carbon", y = "Percent Change ET", title = "Percent Change Stand Carbon vs Evapotranspiration") +
    scale_colour_discrete(labels = c("Prescribed Fire", "Overstory Thinning", "Understory Thinning")) +
    theme_light() 
  
  npp_x_et = pct_chg_DT %>%
    filter(pct_chg_npp >= -1 & pct_chg_npp <= 1) %>%
    ggplot() +
    aes(x = pct_chg_npp, y = pct_chg_et, colour = `Treatment Type`) +
    geom_point(size = 0.5) +
    scale_x_continuous(labels = scales::percent) +
    scale_y_continuous(labels = scales::percent) +
    labs(x = "Percent Change NPP", y = "Percent Change ET", title = "Percent Change Net Primary Productivity vs Evapotranspiration") +
    scale_colour_discrete(labels = c("Prescribed Fire", "Overstory Thinning", "Understory Thinning")) +
    theme_light()

  ggsave(filename = "plots/PaperFigs/Effect_Interactions.svg", plot = grid.arrange(plantc_x_et,npp_x_et, ncol = 2),
         width = 12, height = 4)
  
  
  }


# -=================================================-
# -============= Effect interactions ==============
# ==================================================-

effect_int = T
if (effect_int) {
  
  options(scipen=999)
  
  vartitle = c("Stand Carbon", "Net Primary Productivity", "Evapotranspiration", "Fire Spread Probability")
  load("output/R_obj/pct_chg.rdata")
  
  plantc = pct_chg[[1]]
  npp = pct_chg[[2]]
  et = pct_chg[[3]]
  pspread = pct_chg[[4]]
  
  
  p_npp1 = npp %>%
    filter(pct_chg >= -1 & pct_chg <= 1) %>%
    ggplot() +
    aes(x = `Soil Depth`, y = pct_chg, fill = `Treatment Type`) +
    geom_boxplot() +
    theme_light() +
    scale_fill_hue(labels = c("Prescribed Fire", "Overstory Thinning", "Understory Thinning")) +
    scale_y_continuous(labels = scales::percent) +
    scale_x_discrete(labels = c("High", "Medium", "Low")) +
    labs(y = "Percent Change Net Primary Productivity", x = "Soil Water Capcity")
  
  #, title = "Percent Change in Net Primary Productivity by Treatment Type and Soils"
  
  # p_et1 = ggplot(et) +
  #   aes(x = `Vegetation Type`, y = pct_chg, fill = `Treatment Type`) +
  #   geom_boxplot() +
  #   scale_fill_hue(labels = c("Prescribed Fire", "Overstory Thinning", "Understory Thinning")) +
  #   scale_y_continuous(labels = scales::percent) +
  #   labs(y = "Percent Change ET", title = "Percent Change in Evapotranspiration by Treatment Type and Vegetation Type") +
  #   theme_light()
  # 
  
  ggsave(filename = "plots/PaperFigs/npp_veg_x_soil.svg", p_npp1, width = 9, height = 6)
  ggsave(filename = "plots/PaperFigs/npp_veg_x_soil.tiff", p_npp1, width = 9, height = 6, dpi = 300)
  
}






# -=================================================-
# -============= Min Depth gen + plots ==============
# ==================================================-
# get minimum depth and produce the plots
ex_rf = T
if (ex_rf) {
  
  p = list()
  vars = c("plantc", "psn", "fire_et", "pspread")
  
  for (i in seq_along(vars)) {
    
    load(paste0("output/random_forests/rf_", vars[i], ".rdata"))
    
    cat(vars[i], ": calc min depth\n")
    mindep = min_depth_distribution(rf)
    save(mindep, file = paste0("output/random_forests/mindep_", vars[i], ".rdata"))
    
    cat(vars[i], ": make min depth plot\n")
    p[[i]] = plot_min_depth_distribution(mindep)
    
  }
  
  save(p, file = "output/random_forests/plots_mindep.rdata")
  
  
  load("output/random_forests/rf_h_shrub.rdata")
  mindep_h_shrub = min_depth_distribution(rf_h_shrub)
  save(mindep_h_shrub, file = "output/random_forests/mindep_h_shrub.rdata")
  p_h_shrub = plot_min_depth_distribution(mindep_h_shrub)
  save(p_h_shrub, file = "output/random_forests/plots_mindep_h_shrub.rdata")
  
  load("output/random_forests/rf_h_conifer_diff.rdata")
  mindep_h_conifer_diff = min_depth_distribution(rf_h_conifer_diff)
  save(mindep_h_conifer_diff, file = "output/random_forests/mindep_h_conifer_diff.rdata")
  p_h_conifer_diff = plot_min_depth_distribution(mindep_h_conifer_diff)
  save(p_h_conifer_diff, file = "output/random_forests/plots_mindep_h_conifer.rdata")
  
  beepr::beep(sound = 2)
}

# i guess doesn't work
# explain_forest(rf, interactions = TRUE)

# format plots
mindep_plots = T
if (mindep_plots) {
  #load("output/random_forests/plots_mindep.rdata")
  vartitle = c("Stand Carbon", "Net Primary Productivity", "Evapotranspiration", "Fire Spread Probability")
  vars = c("plantc", "psn", "fire_et", "pspread")
  
  
  
  for (i in seq_along(vars)) {
    
    load(file = paste0("output/random_forests/mindep_", vars[i], ".rdata"))
    
    p[[i]] = plot_min_depth_distribution(mindep)
    
    p[[i]] = p[[i]] + 
      ggtitle(vartitle[i]) +
      theme(axis.title.x=element_blank(), axis.title.y=element_blank())
    
  }
  
  mylegend<-g_legend(p[[1]])
  
  p_mndp1 = grid.arrange(arrangeGrob(p[[1]] + theme(legend.position="none"), 
                                     p[[2]] + theme(legend.position="none"),
                                     p[[3]] + theme(legend.position="none"),
                                     p[[4]] + theme(legend.position="none"), ncol = 2,
                                     bottom=textGrob("Number of Trees", gp=gpar(fontsize=11))),
                                     mylegend, ncol = 2, widths = c(8,2))
  
  ggsave(filename = "plots/PaperFigs/mindep.svg", plot = p_mndp1, width = 12, height = 9)
  ggsave(filename = "plots/PaperFigs/mindep.tiff", plot = p_mndp1, width = 9, height = 6, dpi = 300)
  
  
  load(file = "output/random_forests/mindep_h_shrub.rdata")
  load(file = "output/random_forests/mindep_h_conifer_diff.rdata")
  p_h_shrub = plot_min_depth_distribution(mindep_h_shrub)
  p_h_conifer_diff = plot_min_depth_distribution(mindep_h_conifer_diff)
  
  p_h_shrub = p_h_shrub +  
    ggtitle("Shrub Height") +
    theme(axis.title.x=element_blank(), axis.title.y=element_blank())
  
  p_h_conifer_diff = p_h_conifer_diff +  
    ggtitle("Conifer - Shrub Height Difference") +
    theme(axis.title.x=element_blank(), axis.title.y=element_blank())
  
  mylegend = g_legend(p_h_conifer_diff)
  
  p_mndp2 = grid.arrange(arrangeGrob(p_h_shrub + theme(legend.position="none"), 
                                     p_h_conifer_diff + theme(legend.position="none"), ncol = 2,
                                     bottom=textGrob("Number of Trees", gp=gpar(fontsize=11))),
                         mylegend, ncol = 2, widths = c(8,2))
  
  ggsave(filename = "plots/PaperFigs/MinDep_FE.svg", plot = p_mndp2, width = 9, height = 3)
  ggsave(filename = "plots/PaperFigs/MinDep_FE.tiff", plot = p_mndp2, width = 9, height = 3, dpi = 300)
  
  beepr::beep(2)
  
}



# -=================================================-
# -======== Lineplot - plantc : Stand Carbon ========
# ==================================================-

path = "output/R_obj/plantc.rdata"
DT = read_robj(path)

plantc_mn = mean2mn(DT, scen)
rm(DT); gc()

plantc_mn_selected = plantc_mn[(`Treatment Method & Intensity` == "no treatment" | 
                                  `Treatment Method & Intensity` == "under_0.4 1mo_rm100litter" & 
                                  #`Treatment Interval` %in% c("5", "10", "20")) &
                                  `Treatment Interval` %in% c("10")) &
                                  `Climatic Period` == "wet" & 
                                 (`Soil Depth` %in% c("deep")) & 
                                 (`Root Sharing Coefficient` %in% c("0.5"))  &
                                 `Vegetation Type` == "shrub" & 
                                 Aspect == "N" & 
                                 `Climate Change` == FALSE, ]

p_pc = ggplot(plantc_mn_selected) +
  aes(x = yr_mn, y = V1, color = `Treatment Method & Intensity`) +
  geom_line(size = 2L) +
  labs(x = "Time Since Start (Months)", y = ~ "Stand Carbon " (kgC/m^2), title = "Monthly Stand Carbon", color = "Treatment Type",
       subtitle = "Shrub vegetation treated at 10 year intervals") +
  #     subtitle = get_fixed_vars(plantc_mn_selected), 
  geom_vline(xintercept = seq.int(7, 12*30, 12*10)) +
  scale_color_hue(labels = c("No Treatment", "Understory Treatment")) +
  scale_y_continuous(limits = c(0, 6) ) +
  theme_light() 
  #facet_wrap(`Soil Depth`~`Climatic Period`)


ggsave(filename = "plots/lineplot/standc_shrub.svg", plot = p_pc, width = 9, height = 6)
ggsave(filename = "plots/lineplot/standc_shrub.tiff", plot = p_pc, width = 9, height = 6, dpi = 300)


trt = mean(plantc_mn_selected$V1[plantc_mn_selected$`Treatment Method & Intensity` == "under_0.4 1mo_rm100litter"])
notrt = mean(plantc_mn_selected$V1[plantc_mn_selected$`Treatment Method & Intensity` == "no treatment"])
trtdiff = notrt - trt
trtdiffpct = trtdiff/notrt

# -=================================================-
# -================= Precip x NPP ===================
# ==================================================-

path = "output/R_obj/psn.rdata"
npp_yr_avg = mean2yr(read_robj(path), scen)

npp_yr_chg = treat_change2(npp_yr_avg)
# there's an outlier so get rid of that later
# npp_out = npp_yr_chg[ abs(npp_yr_chg$treat_chg) < 15, ]


# scatterplot - Annual (wy) NPP change vs precip
clim = as.data.table(read_clim("clim/upGGmod.rain"))

pcp_sum_yr = clim[, .(pcp_sum = sum(rain)), by = c("wy")]

start_yr = c(1977, 1943, 2007)
names(start_yr) = c("wet","variable","dry")

npp_yr_chg[, year := sim_yr + unname(start_yr[npp_yr_chg$`Climatic Period`]) ]

npp_yr_chg = merge.data.table(npp_yr_chg, pcp_sum_yr, by.x = "year", by.y = "wy")

npp_yr_chg = npp_yr_chg[ abs(npp_yr_chg$treat_chg) < 15, ]

p_npp =  ggplot(npp_yr_chg) +
  geom_point(aes(x = pcp_sum, y = treat_chg)) +
  #labs(x = "Percent Change", title = vartitle[i]) +
  theme_light()

ggsave(filename = "plots/npp_pcpsum_scatter.svg", p_npp, width = 16, height = 9)

# -=================================================-
# -===== debug - patch & strata height plantc =======
# ==================================================-

# -- NO SHADE DEBUG --
# get the data
read_in = fread("output/noshade/conifer/height") 

# set colnames to just the nums
colnames(read_in)[startsWith(colnames(read_in), "run")] = scen_noshade$run

# add dates to data
dt = add_dates(read_in)

DT_all = melt(data = dt, measure.vars = which(!colnames(dt) %in% c("day", "month", "year", "basinID", "hillID", "zoneID", "patchID", "date", "stratumID", "wy", "yd")),  variable.name = "run", variable.factor = FALSE)
DT_all[, run := as.integer(run), ]

days = length(DT_all$run[DT_all$run == 1])
DT_all[, sim_day :=  rep(c(1:days), nrow(DT_all)/days), ]
DT_all[, mn := lubridate::month(date), ]
DT_all[, yr := lubridate::year(date), ]
DT_all[, wy := data.table::fifelse(mn >= 10, yr + 1, yr), ]
DT_all[, wym := data.table::fifelse(mn >= 10, mn - 9, mn + 3), ]
DT_all[, sim_yr := wy - min(DT_all$yr), ]

hgt_allp_mn = DT_all[, .(mean(value)), by = c("run", "sim_yr","patchID", "stratumID", "wym")]
hgt_allp_mn = merge.data.table(hgt_allp_mn, scen_noshade, by = "run")
hgt_allp_mn$yr_mn = (hgt_allp_mn$wym + (hgt_allp_mn$sim_yr - min(hgt_allp_mn$sim_yr)) * 12)


# -- INDIVIDUAL RUN DEBUG --
# read from full runs
rh_in = readin_rhessys_output("output/debug/thinx1_overstory_fulltime/conifer")
hgt_debug = rh_in$cd[ , .(height, date, patchID, stratumID, month, wy)]
hgt_debug[, wym := data.table::fifelse(month >= 10, month - 9, month + 3), ]
hgt_debug[, sim_yr := wy - min(hgt_debug$wy), ]
hgt_debug_mn = hgt_debug[, .(mean(height)), by = c("wym", "wy","sim_yr", "patchID", "stratumID")]
hgt_debug_mn$yr_mn = (hgt_debug_mn$wym + (hgt_debug_mn$sim_yr - min(hgt_debug_mn$sim_yr)) * 12)

p_hgt_debug = ggplot(hgt_debug_mn[1:1000,]) +
  aes(x = yr_mn, y = V1, color = as.factor(patchID), linetype = as.factor(stratumID)) +
  geom_line(size = 1L) +
  labs(x = "Month", y = "Height", title = "Monthly Height - No Shading", 
       subtitle = "Conifer, 1943 - 2010, deep soils, 0.5 root sharing, 1 treat at start, vertical lines marking treatments") +
  geom_vline(xintercept = c(1)) +
  #ylim(0, 4.5) +
  theme_light()
  #facet_wrap(vars(`Treatment Method & Intensity`))

ggsave(filename = "plots/heights_debug_allp_longseries.svg", plot = p_hgt_debug , width = 16, height = 9)


# general version
rh_in = readin_rhessys_output("output/debug/rh_out_2020-04-14_23.36.11/conifer")

debug = rh_in$cd[ , .(height, date, patchID, stratumID, month, wy)]
debug[, wym := data.table::fifelse(month >= 10, month - 9, month + 3), ]
debug[, sim_yr := wy - min(hgt_debug$wy), ]

debug_mn = debug[, .(mean(height)), by = c("wym", "wy","sim_yr", "patchID", "stratumID")]

debug_mn$yr_mn = (debug_mn$wym + (debug_mn$sim_yr - min(debug_mn$sim_yr)) * 12)

p_debug = ggplot(debug_mn[1:1000,]) +
  aes(x = yr_mn, y = V1, color = as.factor(patchID), linetype = as.factor(stratumID)) +
  geom_line(size = 1L) +
  labs(x = "Month", y = "C", title = "Monthly Height - No Shading", 
       subtitle = "Conifer, 1943 - 2010, deep soils, 0.5 root sharing, 1 treat at start, vertical lines marking treatments") +
  geom_vline(xintercept = c(1)) +
  #ylim(0, 4.5) +
  theme_light()


ggsave(filename = "", plot = p_debug , width = 16, height = 9)


# just explore the data
pdg = rh_in$pdg
pd = rh_in$pd
cdg = rh_in$cdg
cd = rh_in$cd


p1 = ggplot(cdg) +
  aes(x = date, y = plantc, colour = as.factor(stratumID)) +
  geom_line(size = 1L) +
  scale_color_hue() +
  theme_light() +
  ggtitle("Plant Carbon") +
  facet_wrap(vars(patchID))

p2 = ggplot(cd) +
  aes(x = date, y = height, colour = as.factor(stratumID)) +
  geom_line(size = 1L) +
  scale_color_hue() +
  theme_light() +
  ggtitle("Height") +
  facet_wrap(vars(patchID))

p3 = ggplot(cdg) +
  aes(x = date, y = woodc, colour = as.factor(stratumID)) +
  geom_line(size = 1L) +
  scale_color_hue() +
  theme_light() +
  ggtitle("Wood Carbon") +
  facet_wrap(vars(patchID))

p4 = ggplot(cdg) +
  aes(x = date, y = live_stemc, colour = as.factor(stratumID)) +
  geom_line(size = 1L) +
  scale_color_hue() +
  theme_light() +
  ggtitle("live_stemc") +
  facet_wrap(vars(patchID))

p5 = ggplot(cdg) +
  aes(x = date, y = dead_stemc, colour = as.factor(stratumID)) +
  geom_line(size = 1L) +
  scale_color_hue() +
  theme_light() +
  ggtitle("dead_stemc") +
  facet_wrap(vars(patchID))

p6 = ggplot(cdg) +
  aes(x = date, y = live_crootc, colour = as.factor(stratumID)) +
  geom_line(size = 1L) +
  scale_color_hue() +
  theme_light() +
  ggtitle("live_crootc") +
  facet_wrap(vars(patchID))

p7 = ggplot(cdg) +
  aes(x = date, y = dead_crootc, colour = as.factor(stratumID)) +
  geom_line(size = 1L) +
  scale_color_hue() +
  theme_light() +
  ggtitle("dead_crootc") +
  facet_wrap(vars(patchID))

p8 = ggplot(cdg) +
  aes(x = date, y = dead_stemn, colour = as.factor(stratumID)) +
  geom_line(size = 1L) +
  scale_color_hue() +
  theme_light() +
  ggtitle("dead_stemn") +
  facet_wrap(vars(patchID))

p9 = ggplot(cdg) +
  aes(x = date, y = cwdc, colour = as.factor(stratumID)) +
  geom_line(size = 1L) +
  scale_color_hue() +
  theme_light() +
  ggtitle("cwdc") +
  facet_wrap(vars(patchID))

p10 = ggplot(cdg) +
  aes(x = date, y = root_depth, colour = as.factor(stratumID)) +
  geom_line(size = 1L) +
  scale_color_hue() +
  theme_light() +
  ggtitle("root_depth") +
  facet_wrap(vars(patchID))

p11 = ggplot(cdg) +
  aes(x = date, y = leafc, colour = as.factor(stratumID)) +
  geom_line(size = 1L) +
  scale_color_hue() +
  theme_light() +
  ggtitle("leafc") +
  facet_wrap(vars(patchID))

p12 = ggplot(cdg) +
  aes(x = date, y = dead_leafc, colour = as.factor(stratumID)) +
  geom_line(size = 1L) +
  scale_color_hue() +
  theme_light() +
  ggtitle("dead_leafc") +
  facet_wrap(vars(patchID))

p13 = ggplot(cdg) +
  aes(x = date, y = leafc_age1, colour = as.factor(stratumID)) +
  geom_line(size = 1L) +
  scale_color_hue() +
  theme_light() +
  ggtitle("leafc_age1") +
  facet_wrap(vars(patchID))

p14 = ggplot(cdg) +
  aes(x = date, y = leafc_age2, colour = as.factor(stratumID)) +
  geom_line(size = 1L) +
  scale_color_hue() +
  theme_light() +
  ggtitle("leafc_age2") +
  facet_wrap(vars(patchID))

p15 = ggplot(cdg) +
  aes(x = date, y = leafc_store, colour = as.factor(stratumID)) +
  geom_line(size = 1L) +
  scale_color_hue() +
  theme_light() +
  ggtitle("leafc_store") +
  facet_wrap(vars(patchID))

pdf("plots/debug/plantc_crash4.16.20.pdf")
p1; p2; p3; p4; p5; p6; p7; p8; p9; p10; p11; p12; p13; p14; p15
dev.off()


# ==================================================-
# -====================== END =======================
# -=================================================-
