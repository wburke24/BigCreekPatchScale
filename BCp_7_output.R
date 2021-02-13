# -=================================================-
# -============= Output_BigCreek_Patch ==============
# ==================================================-

setwdhere()

library(data.table)
library(tidyverse)
library(plotly)
library(chron)
library(RHESSysIOinR)
library(gridExtra)
library(scales)
library(grid)
library(cowplot)
library(randomForestExplainer)

load("scenarios.rdata")
load("scen.rdata")

# -=================================================-
# -===== Helper Functions =================
# ==================================================-

scen = plyr::rename(scen, c("start_period" = "Aridity",
                            "treatments" = "Method & Intensity",
                            "interval" = "Treatment Interval",
                            "soils" = "Soil Water Capacity",
                            "sharing" = "Root Sharing",
                            "treat_type" = "Treatment Type",
                            "treat_intensity" = "Treatment Intensity",
                            "presc_fire" = "Prescribed Fire",
                            "veg" = "Vegetation Type",
                            "aspect" = "Aspect",
                            "clim_chg" = "Climate Warming"))

scenarios = plyr::rename(scenarios, c("start_period" = "Aridity",
                                      "treatments" = "Method & Intensity",
                                      "interval" = "Treatment Interval",
                                      "soils" = "Soil Water Capacity",
                                      "sharing" = "Root Sharing",
                                      "treat_type" = "Treatment Type",
                                      "treat_intensity" = "Treatment Intensity",
                                      "presc_fire" = "Prescribed Fire",
                                      "veg" = "Vegetation Type",
                                      "aspect" = "Aspect",
                                      "climate" = "Climate Warming"))

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
  dt_merge = merge.data.table(data[data$"Method & Intensity" != "no treatment",],
                              data[data$"Method & Intensity" == "no treatment", 
                                   c("Aridity", 
                                     "Soil Water Capacity",
                                     "Root Sharing",
                                     "Vegetation Type", 
                                     "Aspect", 
                                     "Climate Warming",
                                     "V1")],
                              by = c("Aridity", 
                                     "Soil Water Capacity",
                                     "Root Sharing",
                                     "Vegetation Type", 
                                     "Aspect", 
                                     "Climate Warming"))
  
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

# autogen subset and output text
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
# -===== Gen Stats: Fire Risk, Random Forest ======
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

genpctchg = F
if (genpctchg) {
  vars = c("plantc", "psn", "fire_et", "pspread")
  paths = file.path("output/R_obj/", paste0(vars, ".rdata"))
  pct_chg = list()
  
  for (i in seq_along(vars)) {
    cat(vars[i], ": reading data\n")
    read_dt = read_robj(paths[i])
    cat("NAs:",sum(is.na(read_dt$value)), "\n")
    if (vars[i] == "pspread") {
      # get annual max
      read_dt = max2yr(read_dt, scen)
      colnames(read_dt)[colnames(read_dt) == "V1"] = "value"
    }
    
    cat(vars[i], ": avg data over sim\n")
    sim_avg = mean2sim(read_dt, scen)
    cat(vars[i], ": getting treatment change over sim\n")
    sim_chg = treat_change(sim_avg)
    cat(vars[i], ": getting pct change over sim\n")
    sim_chg[, pct_chg := treat_chg/value]
    pct_chg[[i]] = sim_chg
  }
  
  rm(read_dt)
  gc()
  
  # fire severity, shrub and conifer seperate
  path = "output/R_obj/height.rdata"
  DT = read_robj(path)
  
  h_sim = mean2simCS(DT, scen)
  rm(DT); gc()
  
  h_shrub_chg = treat_change(h_sim[`Vegetation Type` == "shrub", ])
  h_shrub_chg[, pct_chg := treat_chg/value]
  pct_chg[[5]] = h_shrub_chg
  
  h_sim_conifer_diff = merge.data.table(h_sim[stratumID == 11 & `Vegetation Type` == "conifer", ], h_sim[stratumID == 12 & `Vegetation Type` == "conifer", ], 
                                        by = colnames(h_sim)[!colnames(h_sim) %in% c("V1", "stratumID")] )
  h_sim_conifer_diff$V1 = h_sim_conifer_diff$V1.x - h_sim_conifer_diff$V1.y
  h_sim_conifer_diff[,c("V1.x", "V1.y","stratumID.y")] = NULL
  names(h_sim_conifer_diff)[names(h_sim_conifer_diff) == "stratumID.x"] = "stratumID"
  h_conifer_diffchg = treat_change(h_sim_conifer_diff)
  h_conifer_diffchg[, pct_chg := treat_chg/value]
  pct_chg[[6]] = h_conifer_diffchg
  
  names(pct_chg) = c("plantc", "psn", "fire_et", "pspread","severity_shrub", "severity_conifer")
  
  save(pct_chg, file = "output/R_obj/pct_chg.rdata")
  
  beepr::beep(sound = 2)
}

genrf = F
if (genrf) {
  
  library(randomForest)
  library(randomForestExplainer)
  
  vars = c("plantc", "psn", "fire_et", "pspread")
  paths = file.path("output/R_obj/", paste0(vars, ".rdata"))
  load("output/R_obj/pct_chg.rdata")
  
  for (i in seq_along(vars)) {
    cat(vars[i], ": calc random forest\n")
    rf = randomForest(x = sim_chg[,!c("treat_chg", "run", "value", "notreat_value", "Treatment Type", 
                                      "Treatment Intensity", "Prescribed Fire", "pct_chg")], y = sim_chg$treat_chg, localImp = TRUE)
    save(rf, file = paste0("output/random_forests/rf_", vars[i], ".rdata"))
  }
  gc()

  beepr::beep(sound = 2)
}

# -=================================================-
# -===== Fig 3 Monthly Stand Carbon ===========
# ==================================================-

plot_standc = T
if (plot_standc) {
  
  genplantc = F
  if (genplantc) {
    path = "output/R_obj/plantc.rdata"
    DT = read_robj(path)
    plantc_mn = mean2mn(DT, scen)
    rm(DT); gc()
    save(plantc_mn, file = "output/R_obj/plantc_mn.rdata")
  }
  
  load("output/R_obj/plantc_mn.rdata")
  
  plantc_mn_selected = plantc_mn[(`Method & Intensity` == "no treatment" | 
                                    `Method & Intensity` %in% c("under_0.4 1mo_rm100litter", "over_0.4 carbon_removed") & 
                                    #`Treatment Interval` %in% c("5", "10", "20")) &
                                    `Treatment Interval` %in% c("10")) &
                                   `Aridity` == "wet" & 
                                   (`Soil Water Capacity` %in% c("shallow")) & 
                                   (`Root Sharing` %in% c("0.5"))  &
                                   `Vegetation Type` == "conifer" & 
                                   Aspect == "N" & 
                                   `Climate Warming` == FALSE, ]
  
  fig3 = ggplot(plantc_mn_selected) +
    aes(x = yr_mn, y = V1, color = `Method & Intensity`) +
    geom_line(size = 1L,position = position_dodge(width = 1)) +
    labs(x = "Time Since Start (Months)", y = ~ "Stand Carbon " (kgC/m^2), title = "Monthly Total Stand Carbon", color = "Treatment Type",
         subtitle = "Vegetation treated at 10 year intervals") +
    #     subtitle = get_fixed_vars(plantc_mn_selected), 
    geom_vline(xintercept = seq.int(7, 12*30, 12*10), alpha = 0.5) +
    scale_color_hue(labels = c("No Treatment","Overstory Treatment", "Understory Treatment")) +
    scale_y_continuous(limits = c(0, 3.5) ) +
    theme_cowplot()
  #facet_wrap(`Soil Depth`~`Climatic Period`)
  
  ggsave(filename = "plots/Final_figs/f3_standc_mn.svg", plot = fig3, width = 9, height = 6)
  ggsave(filename = "plots/Final_figs/f3_standc_mn.tiff", plot = fig3, width = 9, height = 6, dpi = 300)
  ggsave(filename = "plots/Final_figs/f3_standc_mn_wide.tiff", plot = fig3, width = 10, height = 6, dpi = 300, compression = "lzw")
  ggsave(filename = "plots/Final_figs/f3_standc_mn_wide2.tiff", plot = fig3, width = 11, height = 6, dpi = 300, compression = "lzw")
  
}

# -=================================================-
# -===== Fig 4 Effect Sizes ===================
# ==================================================-

# plot effect sizes
plot_effect = T
if (plot_effect) {

  vartitle = c("Stand Carbon", "Net Primary Productivity", "Evapotranspiration", "Fire Spread Probability", "Shrub Fuel Height", "Conifer Canopy Fuel Gap")
  vars = c("plantc", "psn", "fire_et", "pspread", "severity_shrub", "severity_conifer")
  load("output/R_obj/pct_chg.rdata")
  
  p_effect = list()
  for (i in seq_along(vars)) {
    
    pfun = function(X) {
      Y = X / length(pct_chg[[i]]$pct_chg)
      Y = round(Y,digits = 2)
      return(Y)
    }

    p_effect[[i]] =  ggplot(pct_chg[[i]]) +
      aes(x = pct_chg, color = `Treatment Type`) +
      geom_freqpoly(bins = 40L, position = "dodge", size = 1) +
      labs(title = vartitle[i]) +
      scale_x_continuous(labels = scales::percent) +
      scale_color_hue(labels = c("Prescribed Fire", "Overstory Thinning", "Understory Thinning")) +
      scale_y_continuous(labels = pfun) +
      theme_minimal_hgrid() +
      theme(axis.title.x=element_blank(), axis.title.y=element_blank(), plot.title = element_text(size = 13))
    
    if (vars[i] == "psn") {
      g <- ggplot_build(p_effect[[i]])
      color_cols = unique(g$data[[1]]["colour"])[c(1,3),1]
      p_effect[[i]] = p_effect[[i]] + scale_x_continuous(labels = scales::percent, limits = c(-1.5, .75) )
      
      mylegend<-g_legend(p_effect[[i]])
    }
    if (vars[i] == "severity_shrub") {
      p_effect[[i]] = p_effect[[i]] + scale_color_manual(values = c(none = color_cols[1], under = color_cols[2]), 
                                                         labels = c("Prescribed Fire", "Understory Thinning"))
    }
    
    p_effect[[i]] = p_effect[[i]] + theme(legend.position="none")
    
  }

  fig4 = plot_grid(plotlist = p_effect, labels = "AUTO", label_size = 12)
  fig4 = plot_grid(fig4, mylegend, rel_widths = c(3, .75))
  fig4 = grid.arrange(arrangeGrob(fig4,  
                                  left=textGrob("Proportion of Scenarios", gp=gpar(fontsize=11), rot=90),
                                  bottom=textGrob("Percent Change", gp=gpar(fontsize=11))))
                                 
  ggsave(filename = "plots/Final_figs/f4_hist_effects.svg", fig4, width = 9, height = 6, scale = 1.2)
  ggsave(filename = "plots/Final_figs/f4_hist_effects.tiff", fig4, width = 9, height = 6, dpi = 300, scale = 1.2)
  
  #beepr::beep(sound = 2)
  
  #ggplotly(p_effect[[6]])
  
}

plot_effect_veg = F
if (plot_effect_veg) {
  
  vartitle = c("Stand Carbon", "Net Primary Productivity", "Evapotranspiration", "Fire Intensity Index", "Fire Severity (Shrub)", "Fire Severity (Conifer+Shrub)")
  vars = c("plantc", "psn", "fire_et", "pspread", "severity_shrub", "severity_conifer")
  load("output/R_obj/pct_chg.rdata")
  
  p_effect = list()
  for (i in seq_along(vars)) {
    p_effect[[i]] =  ggplot(pct_chg[[i]]) +
      aes(x = pct_chg, color = `Vegetation Type`) +
      geom_freqpoly(bins = 40L, position = "dodge", size = 1) +
      labs(subtitle = vartitle[i]) +
      scale_x_continuous(labels = scales::percent) +
      scale_fill_hue(labels = c("Prescribed Fire", "Overstory Thinning", "Understory Thinning")) +
      theme_light() +
      theme(axis.title.x=element_blank(), axis.title.y=element_blank(),)
    
    if (vars[i] == "psn") {
      g <- ggplot_build(p_effect[[i]])
      color_cols = unique(g$data[[1]]["colour"])[c(1,3),1]
      
      p_effect[[i]] = p_effect[[i]] + scale_x_continuous(labels = scales::percent, limits = c(-1.5, .75) )
      
      mylegend<-g_legend(p_effect[[i]])
    }
    
    # if (vars[i] == "severity_shrub") {
    #   p_effect[[i]] = p_effect[[i]] + scale_color_manual(values = c(none = color_cols[1], under = color_cols[2]), 
    #                                                      labels = c("Prescribed Fire", "Understory Thinning"))
    # }
    
  }
  
  p_ef_veg = grid.arrange(arrangeGrob(p_effect[[1]] + theme(legend.position="none"), 
                                  p_effect[[2]] + theme(legend.position="none"),
                                  p_effect[[3]] + theme(legend.position="none"),
                                  p_effect[[4]] + theme(legend.position="none"),
                                  p_effect[[5]] + theme(legend.position="none"),
                                  p_effect[[6]] + theme(legend.position="none"),
                                  ncol = 3,
                                  bottom=textGrob("Percent Change", gp=gpar(fontsize=11)),
                                  left=textGrob("Scenario Count", gp=gpar(fontsize=11), rot=90)),
                      mylegend, ncol = 2, widths = c(8,2))
  
  ggsave(filename = "plots/PaperFigs/hist_effectsizes_veg.svg", p_ef_veg, width = 9, height = 6)
  ggsave(filename = "plots/PaperFigs/hist_effectsizes_veg.tiff", p_ef_veg, width = 9, height = 6, dpi = 300)
  
  beepr::beep(sound = 2) 
  
}

# -=================================================-
# -===== Fig 5 + 6 Effect interactions ============
# ==================================================-

effect_int = T
if (effect_int) {
  options(scipen=999)
  
  load("output/R_obj/pct_chg.rdata")
  #[`Vegetation Type` == "conifer"]
  npp = pct_chg$psn
  et = pct_chg$fire_et
  sev_shrub = pct_chg$severity_shrub
  sev_conifer = pct_chg$severity_conifer
  fsp = pct_chg$pspread
  
  # FIG 5 but with conifer only, and fsp instead of shrub sev
  f5_plots = list()
  f5_plots[[1]] = npp %>%
    filter(`Vegetation Type` == "conifer") %>%
    filter(pct_chg >= -1 & pct_chg <= 1) %>%
    ggplot() + aes(x = `Soil Water Capacity`, y = pct_chg, fill = `Treatment Type`) + geom_boxplot() + theme_cowplot() +
    scale_fill_hue(labels = c("Prescribed Fire", "Overstory Thinning", "Understory Thinning")) + 
    scale_y_continuous(labels = scales::percent) + scale_x_discrete(labels = c("High", "Medium", "Low")) +
    labs(y = "Change NPP", title = "Net Primary Productivity") + 
    theme(axis.title.x=element_blank(), plot.title = element_text(size = 16))
  mylegend<-g_legend(f5_plots[[1]])
  f5_plots[[1]] = f5_plots[[1]] + theme(legend.position="none")
  
  f5_plots[[2]] = et %>%
    filter(`Vegetation Type` == "conifer") %>%
    ggplot() + aes(x = `Soil Water Capacity`, y = pct_chg, fill = `Treatment Type`) + geom_boxplot() + theme_cowplot() +
    scale_fill_hue(labels = c("Prescribed Fire", "Overstory Thinning", "Understory Thinning")) +
    scale_y_continuous(labels = scales::percent) +  scale_x_discrete(labels = c("High", "Medium", "Low")) +
    labs(y = "Change in ET", title = "Evapotranspiration") + 
    theme(axis.title.x=element_blank(), legend.position="none",, plot.title = element_text(size = 16))
  f5_plots[[3]] = sev_conifer %>%
    ggplot() +  aes(x = `Soil Water Capacity`, y = pct_chg, fill = `Treatment Type`) + geom_boxplot() + theme_cowplot() +
    scale_fill_hue(labels = c("Prescribed Fire", "Overstory Thinning", "Understory Thinning")) +
    scale_y_continuous(labels = scales::percent) + scale_x_discrete(labels = c("High", "Medium", "Low")) +
    labs(y = "Change in Conifer Canopy Fuel Gap", title = "Conifer Canopy Fuel Gap") +
    theme(axis.title.x=element_blank(), legend.position="none", plot.title = element_text(size = 16))
  f5_plots[[4]] = fsp %>% filter(`Vegetation Type` == "conifer") %>%
    ggplot() + aes(x = `Soil Water Capacity`, y = pct_chg, fill = `Treatment Type`) + geom_boxplot() + theme_cowplot() +
    scale_fill_hue(labels = c("Prescribed Fire", "Overstory Thinning", "Understory Thinning")) +
    scale_y_continuous(labels = scales::percent) + scale_x_discrete(labels = c("High", "Medium", "Low")) +
    labs(y = "Change in Fire Spread Probability", title = "Fire Spread Probability") +
    theme(axis.title.x=element_blank(), legend.position="none", plot.title = element_text(size = 16))
  
  fig5 = plot_grid(plotlist = f5_plots, labels = "AUTO", label_size = 12)
  fig5 = plot_grid(fig5, mylegend, rel_widths = c(3, .75))
  fig5 = grid.arrange(arrangeGrob(fig5, bottom=textGrob("Plant Accessible Water Storage Capacity", gp=gpar(fontsize=11))))
  
  # ggplotly(f5_plots[[1]])
  # ggplotly(f5_plots[[2]])
  # ggplotly(f5_plots[[3]])
  # ggplotly(f5_plots[[4]])

  # FIG 6 conifer severity
  f6_plots = list()
  f6_plots[[1]] = sev_conifer %>%
    ggplot() + aes(x = `Climate Warming`, y = pct_chg, fill = `Treatment Type`) +  geom_boxplot() + theme_cowplot() +
    scale_fill_hue(labels = c("Prescribed Fire", "Overstory Thinning", "Understory Thinning")) +
    scale_y_continuous(labels = scales::percent) + scale_x_discrete(labels = c("Baseline", "+2C")) +
    labs( x = "Climate Warming", title = "Climate Warming")+ 
    theme(axis.title.y=element_blank(), plot.title = element_text(size = 16))
  mylegend<-g_legend(f6_plots[[1]])
  f6_plots[[1]] = f6_plots[[1]] + theme(legend.position="none")
  
  f6_plots[[2]] = sev_conifer %>%
    ggplot() + aes(x = `Aridity`, y = pct_chg, fill = `Treatment Type`) + geom_boxplot() + theme_cowplot() +
    scale_fill_hue(labels = c("Prescribed Fire", "Overstory Thinning", "Understory Thinning")) +
    scale_y_continuous(labels = scales::percent) + scale_x_discrete(labels = c("Dry", "Variable", "Wet")) +
    labs( x = "Aridity", title = "Aridity") + 
    theme(axis.title.y=element_blank(), legend.position="none", plot.title = element_text(size = 16))
  f6_plots[[3]] = sev_conifer %>%
    ggplot() + aes(x = `Treatment Intensity`, y = pct_chg, fill = `Treatment Type`) + geom_boxplot() + theme_cowplot() +
    scale_fill_hue(labels = c("Prescribed Fire", "Overstory Thinning", "Understory Thinning")) +
    scale_y_continuous(labels = scales::percent) + labs( x = "Treatment Intensity", title = "Treatment Intensity")+
    theme(axis.title.y=element_blank(), legend.position="none", plot.title = element_text(size = 16))
  f6_plots[[4]] = sev_conifer %>%
    ggplot() + aes(x = reorder(`Treatment Interval`, pct_chg), y = pct_chg, fill = `Treatment Type`) + geom_boxplot() + 
    theme_cowplot() + scale_fill_hue(labels = c("Prescribed Fire", "Overstory Thinning", "Understory Thinning")) +
    scale_y_continuous(labels = scales::percent) + scale_x_discrete(labels = c("5", "10", "30")) +
    labs( x = "Treatment Interval", title = "Treatment Interval")+ 
    theme(axis.title.y=element_blank(), legend.position="none", plot.title = element_text(size = 16))
  
  ggplotly(f6_plots[[1]])
  ggplotly(f6_plots[[2]])
  ggplotly(f6_plots[[3]])
  ggplotly(f6_plots[[4]])
  
  fig6 = plot_grid(plotlist = f6_plots, labels = "AUTO", label_size = 12)
  fig6 = plot_grid(fig6, mylegend, rel_widths = c(3, .75))
  # f6_title = ggdraw() + draw_label("Conifer Canopy Fuel Gap", fontface = 'bold', x = 0, hjust = 0, size = 16) +
  #   theme(plot.margin = margin(0, 0, 0, 7))
  # fig6 = plot_grid(f6_title, fig6, ncol = 1, rel_heights = c(0.05, 1))
  fig6 = grid.arrange(arrangeGrob(fig6, left=textGrob("Change in Conifer Canopy Fuel Gap", gp=gpar(fontsize=11), rot=90)))
  
  # fig 5
  # scaling might need to be different for svg vs tiff
  ggsave(filename = "plots/Final_figs/f5_boxplots1.svg", fig5, width = 12, height = 8)
  ggsave(filename = "plots/Final_figs/f5_boxplots1.tiff", fig5, width = 12, height = 8, dpi = 300)
  
  # fig 6
  ggsave(filename = "plots/Final_figs/f6_boxplots2.svg", fig6, width = 12, height = 8)
  ggsave(filename = "plots/Final_figs/f6_boxplots2.tiff", fig6, width = 12, height = 8, dpi = 300)
  
}

# -=================================================-
# -===== Fig 7 Min Depth gen + plots =========
# ==================================================-
calc_mindep = F
if (calc_mindep) {
  
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
  mindep = min_depth_distribution(rf_h_shrub)
  save(mindep, file = "output/random_forests/mindep_h_shrub.rdata")
  p_h_shrub = plot_min_depth_distribution(mindep)
  save(p_h_shrub, file = "output/random_forests/plots_mindep_h_shrub.rdata")
  
  load("output/random_forests/rf_h_conifer_diff.rdata")
  mindep = min_depth_distribution(rf_h_conifer_diff)
  save(mindep, file = "output/random_forests/mindep_h_conifer_diff.rdata")
  p_h_conifer_diff = plot_min_depth_distribution(mindep)
  save(p_h_conifer_diff, file = "output/random_forests/plots_mindep_h_conifer.rdata")
  
  beepr::beep(sound = 2)
}
# i guess doesn't work
# explain_forest(rf, interactions = TRUE)

# format plots
mindep_plots = T
if (mindep_plots) {
  #load("output/random_forests/plots_mindep.rdata")
  vartitle = c("Stand Carbon",
                "Net Primary Productivity",
                "Evapotranspiration",
                "Fire Spread Probability",
                "Shrub Fuel Height",
                "Conifer Canopy Fuel Gap")
  
  vars = c("plantc", "psn", "fire_et", "pspread", "h_shrub", "h_conifer_diff")
  p=list()
  
  for (i in seq_along(vars)) {
    load(file = paste0("output/random_forests/mindep_", vars[i], ".rdata"))
    
    mindep$variable = gsub("Treatment Method & Intensity",  "Method & \nIntensity", mindep$variable)
    mindep$variable = gsub("Soil Depth", "PAWSC", mindep$variable)
    mindep$variable = gsub("Root Sharing Coefficient","Root\nSharing", mindep$variable)
    mindep$variable = gsub("Climatic Period", "Aridity", mindep$variable)
    mindep$variable = gsub("Treatment Interval", "Interval", mindep$variable)
    mindep$variable = gsub("Climate Change", "Climate\nWarming", mindep$variable)
    mindep$variable = gsub("Vegetation Type", "Veg Type", mindep$variable)
    # "Aspect"
    
    p[[i]] = plot_min_depth_distribution(mindep)
    p[[i]] = p[[i]] + 
      theme_cowplot() +
      ggtitle(vartitle[i]) +
      theme(axis.title.x=element_blank(), axis.title.y=element_blank(), plot.title = element_text(size = 18))
    # title is bigger since the whole plot is larger
    if (i == 1) {
      mylegend<-g_legend(p[[1]])
    }
    p[[i]] = p[[i]] + theme(legend.position="none")
    
    test = p[[1]] + theme()
    geom_label(label.size = 0, )
    test$layers[[2]]$aes_params$size = 1
    test$theme
  }
  
  fig7 = plot_grid(plotlist = p, labels = "AUTO", label_size = 12)
  fig7 = plot_grid(fig7, mylegend, rel_widths = c(3, .5))
  fig7 = grid.arrange(arrangeGrob(fig7, bottom=textGrob("Number of Trees", gp=gpar(fontsize=11))))
  
  ggsave(filename = "plots/Final_figs/f7_mindep.svg", plot = fig7, width = 12, height = 9)
  ggsave(filename = "plots/Final_figs/f7_mindep_nocomp.tiff", plot = fig7, width = 16, height = 12, dpi = 300, scale = 1)
  ggsave(filename = "plots/Final_figs/f7_mindep.tiff", plot = fig7, width = 16, height = 12, dpi = 300, scale = 1, compression = "lzw")
  #beepr::beep(2)
}

# -=================================================-
# -===== Output  ===================
# ==================================================-

make_pdf = F
if (make_pdf) {
  pdf(file = "plots/figs3-7.pdf")
  fig3
  plot(fig4)
  plot(fig5)
  plot(fig6)
  plot(fig7)
  dev.off()
}

# -=================================================-
# -===== Effect Size Stats ===================
# ==================================================-

effect_stats = F
if (effect_stats) {
  
  options(scipen=999)
  
  pctgt0 = function(X) sum(X > 0)/length(X)
  pctlt0 = function(X) sum(X < 0)/length(X)
  
  vartitle = c("Stand Carbon", "Net Primary Productivity", "Evapotranspiration", "Fire Intensity Index")
  load("output/R_obj/pct_chg.rdata")
  
  
  stats = data.frame(vartitle)
  for (i in 1:4) {
    stats[i,2:7] = summary(pct_chg[[i]]$pct_chg)
    stats$pctgt0[i] = pctgt0(pct_chg[[i]]$pct_chg)
    stats$pctlt0[i] = pctlt0(pct_chg[[i]]$pct_chg)
  }
  colnames(stats)[2:7] = names(summary(pct_chg[[1]]$pct_chg))
}

old_code = T
if (old_code) {
  # NPP, ET, both fire severities by PAWSC and treat type
  f5_plots = list()
  f5_plots[[1]] = npp %>%
    filter(pct_chg >= -1 & pct_chg <= 1) %>%
    ggplot() + aes(x = `Soil Water Capacity`, y = pct_chg, fill = `Treatment Type`) + geom_boxplot() + theme_cowplot() +
    scale_fill_hue(labels = c("Prescribed Fire", "Overstory Thinning", "Understory Thinning")) + 
    scale_y_continuous(labels = scales::percent) + scale_x_discrete(labels = c("High", "Medium", "Low")) +
    labs(y = "Change NPP", title = "Net Primary Productivity") + theme(axis.title.x=element_blank())
  mylegend<-g_legend(f5_plots[[1]])
  f5_plots[[1]] = f5_plots[[1]] + theme(legend.position="none")
  
  f5_plots[[2]] = et %>%
    ggplot() + aes(x = `Soil Water Capacity`, y = pct_chg, fill = `Treatment Type`) + geom_boxplot() + theme_cowplot() +
    scale_fill_hue(labels = c("Prescribed Fire", "Overstory Thinning", "Understory Thinning")) +
    scale_y_continuous(labels = scales::percent) +  scale_x_discrete(labels = c("High", "Medium", "Low")) +
    labs(y = "Change in ET", title = "Evapotranspiration") + theme(axis.title.x=element_blank(), legend.position="none")
  f5_plots[[3]] = sev_conifer %>%
    ggplot() +  aes(x = `Soil Water Capacity`, y = pct_chg, fill = `Treatment Type`) + geom_boxplot() + theme_cowplot() +
    scale_fill_hue(labels = c("Prescribed Fire", "Overstory Thinning", "Understory Thinning")) +
    scale_y_continuous(labels = scales::percent) + scale_x_discrete(labels = c("High", "Medium", "Low")) +
    labs(y = "Change in Conifer Canopy Fuel Gap", title = "Conifer Canopy Fuel Gap") +
    theme(axis.title.x=element_blank(), legend.position="none")
  g <- ggplot_build(f5_plots[[3]])
  fill_cols = unique(g$data[[1]]["fill"])[c(1,3),1]
  f5_plots[[4]] = sev_shrub %>%
    ggplot() + aes(x = `Soil Water Capacity`, y = pct_chg, fill = `Treatment Type`) + geom_boxplot() + theme_cowplot() +
    scale_fill_hue(labels = c("Prescribed Fire", "Overstory Thinning", "Understory Thinning")) +
    scale_y_continuous(labels = scales::percent) + scale_x_discrete(labels = c("High", "Medium", "Low")) +
    scale_fill_manual(values = c(none = fill_cols[1], under = fill_cols[2]), labels = c("Prescribed Fire", "Understory Thinning"))+
    labs(y = "Change in Shrub Fuel Height",  title = "Shrub Fuel Height") + 
    theme(axis.title.x=element_blank(), legend.position="none")
  
  fig5 = plot_grid(plotlist = f5_plots, labels = "AUTO", label_size = 12)
  fig5 = plot_grid(fig5, mylegend, rel_widths = c(3, .75))
  fig5 = grid.arrange(arrangeGrob(fig5, bottom=textGrob("Plant Accessible Water Storage Capacity", gp=gpar(fontsize=11))))
}
# ==================================================-
# -====================== END =======================
# -=================================================-
