getwd()
setwd("F:/R/Meta-analysis-crop/MBC_analysis")

library(readxl)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(ggpmisc)
library(stringr)
library(forcats)
library(metagear)
library(Matrix)
library(metadat)
library(numDeriv)
library(metafor)
library(cowplot)
library(gridExtra)
library(maps)
library(ggpolypath)
library(patchwork)
library(orchaRd)
library(caret)
library(randomForest)
library(patchwork)

# data with extracted agroecological zones https://gaez.fao.org

dat_sta_fin <- read.csv("F:/R/Meta-analysis-crop/Soil C stock following perennialization new.csv", na = "NA")


### impute missing SD ###
# https://cran.r-project.org/web/packages/metagear/metagear.pdf

impute_missingness(dat_sta_fin)

# creating ids of the levels for the random component

dat_sta_fin <- dat_sta_fin %>%
  mutate(
    SOC.annual.SD_imputed = if_else(is.na(SOC.annual.SD), "yes", "no"),
    effect_size_id = row_number(),
    site_lat.log.elevation = interaction(lat.dec, long.dec, elevation)
  ) %>%
  group_by(site_lat.log.elevation) %>%
  mutate(site_id = cur_group_id()) %>%
  ungroup() %>%
  group_by(author) %>%
  mutate(study_id = cur_group_id()) %>%
  ungroup()

# soil pH, crop age of perennial
dat_sta_fin <- dat_sta_fin %>% mutate(
  soil.pH.classification = case_when(
    soil.pH < 6.5 ~ "<6.5",
    soil.pH > 7.5 ~ ">7.5",
    TRUE ~ "6.5-7.5"
  ),
  crop.age.classification = case_when(
    age.peren < 5 ~ "Age<5",
    age.peren >= 5 & age.peren <= 10 ~ "Age 5-10",
    TRUE ~ "Age>10"
  ),
  climatic.zone_1 = sub(",.*$","", climatic.zone),
  textural.class = case_when(
    textural.class == "Silt loam" ~ "Silty loam",
    TRUE ~ textural.class
  )
  
)

set.seed(1234)

dat_sta_fin <- impute_SD(as.data.frame(dat_sta_fin),
                                 columnSDnames = c("SOC.annual.SD", "SOC.peren.SD"), columnXnames = c("soc.annual.mean_30cm", "soc.peren.mean_30cm"),
                                 method = "HotDeck",
                                 range = 10, M = 1
)

### effect sizes ###
# 1. Log Ratio of Means (ROM)

dat_sta_fin_SOC <- escalc(
  measure = "ROM", n1i = SOC.peren.n, n2i = SOC.annual.n, m1i = soc.peren.mean_30cm, m2i = soc.annual.mean_30cm,
  sd1i = SOC.peren.SD, sd2i = SOC.annual.SD, data = dat_sta_fin, var.names = c("yi", "vi")
)

main.model_imp_SOC <- rma.mv(yi = yi, V = vi, random = ~ 1 | site_id / effect_size_id, data = dat_sta_fin_SOC, method = "REML", test = "t")
summary(main.model_imp_SOC)
main.model_imp_SOC$beta # effect size

# 添加标准误（SE）列
dat_sta_fin_SOC$sei <- sqrt(dat_sta_fin_SOC$vi)

# 精度为 1/SE
dat_sta_fin_SOC$precision <- 1 / dat_sta_fin_SOC$sei

# 如果已有样本量列（SOC.peren.n + SOC.annual.n），可以这样创建总样本量：
dat_sta_fin_SOC$sample_size <- dat_sta_fin_SOC$SOC.peren.n + dat_sta_fin_SOC$SOC.annual.n

# 图 (a)：效应值 vs 样本量
Sample_size_1 <- ggplot(dat_sta_fin_SOC, aes(x = sample_size, y = yi)) +
  geom_point(aes(size = precision), alpha = 0.4, color = "grey30") +
  geom_smooth(method = "lm", se = TRUE, fill = "grey", color = "black") +
  scale_size_continuous(name = "Precision (1/SE)") +
  labs(x = "Sample size (number of comparisons)", y="") +
  theme(aspect.ratio = 0.9,
        panel.background = element_rect(fill = 'white', colour = 'white'),
        plot.title = element_text(hjust = 0.5, size=20, margin = ggplot2::margin(0,0,0.5,0,'cm')),
        axis.title=element_text(colour='black', size=14),
        axis.text.x = element_text(colour = 'black', size = 12,margin = ggplot2::margin(0,0,0.3,0,'cm')),
        axis.text.y = element_text(colour = 'black', size = 12,margin = ggplot2::margin(0,0,0,0.3,'cm')),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.line = element_line(colour = 'black',size=1),
        legend.position = c(0.8,0.25))

# 图 (b)：效应值 vs 标准误
SE_1 <- ggplot(dat_sta_fin_SOC, aes(x = sei, y = yi)) +
  geom_point(aes(size = precision), alpha = 0.4, color = "grey30") +
  geom_smooth(method = "lm", se = TRUE, fill = "grey", color = "black") +
  scale_size_continuous(name = "Precision (1/SE)") +
  labs(x = "Standard error", y="") +
  theme(aspect.ratio = 0.9,
        panel.background = element_rect(fill = 'white', colour = 'white'),
        plot.title = element_text(hjust = 0.5, size=20, margin = ggplot2::margin(0,0,0.5,0,'cm')),
        axis.title=element_text(colour='black', size=14),
        axis.text.x = element_text(colour = 'black', size = 12,margin = ggplot2::margin(0,0,0.3,0,'cm')),
        axis.text.y = element_text(colour = 'black', size = 12,margin = ggplot2::margin(0,0,0,0.3,'cm')),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.line = element_line(colour = 'black',size=1),
        legend.position = c(0.8,0.25))
  

filename <- paste('Result_Regression/',"SOC_funnel",'.jpeg', sep='')
jpeg(filename,width = 3000,height=1500,res = 300)
SOC_funnel <- Sample_size_1 + SE_1
print(SOC_funnel)
graphics.off()

summary(lm(yi ~ sample_size, data = dat_sta_fin_SOC))
summary(lm(yi ~ sei, data = dat_sta_fin_SOC))



# checks
par(mfrow=c(1,1))
forest(main.model_imp_SOC, header = TRUE, addpred = TRUE, showweights = TRUE, cex = 0.75)
profile(main.model_imp_SOC)

qqnorm(residuals(main.model_imp_SOC), main = "QQ plot: residuals")
qqline(residuals(main.model_imp_SOC), col = "red")

# variance components
main.model_imp_SOC$sigma2
# total heterogeneity
sum(main.model_imp_SOC$sigma2)

# obtaining I^2 (relative measurement of heterogeneity)
i2_ml(main.model_imp_SOC, method = c("ratio"))
i2_ml(main.model_imp_SOC, method = c("matrix"))

main.model_imp_SOC_I2 <- i2_ml(main.model_imp_SOC)
main.model_imp_SOC_I2 * 100 # (% of heterogeneity)


tiff("F:/R/Meta-analysis-crop/MBC_analysis/control_SOC_main_model.tiff", height = 6, width = 12, units = "cm", compression = "lzw", res = 600)
main.plot_SOC <- orchard_plot(main.model_imp_SOC,mod = "1",group="study.number",xlab = "Effect size (lnRR)", transfm = "none", g = FALSE) +
  scale_colour_manual(
    values = "#a55213",
    aesthetics = c("colour", "fill")
  ) +
  labs(y = "LnRR-Soil Organic Carbon") +
  theme(legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.position = c(1, -0.03),
    panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank()
  ) +
  labs(subtitle = "Perennial vs SOC")
main.plot_SOC$layers[[4]] <- NULL
summary_SOC <- summary(main.model_imp_SOC)
estimate <- summary_SOC$beta[1]    # 点的数值
lowerCL <- summary_SOC$ci.lb       # 置信区间下限
upperCL <- summary_SOC$ci.ub       # 置信区间上限
point_data <- data.frame(
  condition = "Intrcpt",       # x轴标签
  estimate = estimate,
  lowerCL = lowerCL,
  upperCL = upperCL          # 颜色值（和之前orchard一致）
)
main.plot_SOC <- main.plot_SOC +
  geom_point(
    data = point_data,  # 你自定义的数据
    aes(x = condition, y = estimate),
    shape = 21,
    size = 3,    
    fill = "#a55213",  # ✅ 直接指定颜色，不映射
    color = "black"
  )
main.plot_SOC
dev.off()

#percentage of change for the overall effect and CIs
100 * (exp(as.numeric(main.model_imp_SOC$beta)) - 1)
100 * (exp(as.numeric(main.model_imp_SOC$ci.lb)) - 1)
100 * (exp(as.numeric(main.model_imp_SOC$ci.ub)) - 1)

# 敏感性分析代码

dat_sta_fin_SOC_sen <- dat_sta_fin %>%
  filter(SOC.annual.SD_imputed == "no")

# 进行随机效应模型
dat_sta_fin_SOC_sen1 <- escalc(
  measure = "ROM", n1i = SOC.peren.n, n2i = SOC.annual.n, m1i = soc.peren.mean_30cm, m2i = soc.annual.mean_30cm,
  sd1i = SOC.peren.SD, sd2i = SOC.annual.SD, data = dat_sta_fin_SOC_sen, var.names = c("yi", "vi")
)

main.model_imp_SOC_sen <- rma.mv(yi = yi, V = vi, random = ~ 1 | site_id / effect_size_id, data = dat_sta_fin_SOC_sen1, method = "REML", test = "t")

tiff("F:/R/Meta-analysis-crop/MBC_analysis/control_SOC_main_sen_model.tiff", height =8, width = 12, units = "cm", compression = "lzw", res = 600)
main.plot_SOC_sen <- orchard_plot(main.model_imp_SOC_sen,mod = "1",group="study.number",xlab = "Effect size (lnRR)", transfm = "none") +
  scale_colour_manual(
    values = "#a55213",
    aesthetics = c("colour", "fill")
  ) +
  labs(y = "LnRR-Soil Organic Carbon") +
  theme(legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.position = c(1, 0.01),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank()
  ) +
  labs(subtitle = "Perennial vs SOC")
main.plot_SOC_sen
dev.off()

#########################################################################
##function to convert to percentage of change (100 * (exp(estimate) - 1))
#########################################################################

perc_change <- function(model_moderator, moderator_variable){
  
  # 使用 mod_results() 一次性提取预测和估计值
  model_moderator_pre_mod <- mod_results(model_moderator, mod = moderator_variable, group = "site_id")
  model_moderator_pre_mod <- model_moderator_pre_mod$mod_table
  
  # 计算百分比变化
  model_moderator_pre_mod <- model_moderator_pre_mod %>%
    mutate(
      estimate_perc = 100 * (exp(estimate) - 1),
      lowerCL_perc = 100 * (exp(lowerCL) - 1),
      upperCL_perc = 100 * (exp(upperCL) - 1),
      lowerPR_perc = 100 * (exp(lowerPR) - 1),
      upperPR_perc = 100 * (exp(upperPR) - 1)
    ) %>%
    rename(
      yi = estimate_perc,
      moderator = name
    )
  
  return(model_moderator_pre_mod)
}

##################################################################
# first categorical variables with MBC
##################################################################


# climatic zone
main.model_SOC_climatic_zone_1 <- rma.mv(yi = yi, V = vi, mod = ~ climatic.zone_1 - 1, random = ~ 1 | site_id / effect_size_id, data = dat_sta_fin_SOC, method = "REML", test = "t")
main.model_SOC_climatic_zone_1
anova(main.model_SOC_climatic_zone_1)

lab <- c(
  "Climatic.zone_1Subtropics" = "Subtropics",
  "Climatic.zone_1Temperate" = "Temperate",
  "Climatic.zone_1Tropics" = "Tropics"
)

tiff("F:/R/Meta-analysis-crop/MBC_analysis/control_monoculture_main_model_moderator_climate_zone_SOC.tiff", height = 7, width = 10, units = "cm", compression = "lzw", res = 600)
main.plot_SOC_climatic_zone <- orchard_plot(main.model_SOC_climatic_zone_1, mod = "climatic.zone_1",group="study.number",xlab = "Effect size (lnRR)", transfm = "none") +
  scale_x_discrete(limits = c("Temperate", "Subtropics", "Tropics"))+
  scale_y_continuous(limits = c(-1.6, 1.8),breaks = seq(-1, 1.5, 1))+
  scale_colour_manual(values = c("#36c179", "#3ed739", "#225a27"),aesthetics = c("colour", "fill")) + 
  labs(y = "") +  
  theme(
    axis.text.y = element_text(hjust = 0, angle = 0, colour = "black"),
    legend.direction = "horizontal", panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(),
    legend.position="none"
  ) + 
  annotate("text", x = 0.55, y = 1.2, label = "F = 11.7577, p < 0.0001", size = 2.5, fontface = "bold")
main.plot_SOC_climatic_zone
dev.off()

#using the function to calculate percentage of change
main.model_SOC_climatic_zone_1_perc <- perc_change(main.plot_SOC_climatic_zone ,"climatic.zone_1")
main.model_SOC_climatic_zone_1_perc
write.csv(main.model_MBC_climatic_zone_1_perc, file = "climatic_zone_1_perc.csv", row.names = FALSE)


# crop.age.classification
##############################

dat_sta_fin_SOC$crop.age.classification <- as.factor(dat_sta_fin_SOC$crop.age.classification)
dat_sta_fin_SOC$crop.age.classification <- factor(dat_sta_fin_SOC$crop.age.classification, levels = c("Age>10", "Age 5-10", "Age<5"))


main.model_SOC_age_peren <- rma.mv(yi = yi, V = vi, mod = ~ crop.age.classification - 1, random = ~ 1 | site_id / effect_size_id, data = dat_sta_fin_SOC, method = "REML", test = "t")
main.model_SOC_age_peren
anova(main.model_SOC_age_peren)

lab <- c(
  "Crop.age.classification>10" = "Age>10",
  "Crop.age.classification5-10" = "Age 5-10",
  "Crop.age.classification<5" = "Age<5"
)

tiff("F:/R/Meta-analysis-crop/MBC_analysis/control_monoculture_main_model_moderator_age_perennial_SOC.tiff", height = 9, width = 10, units = "cm", compression = "lzw", res = 600)
main.plot_SOC_age_peren <- orchard_plot(main.model_SOC_age_peren, mod = "crop.age.classification",group="study.number", xlab = "Effect size (lnRR)", transfm = "none") +
  scale_colour_manual(values = c("#7e1d45","#c05579","#e6778b"), aesthetics = c("colour", "fill")) + 
  scale_y_continuous(limits = c(-1.6, 1.8),breaks = seq(-1, 1.5, 1)) +
  labs(y = "") + 
  theme(
    axis.text.y = element_text(hjust = 0, angle = 0, colour = "black"),
    legend.direction = "horizontal", panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(),
    legend.position="bottom"
  ) + 
  annotate("text", x = 0.55, y = 1.2, label = "F = 30.6877, p < 0.0001", size = 2.5, fontface = "bold")
main.plot_SOC_age_peren
dev.off()

# pH area
######################################
dat_sta_fin_SOC$soil.pH.classification <- as.factor(dat_sta_fin_SOC$soil.pH.classification)
dat_sta_fin_SOC$soil.pH.classification <- factor(dat_sta_fin_SOC$soil.pH.classification, levels = c(">7.5", "6.5-7.5", "<6.5"))

main.model_SOC_soil_pH <- rma.mv(yi = yi, V = vi, mod = ~ soil.pH.classification - 1, random = ~ 1 | site_id / effect_size_id, data = dat_sta_fin_SOC, method = "REML", test = "t")
main.model_SOC_soil_pH
anova(main.model_SOC_soil_pH)

lab <- c(
  "Soil.pH.classification<6.5" = "< 6.5",
  "Soil.pH.classification>7.5" = "> 7.5",
  "Soil.pH.classification6.5-7.5" = "6.5 - 7.5"
)

tiff("F:/R/Meta-analysis-crop/MBC_analysis/control_monoculture_main_model_moderator_soil_pH_SOC.tiff", height = 7, width = 10, units = "cm", compression = "lzw", res = 600)
main.plot_SOC_soil_pH <- orchard_plot(main.model_SOC_soil_pH,mod = "soil.pH.classification",group="study.number",xlab = "Effect size (lnRR)", transfm = "none") +
  scale_y_continuous(limits = c(-1.6, 1.8),breaks = seq(-1, 1.5, 1))+
  scale_colour_manual(values = c("#69176d","#c04ab5","#dc8ce5"), aesthetics = c("colour", "fill")) + 
  labs(y = "") +
  theme(
    axis.text.y = element_text(hjust = 0, angle = 0, colour = "black"),
    legend.direction = "horizontal", panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(),
    legend.position="none"
  )+
  annotate("text", x = 0.55, y = 1.2, label = "F = 11.6453, p < 0.0001", size = 2.5, fontface = "bold")
main.plot_SOC_soil_pH
dev.off()

#####Plant Functional Type - Perennial
######################################

dat_sta_fin_SOC$PFT.peren <- as.factor(dat_sta_fin_SOC$PFT.peren)
dat_sta_fin_SOC$PFT.peren <- factor(dat_sta_fin_SOC$PFT.peren, levels = c("C4-C3", "C4", "C3"))

main.model_SOC_PFT_peren <- rma.mv(yi = yi, V = vi, mod = ~ PFT.peren - 1, random = ~ 1 | site_id / effect_size_id, data = dat_sta_fin_SOC, method = "REML", test = "t")
main.model_SOC_PFT_peren
anova(main.model_SOC_PFT_peren)

lab <- c(
  "PFT.perenC3" = "C3",
  "PFT.perenC4" = "C4",
  "PFT.perenC4-C3" = "C3 - C4"
)


tiff("F:/R/Meta-analysis-crop/MBC_analysis/control_monoculture_main_model_moderator_PFT_perennial_SOC.tiff", height = 7, width = 10, units = "cm", compression = "lzw", res = 600)
main.plot_SOC_PFT_peren <- orchard_plot(main.model_SOC_PFT_peren, mod = "PFT.peren",group="study.number",xlab = "Effect size (lnRR)", transfm = "none") +
  scale_colour_manual(values = c("#1d3683","#3b67b3","#62a5d8"), aesthetics = c("colour", "fill")) + 
  scale_y_continuous(limits = c(-1.6, 1.8),breaks = seq(-1, 1.5, 1)) +
  labs(y = "") + 
  theme(
    axis.text.y = element_text(hjust = 0, angle = 0, colour = "black"),
    legend.direction = "horizontal", panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(),
    legend.position="none"
  ) + 
  annotate("text", x = 0.55, y = 1.2, label = "F = 11.0547, p < 0.0001", size = 2.5, fontface = "bold")
main.plot_SOC_PFT_peren
dev.off()

# leg.peren#####2
################################################################################

main.model_SOC_leg_peren <- rma.mv(yi = yi, V = vi, mod = ~ leg.peren - 1, random = ~ 1 | site_id / effect_size_id, data = dat_sta_fin_SOC, method = "REML", test = "t")
main.model_SOC_leg_peren
anova(main.model_SOC_leg_peren)

lab <- c(
  "No" = "Non-legume",
  "Yes" = "Legume "
)

tiff("F:/R/Meta-analysis-crop/MBC_analysis/control_monoculture_main_model_moderator_legume_perennial_SOC.tiff", height = 7.5, width = 11, units = "cm", compression = "lzw", res = 600)
main.plot_SOC_leg_peren <- orchard_plot(main.model_SOC_leg_peren, mod = "leg.peren",group="study.number", xlab = "Effect size (lnRR)", transfm = "none") +
  scale_x_discrete(labels = c("Yes"="Legume","No"="Non-legume"))+
  scale_colour_manual(values = c("#e2a06b","#c16528"), aesthetics = c("colour", "fill")) + 
  scale_y_continuous(limits = c(-1.6, 1.8),breaks = seq(-1, 1.5, 1)) +
  labs(y = "") + 
  theme(
    axis.text.y = element_text(hjust = 0, angle = 0, colour = "black"),
    legend.direction = "horizontal", panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(),
    legend.position="bottom",
    legend.key.size = unit(0.4, "cm"), 
    legend.title = element_text(size = 8), 
    legend.text = element_text(size = 7),
  ) + 
  annotate("text", x = 0.55, y = 1.2, label = "F = 16.0758, p < 0.0001", size = 2.5, fontface = "bold")
main.plot_SOC_leg_peren
dev.off()

##############################################
###MBC###
##############################################
dat_MBC <- dat_sta_fin %>%
  filter(!is.na(MBC.annual.n)) %>%
  droplevels()
table(dat_MBC$soil.depth.cm)

# soil pH, crop age of perennial
dat_MBC <- dat_MBC %>% mutate(
  soil.pH.classification = case_when(
    soil.pH < 6.5 ~ "<6.5",
    soil.pH > 7.5 ~ ">7.5",
    TRUE ~ "6.5-7.5"
  ),
  crop.age.classification = case_when(
    age.peren < 5 ~ "Age<5",
    age.peren >= 5 & age.peren <= 10 ~ "Age 5-10",
    TRUE ~ "Age>10"
  ),
  SOC.initial = case_when(
    soc.peren.mean_30cm < 50 ~ "SOC<50",
    soc.peren.mean_30cm >= 50 & soc.peren.mean_30cm <= 100 ~ "SOC 50-100",
    TRUE ~ "SOC>100"
  ),
  climatic.zone_1 = sub(",.*$","", climatic.zone),
  textural.class = case_when(
    textural.class == "Silt loam" ~ "Silty loam",
    TRUE ~ textural.class
  )
  
)


####################################################
#add missing SD
####################################################

impute_missingness(dat_MBC)

dat_MBC <- dat_MBC %>%
  mutate(
    MBC.annual.SD_imputed = if_else(is.na(MBC.annual.SD), "yes", "no"),
    effect_size_id = row_number(),
    site_lat.log.elevation = interaction(lat.dec, long.dec, elevation)
  ) %>%
  group_by(site_lat.log.elevation) %>%
  mutate(site_id = cur_group_id()) %>%
  ungroup() %>%
  group_by(author) %>%
  mutate(study_id = cur_group_id()) %>%
  ungroup()


set.seed(1000)

dat_MBC_imp <- impute_SD(as.data.frame(dat_MBC),
                                 columnSDnames = c("MBC.annual.SD", "MBC.peren.SD"), columnXnames = c("MBC.annual.actual", "MBC.peren.actual"),
                                 method = "HotDeck",
                                 range = 2, M = 1)
write.csv(dat_MBC_imp, file = "dat_SOC_MBC.csv", row.names = FALSE)

##########################################################################
############################################################################

dat_MBC_imp_rom <- escalc(
  measure = "ROM", n1i = MBC.peren.n, n2i = MBC.annual.n, m1i = MBC.peren.actual, m2i = MBC.annual.actual,
  sd1i = MBC.peren.SD, sd2i = MBC.annual.SD, data = dat_MBC_imp, var.names = c("yi", "vi")
)

# 添加标准误（SE）列
dat_MBC_imp_rom$sei <- sqrt(dat_MBC_imp_rom$vi)

# 精度为 1/SE
dat_MBC_imp_rom$precision <- 1 / dat_MBC_imp_rom$sei

# 如果已有样本量列（SOC.peren.n + SOC.annual.n），可以这样创建总样本量：
dat_MBC_imp_rom$sample_size <- dat_MBC_imp_rom$MBC.peren.n + dat_MBC_imp_rom$MBC.annual.n

# 图 (a)：效应值 vs 样本量
Sample_size_2 <- ggplot(dat_MBC_imp_rom, aes(x = sample_size, y = yi)) +
  geom_point(aes(size = precision), alpha = 0.4, color = "grey30") +
  geom_smooth(method = "lm", se = TRUE, fill = "grey", color = "black") +
  scale_size_continuous(name = "Precision (1/SE)") +
  labs(x = "Sample size (number of comparisons)", y="") +
  theme(aspect.ratio = 0.9,
        panel.background = element_rect(fill = 'white', colour = 'white'),
        plot.title = element_text(hjust = 0.5, size=20, margin = ggplot2::margin(0,0,0.5,0,'cm')),
        axis.title=element_text(colour='black', size=14),
        axis.text.x = element_text(colour = 'black', size = 12,margin = ggplot2::margin(0,0,0.3,0,'cm')),
        axis.text.y = element_text(colour = 'black', size = 12,margin = ggplot2::margin(0,0,0,0.3,'cm')),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.line = element_line(colour = 'black',size=1),
        legend.position = c(0.8,0.2))

# 图 (b)：效应值 vs 标准误
SE_2 <- ggplot(dat_MBC_imp_rom, aes(x = sei, y = yi)) +
  geom_point(aes(size = precision), alpha = 0.4, color = "grey30") +
  geom_smooth(method = "lm", se = TRUE, fill = "grey", color = "black") +
  scale_size_continuous(name = "Precision (1/SE)") +
  labs(x = "Standard error", y="") +
  theme(aspect.ratio = 0.9,
        panel.background = element_rect(fill = 'white', colour = 'white'),
        plot.title = element_text(hjust = 0.5, size=20, margin = ggplot2::margin(0,0,0.5,0,'cm')),
        axis.title=element_text(colour='black', size=14),
        axis.text.x = element_text(colour = 'black', size = 12,margin = ggplot2::margin(0,0,0.3,0,'cm')),
        axis.text.y = element_text(colour = 'black', size = 12,margin = ggplot2::margin(0,0,0,0.3,'cm')),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.line = element_line(colour = 'black',size=1),
        legend.position = c(0.8,0.2))


filename <- paste('Result_Regression/',"MBC_funnel",'.jpeg', sep='')
jpeg(filename,width = 3000,height=1500,res = 300)
MBC_funnel <- Sample_size_2 + SE_2
print(MBC_funnel)
graphics.off()

summary(lm(yi ~ sample_size, data = dat_MBC_imp_rom))
summary(lm(yi ~ sei, data = dat_MBC_imp_rom))



# 3 level
main.model_MBC_imp_rom_3l <- rma.mv(yi = yi, V = vi, random = ~ 1 | study_id / site_id / effect_size_id, data = dat_MBC_imp_rom, method = "ML", test = "t")
profile(main.model_MBC_imp_rom_3l) # no issues

# 2 level
main.model_MBC_imp_rom_2l <- rma.mv(yi = yi, V = vi, random = ~ 1 | site_id / effect_size_id, data = dat_MBC_imp_rom, method = "ML", test = "t")
profile(main.model_MBC_imp_rom_2l) # no issues

anova(main.model_MBC_imp_rom_3l,main.model_MBC_imp_rom_2l)

#### main model (REML) - final model-why use 2l model?
main.model_MBC_imp_rom <- rma.mv(yi = yi, V = vi, random = ~ 1 | site_id / effect_size_id, data = dat_MBC_imp_rom, method = "REML", test = "t")
summary(main.model_MBC_imp_rom)
main.model_MBC_imp_rom$beta # effect size

# 计算影响性指标
inf_MBC <-  metafor:::influence.rma.mv(main.model_MBC_imp_rom)

# 作图：标准化残差 vs 杠杆值（帽子值）
plot(inf_SOC, which = "inf")  # 或者 which = 1，默认也行


# checks
par(mfrow=c(1,1))
forest(main.model_MBC_imp_rom, header = TRUE, addpred = TRUE, showweights = TRUE, cex = 0.75)
profile(main.model_MBC_imp_rom)

qqnorm(residuals(main.model_MBC_imp_rom), main = "QQ plot: residuals")
qqline(residuals(main.model_MBC_imp_rom), col = "red")

# variance components
main.model_MBC_imp_rom$sigma2
# total heterogeneity
sum(main.model_MBC_imp_rom$sigma2)

# obtaining I^2 (relative measurement of heterogeneity)
i2_ml(main.model_MBC_imp_rom, method = c("matrix"))

main.model_MBC_imp_rom_I2 <- i2_ml(main.model_MBC_imp_rom)
main.model_MBC_imp_rom_I2 * 100 # (% of heterogeneity)


tiff("F:/R/Meta-analysis-crop/MBC_analysis/control_monoculture_main_model_MBC.tiff", height = 6, width = 12, units = "cm", compression = "lzw", res = 600)
main.plot_MBC <- orchard_plot(main.model_MBC_imp_rom,mod = "1",group="study.number",xlab = "Effect size (lnRR)", transfm = "none") +
  scale_colour_manual(
    values = "#54bc57",
    aesthetics = c("colour", "fill")
  ) +
  labs(y = "LnRR-Microbial Biomass Carbon") +
  theme(legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.position = c(1, -0.03),
        legend.key.size = unit(0.5, 'cm'),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank()
        ) +
  labs(subtitle = "Perennial vs MBC")
main.plot_MBC$layers[[4]] <- NULL
summary_MBC <- summary(main.model_MBC_imp_rom)
estimate <- summary_MBC$beta[1]    # 点的数值
lowerCL <- summary_MBC$ci.lb       # 置信区间下限
upperCL <- summary_MBC$ci.ub       # 置信区间上限
point_data <- data.frame(
  condition = "Intrcpt",       # x轴标签
  estimate = estimate,
  lowerCL = lowerCL,
  upperCL = upperCL          # 颜色值（和之前orchard一致）
)
main.plot_MBC <- main.plot_MBC +
  geom_point(
    data = point_data,  # 你自定义的数据
    aes(x = condition, y = estimate),
    shape = 21,
    size = 3,    
    fill = "#54bc57",  # ✅ 直接指定颜色，不映射
    color = "black"
  )
main.plot_MBC
dev.off()

filename <- paste('Result_Regression/',"SOC_MBC_mon",'.jpeg', sep='')
jpeg(filename,width = 2000,height=2000,res = 300)
main.model_MBC_Line <- ggplot(dat_MBC,aes(x=effect.size.1, y=effect.size))+
  geom_hline(yintercept=0,linetype="longdash",size=0.8,colour="grey65")+
  geom_vline(xintercept=0,linetype="longdash",size=0.8,colour="grey65")+
  scale_x_continuous(limits = c(-0.9, 1.8),breaks = seq(-0.5, 1.5, 0.5))+
  scale_y_continuous(limits = c(-0.7, 1.7),breaks = seq(-0.5, 1.5, 0.5))+
  geom_point(fill="#54bc57", colour="#145742",size=3,shape=21) +
  geom_smooth(method = "lm",span=0.4,se=TRUE, colour= "black",fill="#318b6e",alpha=0.2)+ 
  labs(x="LnRR-Microbial Biomass Carbon" , y="LnRR-Soil Organic Carbon")+
  stat_regline_equation(
    aes(label = paste(..rr.label.., sep = "~~~~")), 
    label.x = min(dat_MBC$effect.size, na.rm = TRUE),
    label.y = max(dat_MBC$effect.size.1, na.rm = TRUE),
    hjust = 0.5,             
    vjust = -1,
    size = 5,
  ) +
  theme(aspect.ratio = 0.9,
        panel.background = element_rect(fill = 'white', colour = 'white'),
        plot.title = element_text(hjust = 0.5, size=20, margin = ggplot2::margin(0,0,0.5,0,'cm')),
        axis.title=element_text(colour='black', size=14),
        axis.text.x = element_text(colour = 'black', size = 12,margin = ggplot2::margin(0,0,0.3,0,'cm')),
        axis.text.y = element_text(colour = 'black', size = 12,margin = ggplot2::margin(0,0,0,0.3,'cm')),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.line = element_line(colour = 'black',size=1))
print(main.model_MBC_Line)
graphics.off()

main.model_mon_MBC_Line <- lm(effect.size ~ effect.size.1, data = dat_MBC)
summary(main.model_mon_MBC_Line)

#percentage of change for the overall effect and CIs
100 * (exp(as.numeric(main.model_MBC_imp_rom$beta)) - 1)
100 * (exp(as.numeric(main.model_MBC_imp_rom$ci.lb)) - 1)
100 * (exp(as.numeric(main.model_MBC_imp_rom$ci.ub)) - 1)

# 敏感性分析代码

dat_MBC_imp_sen <- dat_MBC_imp %>%
  filter(MBC.annual.SD_imputed == "no")

# 进行随机效应模型
dat_MBC_imp_sen1 <- escalc(
  measure = "ROM", n1i = SOC.peren.n, n2i = SOC.annual.n, m1i = soc.peren.mean_30cm, m2i = soc.annual.mean_30cm,
  sd1i = SOC.peren.SD, sd2i = SOC.annual.SD, data = dat_MBC_imp_sen, var.names = c("yi", "vi")
)

main.model_imp_MBC_sen <- rma.mv(yi = yi, V = vi, random = ~ 1 | site_id / effect_size_id, data = dat_MBC_imp_sen1, method = "REML", test = "t")

tiff("F:/R/Meta-analysis-crop/MBC_analysis/control_MBC_main_sen_model.tiff", height =8, width = 12, units = "cm", compression = "lzw", res = 600)
main.plot_MBC_sen <- orchard_plot(main.model_imp_MBC_sen,mod = "1",group="study.number",xlab = "Effect size (lnRR)", transfm = "none") +
  scale_colour_manual(
    values = "#54bc57",
    aesthetics = c("colour", "fill")
  ) +
  labs(y = "LnRR-Microbial Biomass Carbon") +
  theme(legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.position = c(1, 0.01),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank()
  ) +
  labs(subtitle = "Perennial vs MBC")
main.plot_MBC_sen
dev.off()


##################################################################
# first categorical variables with MBC
##################################################################

########sort as different climatic zone
main.model_MBC_climatic_zone_1 <- rma.mv(yi = yi, V = vi, mod = ~ climatic.zone_1 - 1, random = ~ 1 | site_id / effect_size_id, data = dat_MBC_imp_rom, method = "REML", test = "t")
anova(main.model_MBC_climatic_zone_1)

lab <- c(
  "Climatic.zone_1Subtropics" = "Subtropics",
  "Climatic.zone_1Temperate" = "Temperate",
  "Climatic.zone_1Tropics" = "Tropics"
)


tiff("F:/R/Meta-analysis-crop/MBC_analysis/control_monoculture_main_model_moderator_climate_zone.tiff", height = 7, width = 10, units = "cm", compression = "lzw", res = 600)
main.plot_MBC_climatic_zone <- orchard_plot(main.model_MBC_climatic_zone_1, mod = "climatic.zone_1",group="study.number",xlab = "Effect size (lnRR)", transfm = "none") +
  scale_x_discrete(limits = c("Temperate", "Subtropics", "Tropics"))+
  scale_y_continuous(limits = c(-1.2, 1.8),breaks = seq(-1, 1.5, 1))+
  scale_colour_manual(values = c("#36c179", "#3ed739", "#225a27"),aesthetics = c("colour", "fill")) + 
  labs(y = "") +
  theme(
    axis.text.y = element_text(hjust = 0, angle = 0, colour = "black"),
    legend.direction = "horizontal", 
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(), 
    panel.grid.major.y = element_blank(),
    legend.position="none"
    ) + 
  annotate("text", x = 0.55, y = 1.2, label = "F = 10.4132, p < 0.0001", size = 2.5, fontface = "bold")
main.plot_MBC_climatic_zone
dev.off()

#using the function to calculate percentage of change
main.model_MBC_climatic_zone_1_perc <- perc_change(main.model_MBC_climatic_zone_1 ,"climatic.zone_1")
main.model_MBC_climatic_zone_1_perc
write.csv(main.model_MBC_climatic_zone_1_perc, file = "climatic_zone_1_perc.csv", row.names = FALSE)

################################################################
#linear model with climate type
#################################################################

#####Temperate
dat_MBC_imp_rom_Tem <- dat_MBC_imp_rom %>%
  filter(climatic.zone_1 == "Temperate") %>%
  droplevels()
table(dat_MBC_imp_rom_Tem$soil.depth.cm)

#####Subtropics
dat_MBC_imp_rom_Sub <- dat_MBC_imp_rom %>%
  filter(climatic.zone_1 == "Subtropics") %>%
  droplevels()
table(dat_MBC_imp_rom_Sub$soil.depth.cm)

#####Tropics

dat_MBC_imp_rom_Tro <- dat_MBC_imp_rom %>%
  filter(climatic.zone_1 == "Tropics") %>%
  droplevels()

dat_MBC_imp_rom_Tem$clim_zone <- "Temperate"
dat_MBC_imp_rom_Sub$clim_zone <- "Subtropics"
dat_MBC_imp_rom_Tro$clim_zone <- "Tropics"

dat_MBC_combined <- rbind(dat_MBC_imp_rom_Tem, dat_MBC_imp_rom_Sub,dat_MBC_imp_rom_Tro)


filename <- paste('Result_Regression/',"SOC_MBC_Climate",'.jpeg', sep='')
jpeg(filename,width = 2000,height=2000,res = 300)
main.model_MBC_Climate <- ggplot(dat_MBC_combined, aes(x = effect.size.1, y = effect.size, color = clim_zone, fill = clim_zone)) +
  geom_hline(yintercept = 0, linetype = "longdash", color = "grey80", size = 0.8) +
  geom_vline(xintercept = 0, linetype = "longdash", color = "grey80", size = 0.8) +
  geom_point(shape = 21, size = 3, alpha = 0.8, color = "black") +
  geom_smooth(method = "lm", se = TRUE, size = 1.2, alpha = 0.3) +
  stat_poly_eq(
    aes(label = paste("bold(", ..rr.label.., ")")),
    label.x.npc = "left", label.y.npc = "top",
    formula = y ~ x,
    parse = TRUE,
    size = 4,
    fontface = "bold", 
    show.legend = FALSE
  )+
  scale_color_manual(values = c("Temperate" = "#3e7d13", "Subtropics" = "#145742","Tropics" = "#0f5f67")) +
  scale_fill_manual(values  = c("Temperate" = "#86dc4d", "Subtropics" = "#32b671","Tropics" = "#44b9c4") ) +
  labs(
    x = "LRR of MBC",
    y = "LRR of SOC",
    title = "Regression of Different Climate Zones",
    color = "Climate Zone",
    fill = "Climate Zone"
  ) +
  theme_classic(base_size = 14) +
  theme(
    aspect.ratio = 0.9,
    panel.background = element_rect(fill = 'white', colour = 'white'),
    plot.title = element_text(hjust = 0.5, size=20, face=2,margin = ggplot2::margin(0,0,0.5,0,'cm')),
    axis.title=element_text(colour='black', size=14,face=2),
    axis.text.x = element_text(colour = 'black', size = 12,margin = ggplot2::margin(0,0,0.3,0,'cm')),
    axis.text.y = element_text(colour = 'black', size = 12,margin = ggplot2::margin(0,0,0,0.3,'cm')),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    axis.line = element_line(colour = 'black',size=1),
    legend.position = c(0.275, 0.88),
    legend.title = element_blank(),
    legend.text  = element_text(face = "bold", size = 10)
  )
print(main.model_MBC_Climate)
graphics.off()



dat_MBC_combined$ln_MBC <- log(dat_MBC_combined$MBC.actual.peren.sameunit)
dat_MBC_combined <- dat_MBC_combined %>%
  filter(!is.na(ln_MBC)) %>%
  droplevels()

#####Temperate
dat_MBC_imp_inital_Tem <- dat_MBC_combined %>%
  filter(climatic.zone_1 == "Temperate") %>%
  droplevels()
table(dat_MBC_imp_rom_Tem$soil.depth.cm)

#####Subtropics
dat_MBC_imp_inital_Sub <- dat_MBC_combined %>%
  filter(climatic.zone_1 == "Subtropics") %>%
  droplevels()
table(dat_MBC_imp_rom_Sub$soil.depth.cm)

#####Tropics
dat_MBC_imp_inital_Tro <- dat_MBC_combined %>%
  filter(climatic.zone_1 == "Tropics") %>%
  droplevels()

filename <- paste('Result_Regression/',"MBCinital_Tem",'.jpeg', sep='')
jpeg(filename,width = 2000,height=2000,res = 300)
main.model_MBCinital_Tem <- ggplot(dat_MBC_imp_inital_Tem,aes(x=ln_MBC, y=effect.size.1))+
  geom_hline(yintercept=0,linetype="longdash",size=0.8,colour="grey65")+
  geom_vline(xintercept=0,linetype="longdash",size=0.8,colour="grey65")+
  scale_x_continuous(limits = c(-5, 5),breaks = seq(-5, 5, 2.5))+
  scale_y_continuous(limits = c(-1, 1.5),breaks = seq(-1, 1.5, 0.5))+
  geom_point(fill="#86dc4d", colour="#3e7d13",size=3, shape=21) +
  geom_smooth(method = "lm", span=0.4, se=TRUE, colour= "black",fill="#86dc4d",alpha=0.2)+ 
  labs(x="ln(MBC initial)" , y="LnRR-MBC", title = "Temperate")+
  theme(aspect.ratio = 0.9,
        panel.background = element_rect(fill = 'white', colour = 'white'),
        plot.title = element_text(hjust = 0.5, size=20, face=2,margin = ggplot2::margin(0,0,0.5,0,'cm')),
        axis.title=element_text(colour='black', size=14,face=2),
        axis.text.x = element_text(colour = 'black', size = 12,margin = ggplot2::margin(0,0,0.3,0,'cm')),
        axis.text.y = element_text(colour = 'black', size = 12,margin = ggplot2::margin(0,0,0,0.3,'cm')),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.line = element_line(colour = 'black',size=1))
print(main.model_MBCinital_Tem)
graphics.off()

filename <- paste('Result_Regression/',"MBCinital_Sub",'.jpeg', sep='')
jpeg(filename,width = 2000,height=2000,res = 300)
main.model_MBCinital_Sub <- ggplot(dat_MBC_imp_inital_Sub,aes(x=ln_MBC, y=effect.size.1))+
  geom_hline(yintercept=0,linetype="longdash",size=0.8,colour="grey65")+
  geom_vline(xintercept=0,linetype="longdash",size=0.8,colour="grey65")+
  scale_x_continuous(limits = c(-4.5, 6),breaks = seq(-2.5, 5, 2.5))+
  scale_y_continuous(limits = c(-0.6, 1.5),breaks = seq(-0.5, 1.5, 0.5))+
  geom_point(fill="#32b671", colour="#145742",size=3, shape=21) +
  geom_smooth(method = "lm", span=0.4, se=TRUE, colour= "black",fill="#32b671",alpha=0.2)+ 
  labs(x="ln(MBC initial)" , y="LnRR-MBC", title = "Subtropics")+
  theme(aspect.ratio = 0.9,
        panel.background = element_rect(fill = 'white', colour = 'white'),
        plot.title = element_text(hjust = 0.5, size=20, face=2,margin = ggplot2::margin(0,0,0.5,0,'cm')),
        axis.title=element_text(colour='black', size=14,face=2),
        axis.text.x = element_text(colour = 'black', size = 12,margin = ggplot2::margin(0,0,0.3,0,'cm')),
        axis.text.y = element_text(colour = 'black', size = 12,margin = ggplot2::margin(0,0,0,0.3,'cm')),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.line = element_line(colour = 'black',size=1))
print(main.model_MBCinital_Sub)
graphics.off()

filename <- paste('Result_Regression/',"MBCinital_Tro",'.jpeg', sep='')
jpeg(filename,width = 2000,height=2000,res = 300)
main.model_MBCinital_Tro <- ggplot(dat_MBC_imp_inital_Tro,aes(x=ln_MBC, y=effect.size.1))+
  geom_hline(yintercept=0,linetype="longdash",size=0.8,colour="grey65")+
  geom_vline(xintercept=0,linetype="longdash",size=0.8,colour="grey65")+
  scale_x_continuous(limits = c(0.5, 1.3),breaks = seq(0.6, 1.2, 0.3))+
  scale_y_continuous(limits = c(0.2, 1),breaks = seq(0.3, 0.9, 0.3))+
  geom_point(fill="#44b9c4", colour="#0f5f67",size=3, shape=21) +
  geom_smooth(method = "lm", span=0.4, se=TRUE, colour= "black",fill="#44b9c4",alpha=0.2)+ 
  labs(x="ln(MBC initial)" , y="LnRR-MBC", title = "Tropics")+
  theme(aspect.ratio = 0.9,
        panel.background = element_rect(fill = 'white', colour = 'white'),
        plot.title = element_text(hjust = 0.5, size=20, face=2,margin = ggplot2::margin(0,0,0.5,0,'cm')),
        axis.title=element_text(colour='black', size=14,face=2),
        axis.text.x = element_text(colour = 'black', size = 12,margin = ggplot2::margin(0,0,0.3,0,'cm')),
        axis.text.y = element_text(colour = 'black', size = 12,margin = ggplot2::margin(0,0,0,0.3,'cm')),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.line = element_line(colour = 'black',size=1))
print(main.model_MBCinital_Tro)
graphics.off()


################################################################################
# crop age perennial####3
################################################################################


# crop.age.classification
dat_MBC_imp_rom$crop.age.classification <- as.factor(dat_MBC_imp_rom$crop.age.classification)
dat_MBC_imp_rom$crop.age.classification <- factor(dat_MBC_imp_rom$crop.age.classification, levels = c("Age>10", "Age 5-10", "Age<5"))


main.model_MBC_age_peren <- rma.mv(yi = yi, V = vi, mod = ~ crop.age.classification - 1, random = ~ 1 | site_id / effect_size_id, data = dat_MBC_imp_rom, method = "REML", test = "t")
main.model_MBC_age_peren
anova(main.model_MBC_age_peren)

lab <- c(
  "Crop.age.classification>10" = "Age>10",
  "Crop.age.classification5-10" = "Age 5-10",
  "Crop.age.classification<5" = "Age<5"
)

tiff("F:/R/Meta-analysis-crop/MBC_analysis/control_monoculture_main_model_moderator_age_perennial.tiff", height = 9, width = 10, units = "cm", compression = "lzw", res = 600)
main.plot_MBC_age_peren <- orchard_plot(main.model_MBC_age_peren, mod = "crop.age.classification",group="study.number", xlab = "Effect size (lnRR)", transfm = "none") +
  scale_colour_manual(values = c("#7e1d45","#c05579","#e6778b"), aesthetics = c("colour", "fill")) + 
  scale_y_continuous(limits = c(-1.2, 1.8),breaks = seq(-1, 1.5, 1)) +
  labs(y = "") + 
  theme(
    axis.text.y = element_text(hjust = 0, angle = 0, colour = "black"),
    legend.direction = "horizontal", panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(),
    legend.position="bottom"
  ) + 
  annotate("text", x = 0.55, y = 1.2, label = "F = 19.2635, p < 0.0001", size = 2.5, fontface = "bold")
main.plot_MBC_age_peren
dev.off()

#using the function to calculate percentage of change
main.model_MBC_age_peren_perc <- perc_change(main.model_MBC_age_peren ,"crop.age.classification")
main.model_MBC_age_peren_perc
write.csv(main.model_MBC_age_peren_perc, file = "age_peren_perc.csv", row.names = FALSE)

#####crop.age.classification
dat_MBC_imp_rom_age_10 <- dat_MBC_imp_rom %>%
  filter(crop.age.classification == "Age>10") %>%
  droplevels()

dat_MBC_imp_rom_age_7.5 <- dat_MBC_imp_rom %>%
  filter(crop.age.classification == "Age 5-10") %>%
  droplevels()

dat_MBC_imp_rom_age_5 <- dat_MBC_imp_rom %>%
  filter(crop.age.classification == "Age<5") %>%
  droplevels()


dat_MBC_imp_rom_age_5$crop_age <- "Age<5"
dat_MBC_imp_rom_age_7.5$crop_age <- "Age 5-10"
dat_MBC_imp_rom_age_10$crop_age <- "Age>10"

dat_MBC_combined_age <- rbind(dat_MBC_imp_rom_age_5, dat_MBC_imp_rom_age_7.5,dat_MBC_imp_rom_age_10)

filename <- paste('Result_Regression/',"SOC_MBC_Age",'.jpeg', sep='')
jpeg(filename,width = 2000,height=2000,res = 300)
main.model_MBC_Age <- ggplot(dat_MBC_combined_age, aes(x = effect.size.1, y = effect.size, color = crop_age, fill = crop_age)) +
  geom_hline(yintercept = 0, linetype = "longdash", color = "grey80", size = 0.8) +
  geom_vline(xintercept = 0, linetype = "longdash", color = "grey80", size = 0.8) +
  geom_point(shape = 21, size = 3, alpha = 0.8, color = "black") +
  geom_smooth(method = "lm", se = TRUE, size = 1.2, alpha = 0.3) +
  stat_poly_eq(
    aes(label = paste("bold(", ..rr.label.., ")")),
    label.x.npc = "left", label.y.npc = "top",
    formula = y ~ x,
    parse = TRUE,
    size = 4,
    fontface = "bold", 
    show.legend = FALSE
  )+
  scale_color_manual(values = c("Age<5" = "#a03916", "Age 5-10" = "#a11d49","Age>10" = "#9d1a80")) +
  scale_fill_manual(values  = c("Age<5" = "#e57c59", "Age 5-10" = "#e6778b","Age>10" = "#e25fd9") ) +
  labs(
    x = "LRR of MBC",
    y = "LRR of SOC",
    title = "Regression of Different Crop Age",
    color = "Crop Age",
    fill = "Crop Age"
  ) +
  theme_classic(base_size = 14) +
  theme(
    aspect.ratio = 0.9,
    panel.background = element_rect(fill = 'white', colour = 'white'),
    plot.title = element_text(hjust = 0.5, size=20, face=2,margin = ggplot2::margin(0,0,0.5,0,'cm')),
    axis.title=element_text(colour='black', size=14,face=2),
    axis.text.x = element_text(colour = 'black', size = 12,margin = ggplot2::margin(0,0,0.3,0,'cm')),
    axis.text.y = element_text(colour = 'black', size = 12,margin = ggplot2::margin(0,0,0,0.3,'cm')),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    axis.line = element_line(colour = 'black',size=1),
    legend.position = c(0.275, 0.88),
    legend.title = element_blank(),
    legend.text  = element_text(face = "bold", size = 10)
  )
print(main.model_MBC_Age)
graphics.off()

################################################################################
# soil pH####3
################################################################################

#######set different category color

dat_MBC_imp_rom$soil.pH.classification <- as.factor(dat_MBC_imp_rom$soil.pH.classification)
dat_MBC_imp_rom$soil.pH.classification <- factor(dat_MBC_imp_rom$soil.pH.classification, levels = c(">7.5", "6.5-7.5", "<6.5"))

main.model_MBC_soil_pH <- rma.mv(yi = yi, V = vi, mod = ~ soil.pH.classification - 1, random = ~ 1 | site_id / effect_size_id, data = dat_MBC_imp_rom, method = "REML", test = "t")
main.model_MBC_soil_pH
anova(main.model_MBC_soil_pH)

lab <- c(
  "Soil.pH.classification<6.5" = "< 6.5",
  "Soil.pH.classification>7.5" = "> 7.5",
  "Soil.pH.classification6.5-7.5" = "6.5 - 7.5"
)


tiff("F:/R/Meta-analysis-crop/MBC_analysis/control_monoculture_main_model_moderator_soil_pH.tiff", height = 7, width = 10, units = "cm", compression = "lzw", res = 600)
main.plot_MBC_soil_pH <- orchard_plot(main.model_MBC_soil_pH,mod = "soil.pH.classification",group="study.number",xlab = "Effect size (lnRR)", transfm = "none") +
  scale_y_continuous(limits = c(-1.2, 1.8),breaks = seq(-1, 1.5, 1))+
  scale_colour_manual(values = c("#69176d","#c04ab5","#dc8ce5"), aesthetics = c("colour", "fill")) + 
  labs(y = "") +
  theme(
    axis.text.y = element_text(hjust = 0, angle = 0, colour = "black"),
    legend.direction = "horizontal", panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(),
    legend.position="none"
  )+
  annotate("text", x = 0.55, y = 1.2, label = "F = 9.2662, p < 0.0001", size = 2.5, fontface = "bold")
main.plot_MBC_soil_pH
dev.off()

#using the function to calculate percentage of change
main.model_MBC_soil_pH_perc <- perc_change(main.model_MBC_soil_pH ,"soil.pH.classification")
main.model_MBC_soil_pH_perc
write.csv(main.model_MBC_soil_pH_perc, file = "soil_pH_perc.csv", row.names = FALSE)


#####Soil.pH.classification<6.5
dat_MBC_imp_rom_pH6.5 <- dat_MBC_imp_rom %>%
  filter(soil.pH.classification == "<6.5") %>%
  droplevels()
table(dat_MBC_imp_rom_pH6.5$soil.depth.cm)

filename <- paste('Result_Regression/',"SOC_MBC_pH6.5",'.jpeg', sep='')
jpeg(filename,width = 2000,height=2000,res = 300)
main.model_MBC_pH6.5_Line <- ggplot(dat_MBC_imp_rom_pH6.5,aes(x=effect.size.1, y=effect.size))+
  geom_hline(yintercept=0,linetype="longdash",size=0.8,colour="grey65")+
  geom_vline(xintercept=0,linetype="longdash",size=0.8,colour="grey65")+
  scale_x_continuous(limits = c(-0.7, 1.8),breaks = seq(-0.5, 1.5, 0.5))+
  scale_y_continuous(limits = c(-0.9, 1.8),breaks = seq(-0.5, 1.5, 0.5))+
  geom_point(fill="#c04ab5", colour="#69176d",size=3,shape=21) +
  geom_smooth(method = "lm",span=0.4,se=TRUE, colour= "black",fill="#dc8ce5",alpha=0.2)+ 
  labs(x="LRR of MBC" , y="LRR of SOC", title = "pH<6.5 Regression Analysis")+
  stat_regline_equation(
    aes(label = paste(..rr.label.., sep = "~~~~")), 
    label.x = min(dat_MBC$effect.size, na.rm = TRUE),
    label.y = max(dat_MBC$effect.size.1, na.rm = TRUE),
    hjust = 0.5,             
    vjust = -1.5,
    size = 5,
  ) +
  theme(aspect.ratio = 0.9,
        panel.background = element_rect(fill = 'white', colour = 'white'),
        plot.title = element_text(hjust = 0.5, size=20, face=2,margin = ggplot2::margin(0,0,0.5,0,'cm')),
        axis.title=element_text(colour='black', size=14,face=2),
        axis.text.x = element_text(colour = 'black', size = 12,margin = ggplot2::margin(0,0,0.3,0,'cm')),
        axis.text.y = element_text(colour = 'black', size = 12,margin = ggplot2::margin(0,0,0,0.3,'cm')),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.line = element_line(colour = 'black',size=1))
print(main.model_MBC_pH6.5_Line)
graphics.off()


#####Soil.pH.classification6.5-7.5
dat_MBC_imp_rom_pH7.0 <- dat_MBC_imp_rom %>%
  filter(soil.pH.classification == "6.5-7.5") %>%
  droplevels()
table(dat_MBC_imp_rom_pH6.5$soil.depth.cm)

filename <- paste('Result_Regression/',"SOC_MBC_pH6.5-7.5",'.jpeg', sep='')
jpeg(filename,width = 2000,height=2000,res = 300)
main.model_MBC_pH7.0_Line <- ggplot(dat_MBC_imp_rom_pH7.0,aes(x=effect.size.1, y=effect.size))+
  geom_hline(yintercept=0,linetype="longdash",size=0.8,colour="grey65")+
  geom_vline(xintercept=0,linetype="longdash",size=0.8,colour="grey65")+
  scale_x_continuous(limits = c(-1, 1.5),breaks = seq(-0.5, 1.5, 0.5))+
  scale_y_continuous(limits = c(-0.7, 1.8),breaks = seq(-0.5, 1.5, 0.5))+
  geom_point(fill="#c04ab5", colour="#69176d",size=3,shape=21) +
  geom_smooth(method = "lm",span=0.4,se=TRUE, colour= "black",fill="#dc8ce5",alpha=0.2)+ 
  labs(x="LRR of MBC" , y="LRR of SOC", title = "pH=6.5-7.5 Regression Analysis")+
  stat_regline_equation(
    aes(label = paste(..rr.label.., sep = "~~~~")), 
    label.x = min(dat_MBC$effect.size, na.rm = TRUE),
    label.y = max(dat_MBC$effect.size.1, na.rm = TRUE),
    hjust = 0.9,             
    vjust = -1.5,
    size = 5,
  ) +
  theme(aspect.ratio = 0.9,
        panel.background = element_rect(fill = 'white', colour = 'white'),
        plot.title = element_text(hjust = 0.5, size=20, face=2,margin = ggplot2::margin(0,0,0.5,0,'cm')),
        axis.title=element_text(colour='black', size=14,face=2),
        axis.text.x = element_text(colour = 'black', size = 12,margin = ggplot2::margin(0,0,0.3,0,'cm')),
        axis.text.y = element_text(colour = 'black', size = 12,margin = ggplot2::margin(0,0,0,0.3,'cm')),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.line = element_line(colour = 'black',size=1))
print(main.model_MBC_pH7.0_Line)
graphics.off()

#####Soil.pH.classification>7.5
dat_MBC_imp_rom_pH7.5 <- dat_MBC_imp_rom %>%
  filter(soil.pH.classification == ">7.5") %>%
  droplevels()
table(dat_MBC_imp_rom_pH6.5$soil.depth.cm)

filename <- paste('Result_Regression/',"SOC_MBC_pH7.5",'.jpeg', sep='')
jpeg(filename,width = 2000,height=2000,res = 300)
main.model_MBC_pH7.5_Line <- ggplot(dat_MBC_imp_rom_pH7.5,aes(x=effect.size.1, y=effect.size))+
  geom_hline(yintercept=0,linetype="longdash",size=0.8,colour="grey65")+
  geom_vline(xintercept=0,linetype="longdash",size=0.8,colour="grey65")+
  scale_x_continuous(limits = c(-0.5, 1.5),breaks = seq(-0.5, 1.5, 0.5))+
  scale_y_continuous(limits = c(-0.7, 1.5),breaks = seq(-0.5, 1.5, 0.5))+
  geom_point(fill="#c04ab5", colour="#69176d",size=3,shape=21) +
  geom_smooth(method = "lm",span=0.4,se=TRUE, colour= "black",fill="#dc8ce5",alpha=0.2)+ 
  labs(x="LRR of MBC" , y="LRR of SOC", title = "pH>7.5 Regression Analysis")+
  stat_regline_equation(
    aes(label = paste(..rr.label.., sep = "~~~~")), 
    label.x = min(dat_MBC$effect.size, na.rm = TRUE),
    label.y = max(dat_MBC$effect.size.1, na.rm = TRUE),
    hjust = -0.1,             
    vjust = 1,
    size = 5,
  ) +
  theme(aspect.ratio = 0.9,
        panel.background = element_rect(fill = 'white', colour = 'white'),
        plot.title = element_text(hjust = 0.5, size=20, face=2,margin = ggplot2::margin(0,0,0.5,0,'cm')),
        axis.title=element_text(colour='black', size=14,face=2),
        axis.text.x = element_text(colour = 'black', size = 12,margin = ggplot2::margin(0,0,0.3,0,'cm')),
        axis.text.y = element_text(colour = 'black', size = 12,margin = ggplot2::margin(0,0,0,0.3,'cm')),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.line = element_line(colour = 'black',size=1))
print(main.model_MBC_pH7.5_Line)
graphics.off()

#################################################################################
# PFT_peren###3
#################################################################################

dat_MBC_imp_rom$PFT.peren <- as.factor(dat_MBC_imp_rom$PFT.peren)
dat_MBC_imp_rom$PFT.peren <- factor(dat_MBC_imp_rom$PFT.peren, levels = c("C4-C3", "C4", "C3"))

main.model_MBC_PFT_peren <- rma.mv(yi = yi, V = vi, mod = ~ PFT.peren - 1, random = ~ 1 | site_id / effect_size_id, data = dat_MBC_imp_rom, method = "REML", test = "t")
main.model_MBC_PFT_peren
anova(main.model_MBC_PFT_peren)

lab <- c(
  "PFT.perenC3" = "C3",
  "PFT.perenC4" = "C4",
  "PFT.perenC4-C3" = "C3 - C4"
)


tiff("F:/R/Meta-analysis-crop/MBC_analysis/control_monoculture_main_model_moderator_PFT_perennial.tiff", height = 7, width = 10, units = "cm", compression = "lzw", res = 600)
main.plot_MBC_PFT_peren <- orchard_plot(main.model_MBC_PFT_peren, mod = "PFT.peren",group="study.number",xlab = "Effect size (lnRR)", transfm = "none") +
  scale_colour_manual(values = c("#1d3683","#3b67b3","#62a5d8"), aesthetics = c("colour", "fill")) + 
  scale_y_continuous(limits = c(-1.2, 1.8),breaks = seq(-1, 1.5, 1)) +
  labs(y = "") + 
  theme(
    axis.text.y = element_text(hjust = 0, angle = 0, colour = "black"),
    legend.direction = "horizontal", panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(),
    legend.position="none"
    ) + 
  annotate("text", x = 0.55, y = 1.2, label = "F = 9.4817, p < 0.0001", size = 2.5, fontface = "bold")
main.plot_MBC_PFT_peren
dev.off()

#using the function to calculate percentage of change
main.model_MBC_PFT_peren_perc <- perc_change(main.model_MBC_PFT_peren ,"PFT.peren")
main.model_MBC_PFT_peren_perc
write.csv(main.model_MBC_PFT_peren_perc, file = "PFT_peren_perc.csv", row.names = FALSE)


#####PFT_perenC3
dat_MBC_imp_rom_C3 <- dat_MBC_imp_rom %>%
  filter(PFT.peren == "C3") %>%
  droplevels()
table(dat_MBC_imp_rom_pH6.5$soil.depth.cm)

filename <- paste('Result_Regression/',"SOC_MBC_C3",'.jpeg', sep='')
jpeg(filename,width = 2000,height=2000,res = 300)
main.model_MBC_C3_Line <- ggplot(dat_MBC_imp_rom_C3,aes(x=effect.size.1, y=effect.size))+
  geom_hline(yintercept=0,linetype="longdash",size=0.8,colour="grey65")+
  geom_vline(xintercept=0,linetype="longdash",size=0.8,colour="grey65")+
  scale_x_continuous(limits = c(-0.9, 1.6),breaks = seq(-0.5, 1.5, 0.5))+
  scale_y_continuous(limits = c(-0.9, 1.7),breaks = seq(-0.5, 1.5, 0.5))+
  geom_point(fill="#3b67b3", colour="#1d3683",size=3,shape=21) +
  geom_smooth(method = "lm",span=0.4,se=TRUE, colour= "black",fill="#62a5d8",alpha=0.2)+ 
  labs(x="LRR of MBC" , y="LRR of SOC", title = "PFT C3 Regression Analysis")+
  stat_regline_equation(
    aes(label = paste(..rr.label.., sep = "~~~~")), 
    label.x = min(dat_MBC$effect.size, na.rm = TRUE),
    label.y = max(dat_MBC$effect.size.1, na.rm = TRUE),
    hjust = 0.5,             
    vjust = -0.5,
    size = 5,
  ) +
  theme(aspect.ratio = 0.9,
        panel.background = element_rect(fill = 'white', colour = 'white'),
        plot.title = element_text(hjust = 0.5, size=20, face=2,margin = ggplot2::margin(0,0,0.5,0,'cm')),
        axis.title=element_text(colour='black', size=14,face=2),
        axis.text.x = element_text(colour = 'black', size = 12,margin = ggplot2::margin(0,0,0.3,0,'cm')),
        axis.text.y = element_text(colour = 'black', size = 12,margin = ggplot2::margin(0,0,0,0.3,'cm')),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.line = element_line(colour = 'black',size=1))
print(main.model_MBC_C3_Line)
graphics.off()


#####PFT_perenC3-C4
dat_MBC_imp_rom_C3_C4 <- dat_MBC_imp_rom %>%
  filter(PFT.peren == "C4-C3") %>%
  droplevels()
table(dat_MBC_imp_rom_pH6.5$soil.depth.cm)

filename <- paste('Result_Regression/',"SOC_MBC_C3-C4",'.jpeg', sep='')
jpeg(filename,width = 2000,height=2000,res = 300)
main.model_MBC_C3_C4_Line <- ggplot(dat_MBC_imp_rom_C3_C4,aes(x=effect.size.1, y=effect.size))+
  geom_hline(yintercept=0,linetype="longdash",size=0.8,colour="grey65")+
  geom_vline(xintercept=0,linetype="longdash",size=0.8,colour="grey65")+
  scale_x_continuous(limits = c(-0.7, 1.3),breaks = seq(-0.5, 1.0, 0.5))+
  scale_y_continuous(limits = c(-0.9, 1.5),breaks = seq(-0.5, 1.5, 0.5))+
  geom_point(fill="#3b67b3", colour="#1d3683",size=3,shape=21) +
  geom_smooth(method = "lm",span=0.4,se=TRUE, colour= "black",fill="#62a5d8",alpha=0.2)+ 
  labs(x="LRR of MBC" , y="LRR of SOC", title = "PFT C3-C4 Regression Analysis")+
  stat_regline_equation(
    aes(label = paste(..rr.label.., sep = "~~~~")), 
    label.x = min(dat_MBC$effect.size, na.rm = TRUE),
    label.y = max(dat_MBC$effect.size.1, na.rm = TRUE),
    hjust = 0.5,             
    vjust = 0.8,
    size = 5,
  ) +
  theme(aspect.ratio = 0.9,
        panel.background = element_rect(fill = 'white', colour = 'white'),
        plot.title = element_text(hjust = 0.5, size=20, face=2,margin = ggplot2::margin(0,0,0.5,0,'cm')),
        axis.title=element_text(colour='black', size=14,face=2),
        axis.text.x = element_text(colour = 'black', size = 12,margin = ggplot2::margin(0,0,0.3,0,'cm')),
        axis.text.y = element_text(colour = 'black', size = 12,margin = ggplot2::margin(0,0,0,0.3,'cm')),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.line = element_line(colour = 'black',size=1))
print(main.model_MBC_C3_C4_Line)
graphics.off()

#####PFT_perenC3-C4
dat_MBC_imp_rom_C4 <- dat_MBC_imp_rom %>%
  filter(PFT.peren == "C4") %>%
  droplevels()
table(dat_MBC_imp_rom_pH6.5$soil.depth.cm)

filename <- paste('Result_Regression/',"SOC_MBC_C4",'.jpeg', sep='')
jpeg(filename,width = 2000,height=2000,res = 300)
main.model_MBC_C4_Line <- ggplot(dat_MBC_imp_rom_C4,aes(x=effect.size.1, y=effect.size))+
  geom_hline(yintercept=0,linetype="longdash",size=0.8,colour="grey65")+
  geom_vline(xintercept=0,linetype="longdash",size=0.8,colour="grey65")+
  scale_x_continuous(limits = c(-0.7, 1.8),breaks = seq(-0.5, 1.5, 0.5))+
  scale_y_continuous(limits = c(-0.5, 1.5),breaks = seq(-0.5, 1.5, 0.5))+
  geom_point(fill="#3b67b3", colour="#1d3683",size=3,shape=21) +
  geom_smooth(method = "lm",span=0.4,se=TRUE, colour= "black",fill="#62a5d8",alpha=0.2)+ 
  labs(x="LRR of MBC" , y="LRR of SOC", title = "PFT C4 Regression Analysis")+
  stat_regline_equation(
    aes(label = paste(..rr.label.., sep = "~~~~")), 
    label.x = min(dat_MBC$effect.size, na.rm = TRUE),
    label.y = max(dat_MBC$effect.size.1, na.rm = TRUE),
    hjust = 0.3,             
    vjust = 0.7,
    size = 5,
  ) +
  theme(aspect.ratio = 0.9,
        panel.background = element_rect(fill = 'white', colour = 'white'),
        plot.title = element_text(hjust = 0.5, size=20, face=2,margin = ggplot2::margin(0,0,0.5,0,'cm')),
        axis.title=element_text(colour='black', size=14,face=2),
        axis.text.x = element_text(colour = 'black', size = 12,margin = ggplot2::margin(0,0,0.3,0,'cm')),
        axis.text.y = element_text(colour = 'black', size = 12,margin = ggplot2::margin(0,0,0,0.3,'cm')),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.line = element_line(colour = 'black',size=1))
print(main.model_MBC_C4_Line)
graphics.off()




################################################################################
# leg.peren#####2
################################################################################

main.model_MBC_leg_peren <- rma.mv(yi = yi, V = vi, mod = ~ leg.peren - 1, random = ~ 1 | site_id / effect_size_id, data = dat_MBC_imp_rom, method = "REML", test = "t")
main.model_MBC_leg_peren
anova(main.model_MBC_leg_peren)
lab <- c(
  "No" = "Non-legume",
  "Yes" = "Legume "
)

tiff("F:/R/Meta-analysis-crop/MBC_analysis/control_monoculture_main_model_moderator_legume_perennial.tiff", height = 7.5, width = 11, units = "cm", compression = "lzw", res = 600)
main.plot_MBC_leg_peren <- orchard_plot(main.model_MBC_leg_peren, mod = "leg.peren",group="study.number", xlab = "Effect size (lnRR)", transfm = "none") +
  scale_x_discrete(labels = c("Yes"="Legume","No"="Non-legume"))+
  scale_colour_manual(values = c("#e2a06b","#c16528"), aesthetics = c("colour", "fill")) + 
  scale_y_continuous(limits = c(-1.2, 1.8),breaks = seq(-1, 1.5, 1)) +
  labs(y = "") + 
  theme(
    axis.text.y = element_text(hjust = 0, angle = 0, colour = "black"),
    legend.direction = "horizontal", panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(),
    legend.position="bottom",
    legend.key.size = unit(0.4, "cm"), 
    legend.title = element_text(size = 8), 
    legend.text = element_text(size = 7),
  ) + 
  annotate("text", x = 0.55, y = 1.2, label = "F = 14.3354, p < 0.0001", size = 2.5, fontface = "bold")
main.plot_MBC_leg_peren
dev.off()

#using the function to calculate percentage of change
main.model_MBC_leg_peren_perc <- perc_change(main.model_MBC_leg_peren ,"leg.peren")
main.model_MBC_leg_peren_perc


dat_MBC_imp_rom_leg <- dat_MBC_imp_rom %>%
  filter(leg.peren == "Yes") %>%
  droplevels()
table(dat_MBC_imp_rom_pH6.5$soil.depth.cm)

filename <- paste('Result_Regression/',"SOC_MBC_leg",'.jpeg', sep='')
jpeg(filename,width = 2000,height=2000,res = 300)
main.model_MBC_leg_Line <- ggplot(dat_MBC_imp_rom_leg,aes(x=effect.size.1, y=effect.size))+
  geom_hline(yintercept=0,linetype="longdash",size=0.8,colour="grey65")+
  geom_vline(xintercept=0,linetype="longdash",size=0.8,colour="grey65")+
  scale_x_continuous(limits = c(-0.6, 1.5),breaks = seq(-0.5, 1.5, 0.5))+
  scale_y_continuous(limits = c(-0.6, 1.5),breaks = seq(-0.5, 1.5, 0.5))+
  geom_point(fill="#e2a06b", colour="#c16528",size=3,shape=21) +
  geom_smooth(method = "lm",span=0.4,se=TRUE, colour= "black",fill="#e0c37c",alpha=0.2)+ 
  labs(x="LRR of MBC" , y="LRR of SOC", title = "Legume Regression Analysis")+
  stat_regline_equation(
    aes(label = paste(..rr.label.., sep = "~~~~")), 
    label.x = min(dat_MBC$effect.size, na.rm = TRUE),
    label.y = max(dat_MBC$effect.size.1, na.rm = TRUE),
    hjust = 0.1,              
    vjust = 2,              
    size = 5
  ) +
  theme(aspect.ratio = 0.9,
        panel.background = element_rect(fill = 'white', colour = 'white'),
        plot.title = element_text(hjust = 0.5, size=20, face=2,margin = ggplot2::margin(0,0,0.5,0,'cm')),
        axis.title=element_text(colour='black', size=14,face=2),
        axis.text.x = element_text(colour = 'black', size = 12,margin = ggplot2::margin(0,0,0.3,0,'cm')),
        axis.text.y = element_text(colour = 'black', size = 12,margin = ggplot2::margin(0,0,0,0.3,'cm')),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.line = element_line(colour = 'black',size=1))
print(main.model_MBC_leg_Line)
graphics.off()

dat_MBC_imp_rom_noleg <- dat_MBC_imp_rom %>%
  filter(leg.peren == "No") %>%
  droplevels()
table(dat_MBC_imp_rom_pH6.5$soil.depth.cm)

filename <- paste('Result_Regression/',"SOC_MBC_noleg",'.jpeg', sep='')
jpeg(filename,width = 2000,height=2000,res = 300)
main.model_MBC_noleg_Line <- ggplot(dat_MBC_imp_rom_noleg,aes(x=effect.size.1, y=effect.size))+
  geom_hline(yintercept=0,linetype="longdash",size=0.8,colour="grey65")+
  geom_vline(xintercept=0,linetype="longdash",size=0.8,colour="grey65")+
  scale_x_continuous(limits = c(-0.9, 1.7),breaks = seq(-0.5, 1.5, 0.5))+
  scale_y_continuous(limits = c(-0.9, 1.7),breaks = seq(-0.5, 1.5, 0.5))+
  geom_point(fill="#e2a06b", colour="#c16528",size=3,shape=21) +
  geom_smooth(method = "lm",span=0.4,se=TRUE, colour= "black",fill="#e0c37c",alpha=0.2)+ 
  labs(x="LRR of MBC" , y="LRR of SOC", title = "Non-legume Regression Analysis")+
  stat_regline_equation(
    aes(label = paste(..rr.label.., sep = "~~~~")), 
    label.x = min(dat_MBC$effect.size, na.rm = TRUE),
    label.y = max(dat_MBC$effect.size.1, na.rm = TRUE),
    hjust = 0.1,              
    vjust = 2,              
    size = 5
  ) +
  theme(aspect.ratio = 0.9,
        panel.background = element_rect(fill = 'white', colour = 'white'),
        plot.title = element_text(hjust = 0.5, size=20, face=2,margin = ggplot2::margin(0,0,0.5,0,'cm')),
        axis.title=element_text(colour='black', size=14,face=2),
        axis.text.x = element_text(colour = 'black', size = 12,margin = ggplot2::margin(0,0,0.3,0,'cm')),
        axis.text.y = element_text(colour = 'black', size = 12,margin = ggplot2::margin(0,0,0,0.3,'cm')),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.line = element_line(colour = 'black',size=1))
print(main.model_MBC_noleg_Line)
graphics.off()

#################################################################################
# location plot
#################################################################################

# loading dataset from
load(file = "F:/R/Meta-analysis-crop/metaanalysis_soc_perennial_vs_rotation.RData")

# plot

dat <- dat_sta_fin_SOC
dat$crop.age.classification <- as.factor(dat$crop.age.classification)
dat$crop.age.classification <- factor(dat$crop.age.classification, levels = levels(dat$crop.age.classification)[c(3,2,1)])

dat$MBC_presence <- ifelse(is.na(dat$effect.size.1), "No MBC", "Include MBC")
dat$MBC_presence <- factor(dat$MBC_presence, levels=c("Include MBC","No MBC"))


world_map <- map_data("world")

Global.Fig <- ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), colour = "grey70", fill = "grey70") +
  theme(
    axis.line = element_blank(), panel.background = element_rect(fill = "white", colour = "white"), axis.ticks = element_blank(),
    panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()
  ) +
  # scale_x_continuous(limits = c(-179,180))+scale_y_continuous(limits = c(-84,82))+
  ylab("Latitude") +
  xlab("Longitude") +
  theme_bw() +
  theme(legend.title = element_text(size = 9), legend.text = element_text(size = 8), legend.position = c(0.1, 0.5)) +
  theme(axis.title = element_text(size = 12), 
        panel.grid.minor.x = element_blank(), 
        panel.grid.minor.y = element_blank()
        ) +
  geom_point(data = dat, aes(x = long.dec, y = lat.dec, color = MBC_presence, size = age.peren), alpha = 0.7) +
  scale_colour_manual(values = c("Include MBC" = "#54bc57", "No MBC" = "#a55213"))+
  guides(color = guide_legend(title = "Comparison of perennial systems vs.", order = 1), shape = guide_legend(title = "Plant functional type", order = 2), size = guide_legend(title = "Perennial crop age", order = 3))

Global.Fig


# 假设你的数据集为 dat
dat_long <- dat %>%
  rename(LRR.SOC = effect.size,  # 自定义新名称
         LRR.MBC = effect.size.1) %>%
  pivot_longer(cols = c(LRR.SOC, LRR.MBC), 
               names_to = "Comparison_Type",  # 新的分类列
               values_to = "Effect_Size")    # 合并的新列


# 查看新数据集
head(dat_long)



dodge <- position_dodge(width = 0.7)
Dis1 <- ggplot(dat_long, aes(x = Effect_Size, color = Comparison_Type, fill = Comparison_Type)) +
  geom_histogram(aes(y = ..density..), alpha = 0.3, position = "dodge") +
  geom_density(alpha = 0.1, size = 1) +
  theme_bw() +
  theme(axis.title = element_text(size = 12)) +
  theme(
    legend.position = c(0.25, 0.85), legend.title = element_text(size = 5), legend.text = element_text(size = 4),
    legend.background = element_blank(), legend.key.size = unit(0.3, "cm"), 
    panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()
  ) +
  scale_colour_manual(values = c("LRR.MBC" = "#54bc57", "LRR.SOC" = "#a55213")) +
  scale_fill_manual(values = c("LRR.MBC" = "#54bc57", "LRR.SOC" = "#a55213")) +
  ylab("Density") +
  xlab("LnRR Distribution") +
  labs(color = "Comparison of perennial systems vs.", fill = "Comparison of perennial systems vs.")
Dis1



Global.Fig.1 <- ggdraw(Global.Fig) + 
                draw_plot(Dis1, x = 0.05, y = 0.07, width = 0.25, height = 0.32) + 
                draw_plot_label("(A)", size=12)
Dis1.1 <- ggdraw(main.plot_SOC) + draw_plot_label("(B)", size=12)
Dis2.1 <- ggdraw(main.plot_MBC) + draw_plot_label("(C)", size=12)
Dis3.1 <- ggdraw(main.model_MBC_Line) + draw_plot_label("(D)", size=12)

# 右侧三个图竖直排列，调整大小：
right_panel <- plot_grid(Dis1.1, Dis2.1, Dis3.1, 
                         ncol = 1, 
                         rel_heights = c(0.8, 0.8, 2.3)) # 第三个图较大

# 左边全局图，右边三个图：
Glo.Dis.F <- plot_grid(Global.Fig.1, right_panel, 
                       ncol=2, 
                       rel_widths = c(4.9,2))

ggsave("F:/R/Meta-analysis-crop/general_figure_perennial_MBC.tiff", dpi = 100, compression = "lzw", plot = Glo.Dis.F, width = 18.3, height = 7.95)


#################################################################################
###################################
# Outlier/Influence Diagnostics

# # standardized residuals
# plot(resid(main.model_mon))
# # Cook's distances: reflects how the removal of a study influences the predicted effect of all studies
# plot(cooks.distance(main.model_mon)) # based on the effect sizes [A Cook’s Distance is often considered large if Di > 4/n (heuristic) or Di  > 0.45]
# plot(cooks.distance(main.model_mon, cluster = dat_monoculture_imp_rom$site_id)) # based on the sites
# plot(cooks.distance(main.model_mon, cluster = dat_monoculture_imp_rom$study_id)) # based on the studies
# # hatvalues: indication which points have a strong influence on the results
# plot(hatvalues(main.model_mon))
# plot(hatvalues(main.model_mon, cluster = dat_monoculture_imp_rom$study_id))
# # DFBETAS:reflects how the removal of a study influences the estimated value of a particular model coefficient
# plot(dfbetas(main.model_mon)) # based on the effect sizes
# plot(dfbetas(main.model_mon, cluster = dat_monoculture_imp_rom$study_id)) # based on the studies
#
# plot(hatvalues(main.model_mon), cooks.distance(main.model_mon))
# plot(hatvalues(main.model_mon), resid(main.model_mon))

influence_dat <- cbind(dat_monoculture_imp_rom, hatvalues(main.model_mon), cooks.distance(main.model_mon), resid(main.model_mon))
influence_dat <- influence_dat %>% rename(
  hat_value = "hatvalues(main.model_mon)",
  cook_dist = "cooks.distance(main.model_mon)",
  resid = "resid(main.model_mon)"
)

tiff("F:/R/Meta-analysis-crop/control_monoculture_main_model_influence_hatvalues_residuals.tiff", height = 11, width = 11, units = "cm", compression = "lzw", res = 600)
ggplot(influence_dat, aes(x = hat_value, y = resid)) +
  geom_point(aes(size = sqrt(1 / (vi))), shape = 21, color = "black", fill = "grey", alpha = 0.5) +
  labs(x = "Leverage (hat values)", y = "Standardized residuals") +
  guides(size = guide_legend("Precision (1/SE)")) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = c(0.975, 0.925), legend.justification = c(1, 1), legend.direction = "horizontal",
    legend.title = element_text(size = 10.5), legend.text = element_text(size = 10.5),
    panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(),
    axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black")
  ) +
  geom_vline(xintercept = 2 * mean(influence_dat$hat_value), linetype = "dashed", color = "grey50", size = 0.75)+
  labs(tag = "A",
       subtitle = "Perennial vs monoculture")
dev.off()

