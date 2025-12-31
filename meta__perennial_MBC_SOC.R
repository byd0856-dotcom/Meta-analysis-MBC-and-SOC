getwd()
setwd("F:/R/Meta-analysis-crop/MBC_analysis")

library(readxl)
library(readr)
library(dplyr)
library(tidyr)
library(broom) 
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
    soil.pH < 6.5 ~ "pH<6.5",
    soil.pH > 7.5 ~ "pH>7.5",
    TRUE ~ "pH 6.5-7.5"
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
  ) 
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

tiff("F:/R/Meta-analysis-crop/MBC_analysis/control_monoculture_main_model_moderator_climate_zone_SOC.tiff", height = 8, width = 13, units = "cm", compression = "lzw", res = 600)
main.plot_SOC_climatic_zone <- orchard_plot(main.model_SOC_climatic_zone_1, mod = "climatic.zone_1",group="study.number",xlab = "Effect size (lnRR)", transfm = "none",g=FALSE, branch.size=0.8) +
  scale_x_discrete(limits = c("Temperate", "Subtropics", "Tropics"))+
  scale_y_continuous(limits = c(-1.6, 1.8),breaks = seq(-1, 1.5, 1))+
  scale_colour_manual(values = c("#36c179", "#3ed739", "#225a27"),aesthetics = c("colour", "fill")) + 
  labs(y = "") +  
  theme(
    axis.text.y = element_text(hjust = 0, angle = 0, size = 12, colour = "black"),
    legend.direction = "horizontal", panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(), 
    panel.grid.major.y = element_blank(),
    legend.position="none"
  ) + 
  annotate("text", x = 0.55, y = 1.2, label = "F = 11.758, p < 0.001", size = 2.5, fontface = "bold")
main.plot_SOC_climatic_zone$layers[[4]] <- NULL
summary_SOC_climatic <- summary(main.model_SOC_climatic_zone_1)
estimate <- summary_SOC_climatic$beta    # 点的数值
lowerCL <- summary_SOC_climatic$ci.lb       # 置信区间下限
upperCL <- summary_SOC_climatic$ci.ub       # 置信区间上限
point_data_climatic <- data.frame(
  condition = c( "Subtropics","Temperate", "Tropics"),       # x轴标签
  estimate = estimate,
  lowerCL = lowerCL,
  upperCL = upperCL          # 颜色值（和之前orchard一致）
)
main.plot_SOC_climatic_zone <- main.plot_SOC_climatic_zone +
  geom_point(
    data = point_data_climatic,  # 你自定义的数据
    aes(x = condition, y = estimate),
    shape = 21,
    size = 2.4,    
    fill = c("#36c179", "#3ed739", "#225a27"),  # ✅ 直接指定颜色，不映射
    color = "black"
  )
main.plot_SOC_climatic_zone
dev.off()

#using the function to calculate percentage of change
main.model_SOC_climatic_zone_1_perc <- perc_change(main.model_SOC_climatic_zone_1 ,"climatic.zone_1")
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

tiff("F:/R/Meta-analysis-crop/MBC_analysis/control_monoculture_main_model_moderator_age_perennial_SOC.tiff", height = 9.5, width = 13, units = "cm", compression = "lzw", res = 600)
main.plot_SOC_age_peren <- orchard_plot(main.model_SOC_age_peren, mod = "crop.age.classification",group="study.number", xlab = "Effect size (lnRR)", transfm = "none", g=FALSE, branch.size=0.8) +
  scale_colour_manual(values = c("#7e1d45","#c05579","#e6778b"), aesthetics = c("colour", "fill")) + 
  scale_y_continuous(limits = c(-1.6, 1.8),breaks = seq(-1, 1.5, 1)) +
  labs(y = "") + 
  theme(
    axis.text.y = element_text(hjust = 0, angle = 0, size = 12, colour = "black"),
    legend.direction = "horizontal", panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(),
    legend.position="bottom"
  ) + 
  annotate("text", x = 0.55, y = 1.2, label = "F = 30.688, p < 0.001", size = 2.5, fontface = "bold")
main.plot_SOC_age_peren$layers[[4]] <- NULL
summary_SOC_age <- summary(main.model_SOC_age_peren)
estimate <- summary_SOC_age$beta    # 点的数值
lowerCL <- summary_SOC_age$ci.lb       # 置信区间下限
upperCL <- summary_SOC_age$ci.ub       # 置信区间上限
point_data_age <- data.frame(
  condition = c("Age>10", "Age 5-10", "Age<5"),       # x轴标签
  estimate = estimate,
  lowerCL = lowerCL,
  upperCL = upperCL          # 颜色值（和之前orchard一致）
)
main.plot_SOC_age_peren <- main.plot_SOC_age_peren +
  geom_point(
    data = point_data_age,  # 你自定义的数据
    aes(x = condition, y = estimate),
    shape = 21,
    size = 2.4,    
    fill = c("#7e1d45","#c05579","#e6778b"),  # ✅ 直接指定颜色，不映射
    color = "black"
  )
main.plot_SOC_age_peren
dev.off()

#using the function to calculate percentage of change
main.model_SOC_crop_age_perc <- perc_change(main.model_SOC_age_peren ,"crop.age.classification")
main.model_SOC_crop_age_perc

# pH area
######################################
dat_sta_fin_SOC$soil.pH.classification <- as.factor(dat_sta_fin_SOC$soil.pH.classification)
dat_sta_fin_SOC$soil.pH.classification <- factor(dat_sta_fin_SOC$soil.pH.classification, levels = c("pH>7.5", "pH 6.5-7.5", "pH<6.5"))

main.model_SOC_soil_pH <- rma.mv(yi = yi, V = vi, mod = ~ soil.pH.classification - 1, random = ~ 1 | site_id / effect_size_id, data = dat_sta_fin_SOC, method = "REML", test = "t")
main.model_SOC_soil_pH
anova(main.model_SOC_soil_pH)

lab <- c(
  "Soil.pH.classification<6.5" = "< 6.5",
  "Soil.pH.classification>7.5" = "> 7.5",
  "Soil.pH.classification6.5-7.5" = "6.5 - 7.5"
)

tiff("F:/R/Meta-analysis-crop/MBC_analysis/control_monoculture_main_model_moderator_soil_pH_SOC.tiff", height = 8, width = 13, units = "cm", compression = "lzw", res = 600)
main.plot_SOC_soil_pH <- orchard_plot(main.model_SOC_soil_pH,mod = "soil.pH.classification",group="study.number",xlab = "Effect size (lnRR)", transfm = "none", g=FALSE, branch.size=0.8) +
  scale_x_discrete(labels = c("pH>7.5", "pH 6.5–7.5", "pH<6.5"))+
  scale_y_continuous(limits = c(-1.6, 1.8),breaks = seq(-1, 1.5, 1))+
  scale_colour_manual(values = c("#125f86","#4299cb","#7cc7e4"), aesthetics = c("colour", "fill")) + 
  labs(y = "") +
  theme(
    axis.text.y = element_text(hjust = 0, angle = 0, size = 12, colour = "black"),
    legend.direction = "horizontal", panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(),
    legend.position="none"
  )+
  annotate("text", x = 0.55, y = 1.2, label = "F = 11.6453, p < 0.001", size = 2.5, fontface = "bold")
main.plot_SOC_soil_pH$layers[[4]] <- NULL
summary_SOC_pH <- summary(main.model_SOC_soil_pH)
estimate <- summary_SOC_pH$beta    # 点的数值
lowerCL <- summary_SOC_pH$ci.lb       # 置信区间下限
upperCL <- summary_SOC_pH$ci.ub       # 置信区间上限
point_data_age <- data.frame(
  condition = c("PH>7.5", "PH 6.5-7.5", "PH<6.5"),       # x轴标签
  estimate = estimate,
  lowerCL = lowerCL,
  upperCL = upperCL          # 颜色值（和之前orchard一致）
)
main.plot_SOC_soil_pH <- main.plot_SOC_soil_pH +
  geom_point(
    data = point_data_age,  # 你自定义的数据
    aes(x = condition, y = estimate),
    shape = 21,
    size = 2.4,    
    fill = c("#125f86","#4299cb","#7cc7e4"),  # ✅ 直接指定颜色，不映射
    color = "black")
main.plot_SOC_soil_pH
dev.off()

#using the function to calculate percentage of change
main.model_SOC_soil_pH_perc <- perc_change(main.model_SOC_soil_pH ,"soil.pH.classification")
main.model_SOC_soil_pH_perc


# vegetation type
main.model_vegetation_type <- rma.mv(yi = yi, V = vi, mod = ~ vegetation.type - 1, random = ~ 1 | site_id / effect_size_id, data = dat_sta_fin_SOC, method = "REML", test = "t")
main.model_vegetation_type

lab <- c(
  "Vegetation.typeHerbaceous" = "Herbaceous",
  "Vegetation.typeWoody" = "Woody"
)

tiff("F:/R/Meta-analysis-crop/control_monoculture_main_model_moderator_vegetation_type.tiff", height = 11, width = 11, units = "cm", compression = "lzw", res = 600)
main.plot_mon_vegetation_type <- orchard_plot(main.model_vegetation_type,mod = "vegetation.type",group="study.number", xlab = "Effect size (lnRR)", transfm = "none") +
  labs(y = "") + scale_y_discrete(labels = lab) +
  theme(
    axis.text.y = element_text(hjust = 0, angle = 0, colour = "black"),
    legend.direction = "horizontal", panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(),
    legend.position="bottom"
  ) +   scale_colour_manual(values = c("#125f86","#4299cb","#7cc7e4"), aesthetics = c("colour", "fill"))
main.plot_mon_vegetation_type
dev.off()

#using the function to calculate percentage of change
main.model_mon_vegetation_type_perc <- perc_change(main.model_mon_vegetation_type ,"vegetation.type")
main.model_mon_vegetation_type_perc
write.csv(main.model_mon_vegetation_type_perc, file = "vegetation_type_perc.csv", row.names = FALSE)


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
    soil.pH < 6.5 ~ "pH<6.5",
    soil.pH > 7.5 ~ "pH>7.5",
    TRUE ~ "pH 6.5-7.5"
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
  scale_y_continuous(limits = c(-1.5, 1.5), breaks = seq(-1, 1, by = 1)) +
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
  scale_y_continuous(limits = c(-1.5, 1.5), breaks = seq(-1, 1, by = 1)) +
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
main.plot_MBC <- orchard_plot(main.model_MBC_imp_rom,mod = "1",group="study.number",xlab = "Effect size (lnRR)", transfm = "none", g = FALSE) +
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
        ) 
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


tiff("F:/R/Meta-analysis-crop/MBC_analysis/control_monoculture_main_model_moderator_climate_zone.tiff", height = 8, width = 13, units = "cm", compression = "lzw", res = 600)
main.plot_MBC_climatic_zone <- orchard_plot(main.model_MBC_climatic_zone_1, mod = "climatic.zone_1",group="study.number",xlab = "Effect size (lnRR)", transfm = "none", g=FALSE, branch.size=0.8) +
  scale_x_discrete(limits = c("Temperate", "Subtropics", "Tropics"))+
  scale_y_continuous(limits = c(-1.2, 1.8),breaks = seq(-1, 1.5, 1))+
  scale_colour_manual(values = c("#36c179", "#3ed739", "#225a27"),aesthetics = c("colour", "fill")) + 
  labs(y = "") +
  theme(
    axis.text.y = element_text(hjust = 0, angle = 0, size = 12, colour = "black"),
    legend.direction = "horizontal", 
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(), 
    panel.grid.major.y = element_blank(),
    legend.position="none"
    ) + 
  annotate("text", x = 0.55, y = 1.2, label = "F = 10.413, p < 0.001", size = 2.5, fontface = "bold")
main.plot_MBC_climatic_zone$layers[[4]] <- NULL
summary_MBC_climatic <- summary(main.model_MBC_climatic_zone_1)
estimate <- summary_MBC_climatic$beta    # 点的数值
lowerCL <- summary_MBC_climatic$ci.lb       # 置信区间下限
upperCL <- summary_MBC_climatic$ci.ub       # 置信区间上限
point_data_climatic_MBC <- data.frame(
  condition = c( "Subtropics","Temperate", "Tropics"),       # x轴标签
  estimate = estimate,
  lowerCL = lowerCL,
  upperCL = upperCL          # 颜色值（和之前orchard一致）
)
main.plot_MBC_climatic_zone <- main.plot_MBC_climatic_zone+
  geom_point(
    data = point_data_climatic_MBC,  # 你自定义的数据
    aes(x = condition, y = estimate),
    shape = 21,
    size = 2.4,    
    fill = c("#36c179", "#3ed739", "#225a27"),  # ✅ 直接指定颜色，不映射
    color = "black"
  )
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
dat_MBC_combined$clim_zone <- factor(
  dat_MBC_combined$clim_zone,
  levels = c("Temperate", "Subtropics", "Tropics")
)

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
  scale_color_manual(values = c("Temperate" = "#0f5f67", "Subtropics" = "#618d10","Tropics" = "#a11d49")) +
  scale_fill_manual(values  = c("Temperate" = "#44b9c4", "Subtropics" = "#b4db66","Tropics" = "#e6778b") ) +
  labs(
    x = "LnRR-MBC",
    y = "LnRR-SOC",
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

model_stats <- dat_MBC_combined %>%
  group_by(clim_zone) %>%
  do({
    model <- lm(effect.size ~ effect.size.1, data = .)
    data.frame(
      intercept = coef(model)[1],
      slope = coef(model)[2],
      r_squared = summary(model)$r.squared,
      p_value = summary(model)$coefficients[2,4]
    )
  })

print(model_stats)

model_Climate <- lm(effect.size ~ effect.size.1 * clim_zone, data = dat_MBC_combined)
anova(model_Climate)

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

tiff("F:/R/Meta-analysis-crop/MBC_analysis/control_monoculture_main_model_moderator_age_perennial.tiff", height = 9.5, width = 13, units = "cm", compression = "lzw", res = 600)
main.plot_MBC_age_peren <- orchard_plot(main.model_MBC_age_peren, mod = "crop.age.classification",group="study.number", xlab = "Effect size (lnRR)", transfm = "none", g=FALSE, branch.size=0.8) +
  scale_colour_manual(values = c("#7e1d45","#c05579","#e6778b"), aesthetics = c("colour", "fill")) + 
  scale_y_continuous(limits = c(-1.2, 1.8),breaks = seq(-1, 1.5, 1)) +
  labs(y = "") + 
  theme(
    axis.text.y = element_text(hjust = 0, angle = 0, size = 12, colour = "black"),
    legend.direction = "horizontal", panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(),
    legend.position="bottom"
  ) + 
  annotate("text", x = 0.55, y = 1.2, label = "F = 19.264, p < 0.001", size = 2.5, fontface = "bold")
main.plot_MBC_age_peren$layers[[4]] <- NULL
summary_MBC_age <- summary(main.model_MBC_age_peren)
estimate <- summary_MBC_age$beta    # 点的数值
lowerCL <- summary_MBC_age$ci.lb       # 置信区间下限
upperCL <- summary_MBC_age$ci.ub       # 置信区间上限
point_data_age_MBC <- data.frame(
  condition = c("Age>10", "Age 5-10", "Age<5"),       # x轴标签
  estimate = estimate,
  lowerCL = lowerCL,
  upperCL = upperCL          # 颜色值（和之前orchard一致）
)
main.plot_MBC_age_peren <- main.plot_MBC_age_peren +
  geom_point(
    data = point_data_age_MBC,  # 你自定义的数据
    aes(x = condition, y = estimate),
    shape = 21,
    size = 2.4,    
    fill = c("#7e1d45","#c05579","#e6778b"),  # ✅ 直接指定颜色，不映射
    color = "black"
  )
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
dat_MBC_combined_age$crop_age <- factor(
  dat_MBC_combined_age$crop_age,
  levels = c("Age<5", "Age 5-10", "Age>10")
)

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
  scale_color_manual(values = c("Age<5" = "#0f5f67", "Age 5-10" = "#618d10","Age>10" = "#a11d49")) +
  scale_fill_manual(values  = c("Age<5" = "#44b9c4", "Age 5-10" = "#b4db66","Age>10" = "#e6778b") ) +
  labs(
    x = "LnRR-MBC",
    y = "LnRR-SOC",
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

model_stats_age <- dat_MBC_combined_age %>%
  group_by(crop_age) %>%
  do({
    model <- lm(effect.size ~ effect.size.1, data = .)
    data.frame(
      intercept = coef(model)[1],
      slope = coef(model)[2],
      r_squared = summary(model)$r.squared,
      p_value = summary(model)$coefficients[2,4]
    )
  })

print(model_stats_age)

model_Age <- lm(effect.size ~ effect.size.1 * crop_age, data = dat_MBC_combined_age)
anova(model_Age)

################################################################################
# soil pH####3
################################################################################

#######set different category color

dat_MBC_imp_rom$soil.pH.classification <- as.factor(dat_MBC_imp_rom$soil.pH.classification)
dat_MBC_imp_rom$soil.pH.classification <- factor(dat_MBC_imp_rom$soil.pH.classification, levels = c("pH>7.5", "pH 6.5-7.5", "pH<6.5"))

main.model_MBC_soil_pH <- rma.mv(yi = yi, V = vi, mod = ~ soil.pH.classification - 1, random = ~ 1 | site_id / effect_size_id, data = dat_MBC_imp_rom, method = "REML", test = "t")
main.model_MBC_soil_pH
anova(main.model_MBC_soil_pH)

lab <- c(
  "Soil.pH.classification<6.5" = "< 6.5",
  "Soil.pH.classification>7.5" = "> 7.5",
  "Soil.pH.classification6.5-7.5" = "6.5 - 7.5"
)


tiff("F:/R/Meta-analysis-crop/MBC_analysis/control_monoculture_main_model_moderator_soil_pH.tiff", height = 8, width = 13, units = "cm", compression = "lzw", res = 600)
main.plot_MBC_soil_pH <- orchard_plot(main.model_MBC_soil_pH,mod = "soil.pH.classification",group="study.number",xlab = "Effect size (lnRR)", transfm = "none", g=FALSE, branch.size=0.8) +
  scale_x_discrete(labels = c("pH>7.5", "pH 6.5–7.5", "pH<6.5"))+
  scale_y_continuous(limits = c(-1.2, 1.8),breaks = seq(-1, 1.5, 1))+
  scale_colour_manual(values = c("#125f86","#4299cb","#7cc7e4"), aesthetics = c("colour", "fill")) + 
  labs(y = "") +
  theme(
    axis.text.y = element_text(hjust = 0, angle = 0, size = 12, colour = "black"),
    legend.direction = "horizontal", panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(),
    legend.position="none"
  )+
  annotate("text", x = 0.55, y = 1.2, label = "F = 9.2662, p < 0.001", size = 2.5, fontface = "bold")
main.plot_MBC_soil_pH$layers[[4]] <- NULL
summary_MBC_pH <- summary(main.model_MBC_soil_pH)
estimate <- summary_MBC_pH$beta    # 点的数值
lowerCL <- summary_MBC_pH$ci.lb       # 置信区间下限
upperCL <- summary_MBC_pH$ci.ub       # 置信区间上限
point_data_age <- data.frame(
  condition = c("PH>7.5", "PH 6.5-7.5", "PH<6.5"),       # x轴标签
  estimate = estimate,
  lowerCL = lowerCL,
  upperCL = upperCL          # 颜色值（和之前orchard一致）
)
main.plot_MBC_soil_pH <- main.plot_MBC_soil_pH +
  geom_point(
    data = point_data_age,  # 你自定义的数据
    aes(x = condition, y = estimate),
    shape = 21,
    size = 2.4,    
    fill = c("#125f86","#4299cb","#7cc7e4"),  # ✅ 直接指定颜色，不映射
    color = "black")
main.plot_MBC_soil_pH
dev.off()

#using the function to calculate percentage of change
main.model_MBC_soil_pH_perc <- perc_change(main.model_MBC_soil_pH ,"soil.pH.classification")
main.model_MBC_soil_pH_perc
write.csv(main.model_MBC_soil_pH_perc, file = "soil_pH_perc.csv", row.names = FALSE)

#####soil.pH.classification
dat_MBC_imp_rom_pH_7.5 <- dat_MBC_imp_rom %>%
  filter(soil.pH.classification == "pH>7.5") %>%
  droplevels()

dat_MBC_imp_rom_pH_7 <- dat_MBC_imp_rom %>%
  filter(soil.pH.classification == "pH 6.5-7.5") %>%
  droplevels()

dat_MBC_imp_rom_pH_6.5 <- dat_MBC_imp_rom %>%
  filter(soil.pH.classification == "pH<6.5") %>%
  droplevels()


dat_MBC_imp_rom_pH_6.5$crop_pH <- "pH<6.5"
dat_MBC_imp_rom_pH_7$crop_pH <- "pH 6.5-7.5"
dat_MBC_imp_rom_pH_7.5$crop_pH <- "pH>7.5"

dat_MBC_combined_pH <- rbind(dat_MBC_imp_rom_pH_6.5, dat_MBC_imp_rom_pH_7,dat_MBC_imp_rom_pH_7.5)
dat_MBC_combined_pH$crop_pH <- factor(
  dat_MBC_combined_pH$crop_pH,
  levels = c("pH<6.5", "pH 6.5-7.5", "pH>7.5")
)

filename <- paste('Result_Regression/',"SOC_MBC_pH",'.jpeg', sep='')
jpeg(filename,width = 2000,height=2000,res = 300)
main.model_MBC_pH <- ggplot(dat_MBC_combined_pH, aes(x = effect.size.1, y = effect.size, color = crop_pH, fill = crop_pH)) +
  geom_hline(yintercept = 0, linetype = "longdash", color = "grey80", size = 0.8) +
  geom_vline(xintercept = 0, linetype = "longdash", color = "grey80", size = 0.8) +
  geom_point(shape = 21, size = 3, alpha = 0.8, color = "black") +
  geom_smooth(method = "lm", se = TRUE, size = 1.2, alpha = 0.3) +
  scale_color_manual(values = c("pH<6.5" = "#0f5f67", "pH 6.5-7.5" = "#618d10","pH>7.5" = "#a11d49")) +
  scale_fill_manual(values  = c("pH<6.5" = "#44b9c4", "pH 6.5-7.5" = "#b4db66","pH>7.5" = "#e6778b")) +
  stat_poly_eq(
    aes(label = paste("bold(", ..rr.label.., ")")),
    label.x.npc = "left", label.y.npc = "top",
    formula = y ~ x,
    parse = TRUE,
    size = 4,
    fontface = "bold", 
    show.legend = FALSE
  )+
  labs(
    x = "LnRR-MBC",
    y = "LnRR-SOC",
    title = "Regression of Different Soil pH",
    color = "Soil pH",
    fill = "Soil pH"
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
print(main.model_MBC_pH)
graphics.off()

model_stats_pH <- dat_MBC_combined_pH %>%
  group_by(crop_pH) %>%
  do({
    model <- lm(effect.size ~ effect.size.1, data = .)
    data.frame(
      intercept = coef(model)[1],
      slope = coef(model)[2],
      r_squared = summary(model)$r.squared,
      p_value = summary(model)$coefficients[2,4]
    )
  })

print(model_stats_pH)

model_pH <- lm(effect.size ~ effect.size.1 * crop_pH, data = dat_MBC_combined_pH)
anova(model_pH)




# vegetation type
main.model_MBC_vegetation_type <- rma.mv(yi = yi, V = vi, mod = ~ vegetation.type - 1, random = ~ 1 | site_id / effect_size_id, data = dat_MBC_imp_rom, method = "REML", test = "t")
main.model_MBC_vegetation_type

lab <- c(
  "Vegetation.typeHerbaceous" = "Herbaceous",
  "Vegetation.typeWoody" = "Woody"
)

tiff("F:/R/Meta-analysis-crop/control_monoculture_main_model_moderator_vegetation_type.tiff", height = 11, width = 11, units = "cm", compression = "lzw", res = 600)
main.plot_MBC_vegetation_type <- orchard_plot(main.model_MBC_vegetation_type,mod = "vegetation.type",group="study.number", xlab = "Effect size (lnRR)", transfm = "none") +
  labs(y = "") + scale_y_discrete(labels = lab) +
  theme(
    axis.text.y = element_text(hjust = 0, angle = 0, colour = "black"),
    legend.direction = "horizontal", panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(),
    legend.position="bottom"
  ) +   scale_colour_manual(values = c("#125f86","#4299cb","#7cc7e4"), aesthetics = c("colour", "fill"))
main.plot_MBC_vegetation_type
dev.off()

#using the function to calculate percentage of change
main.model_mon_vegetation_type_perc <- perc_change(main.model_mon_vegetation_type ,"vegetation.type")
main.model_mon_vegetation_type_perc
write.csv(main.model_mon_vegetation_type_perc, file = "vegetation_type_perc.csv", row.names = FALSE)
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

ggsave("F:/R/Meta-analysis-crop/general_figure_perennial_MBC.tiff", dpi = 1000, compression = "lzw", plot = Glo.Dis.F, width = 18.3, height = 7.95)



