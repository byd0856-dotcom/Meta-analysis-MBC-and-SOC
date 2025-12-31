
# 加载必要包
library(ggplot2)
library(patchwork)  # 可选：用于拼图显示两个图

# 读取数据
dat_SOC_MBC <- read.csv("dat_SOC_MBC.csv", header = TRUE)
dat_SOC_MBC <- dat_SOC_MBC %>%
  rename(SOC.effect.size = effect.size, MBC.effect.size = effect.size.1)

# 图1: SOC.effect.size vs soc.peren.mean_30cm
p1 <- ggplot(dat_SOC_MBC, aes(x = soc.peren.mean_30cm, y = SOC.effect.size)) +
  geom_point(color = "#dd5353", size = 2.5, alpha = 0.7) +
  geom_smooth(method = "lm", color = "#de2d26", fill="#ef7d90",se = TRUE) +
  geom_hline(yintercept=0,linetype="longdash",size=0.8,colour="grey65")+
  labs(x = "Initial SOC(Mg/ha)",
       y = "LnRR-SOC",
       title = "LnRR-SOC vs SOC initial") +
  theme(
    aspect.ratio = 0.9,
    panel.background = element_rect(fill = 'white', colour = 'white'),
    plot.title = element_text(hjust = 0.5, size=20, face=2,margin = ggplot2::margin(0,0,0.5,0,'cm')),
    axis.title=element_text(colour='black', size=14,face=2),
    axis.text.x = element_text(colour = 'black', size = 12,margin = ggplot2::margin(0,0,0.3,0,'cm')),
    axis.text.y = element_text(colour = 'black', size = 12,margin = ggplot2::margin(0,0,0,0.3,'cm')),
    axis.line = element_line(colour = 'black',size=1),
    legend.position = c(0.275, 0.88),
    legend.title = element_blank(),
    legend.text  = element_text(face = "bold", size = 10)
  )

# 图2: MBC.effect.size vs soc.peren.mean_30cm
p2 <- ggplot(dat_SOC_MBC, aes(x = soc.peren.mean_30cm, y = MBC.effect.size)) +
  geom_point(color = "#6baed6", size = 2.5, alpha = 0.7) +
  geom_smooth(method = "lm", color = "#08519c", fill="#c6dbef",se = TRUE) +
  geom_hline(yintercept=0,linetype="longdash",size=0.8,colour="grey65")+
  labs(x = "Initial SOC(Mg/ha)",
       y = "LnRR-MBC",
       title = "LnRR-MBC vs SOC initial") +
  theme(
    aspect.ratio = 0.9,
    panel.background = element_rect(fill = 'white', colour = 'white'),
    plot.title = element_text(hjust = 0.5, size=20, face=2,margin = ggplot2::margin(0,0,0.5,0,'cm')),
    axis.title=element_text(colour='black', size=14,face=2),
    axis.text.x = element_text(colour = 'black', size = 12,margin = ggplot2::margin(0,0,0.3,0,'cm')),
    axis.text.y = element_text(colour = 'black', size = 12,margin = ggplot2::margin(0,0,0,0.3,'cm')),
    axis.line = element_line(colour = 'black',size=1),
    legend.position = c(0.275, 0.88),
    legend.title = element_blank(),
    legend.text  = element_text(face = "bold", size = 10)
  )

model_soc <- lm(SOC.effect.size ~ soc.peren.mean_30cm, data = dat_SOC_MBC)
r2_soc <- summary(model_soc)$r.squared
p_soc <- summary(model_soc)$coefficients[2,4]

model_mbc <- lm(MBC.effect.size ~ soc.peren.mean_30cm, data = dat_SOC_MBC)
r2_mbc <- summary(model_mbc)$r.squared
p_mbc <- summary(model_mbc)$coefficients[2,4]


p1 <- p1 + annotate("text",
                    x = max(dat_SOC_MBC$soc.peren.mean_30cm, na.rm = TRUE),
                    y = max(dat_SOC_MBC$SOC.effect.size, na.rm = TRUE),
                    label = paste0("R² = ", round(r2_soc, 3), 
                                   "\np = ",ifelse(p_mbc < 1e-3, "0.002", paste0("= ", signif(p_mbc, 3)))),
                    hjust = 1, vjust = 1, size = 5)

p2 <- p2 + annotate("text",
                    x = max(dat_SOC_MBC$soc.peren.mean_30cm, na.rm = TRUE),
                    y = max(dat_SOC_MBC$MBC.effect.size, na.rm = TRUE),
                    label = paste0("R² = ", round(r2_mbc, 3),
                                   "\np ", ifelse(p_mbc < 1e-3, "< 0.001", paste0("= ", signif(p_mbc, 3)))),
                    hjust = 1, vjust = 1, size = 5)
# 拼图显示
filename <- paste('Result_Regression/',"SOC_initial",'.jpeg', sep='')
jpeg(filename,width = 4000,height=2000,res = 300)
SOC_initial <- p1 + p2
print(SOC_initial)
graphics.off()

#####Temperate
dat_MBC_imp_rom_50 <- dat_MBC_imp_rom %>%
  filter(SOC.initial == "SOC<50") %>%
  droplevels()

#####Subtropics
dat_MBC_imp_rom_75 <- dat_MBC_imp_rom %>%
  filter(SOC.initial == "SOC 50-100") %>%
  droplevels()

#####Tropics

dat_MBC_imp_rom_100 <- dat_MBC_imp_rom %>%
  filter(SOC.initial == "SOC>100") %>%
  droplevels()

dat_MBC_imp_rom_50$SOC_initial <- "SOC<50"
dat_MBC_imp_rom_75$SOC_initial <- "SOC 50-100"
dat_MBC_imp_rom_100$SOC_initial <- "SOC>100"

dat_MBC_combined_initial <- rbind(dat_MBC_imp_rom_50, dat_MBC_imp_rom_75,dat_MBC_imp_rom_100)
dat_MBC_combined_initial$SOC_initial <- factor(
  dat_MBC_combined_initial$SOC_initial,
  levels = c("SOC<50", "SOC 50-100", "SOC>100")
)

filename <- paste('Result_Regression/',"SOC_MBC_initialSOC",'.jpeg', sep='')
jpeg(filename,width = 2000,height=2000,res = 300)
main.model_MBC_SOCinitial <- ggplot(dat_MBC_combined_initial, aes(x = effect.size.1, y = effect.size, color = SOC_initial, fill = SOC_initial)) +
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
  scale_color_manual(values = c("SOC<50" = "#0f5f67", "SOC 50-100" = "#618d10","SOC>100" = "#a11d49")) +
  scale_fill_manual(values  = c("SOC<50" = "#44b9c4", "SOC 50-100" = "#b4db66","SOC>100" = "#e6778b") ) +
  labs(
    x = "LnRR of MBC",
    y = "LnRR-SOC",
    title = "Regression of Different Initial SOC",
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
    legend.position = c(0.276, 0.88),
    legend.title = element_blank(),
    legend.text  = element_text(face = "bold", size = 10)
  )
print(main.model_MBC_SOCinitial)
graphics.off()

model_stats_initial <- dat_MBC_combined_initial %>%
  group_by(SOC_initial) %>%
  do({
    model <- lm(effect.size ~ effect.size.1, data = .)
    data.frame(
      intercept = coef(model)[1],
      slope = coef(model)[2],
      r_squared = summary(model)$r.squared,
      p_value = summary(model)$coefficients[2,4]
    )
  })

print(model_stats_initial)

model_initial <- lm(effect.size ~ effect.size.1 * SOC_initial, data = dat_MBC_combined_initial)
anova(model_initial)

library(ggplot2)
library(dplyr)
library(patchwork)

# 假设你数据框叫 dat_SOC_MBC
dat_SOC_MBC$ln_MBC <- log(dat_SOC_MBC$MBC.actual.peren.sameunit)

# --------- 函数：拟合、剔除异常值、绘图 ----------
plot_poly_with_outlier_removal <- function(data, y_var, title_label,
                                           y_label = NULL,
                                           color_point = "#1f77b4",
                                           line_color = "red",
                                           fill_color = "pink") {
  if (is.null(y_label)) y_label <- y_var
  
  # 构建公式
  formula <- as.formula(paste(y_var, "~ poly(ln_MBC, 2, raw=TRUE)"))
  
  # 拟合初始模型
  model <- lm(formula, data = data)
  
  # 提取标准化残差
  resid_std <- rstandard(model)
  
  # 剔除异常值
  data_clean <- data[abs(resid_std) <= 2, ]
  
  # 重新拟合模型
  model_clean <- lm(formula, data = data_clean)
  
  # 提取R²和p值
  r2 <- summary(model_clean)$r.squared
  pval <- summary(model_clean)$coefficients[2, 4]
  
  y_values <- data[[y_var]]  # 新增，避免aes环境下y_var无法引用
  
  # 绘图
  p <- ggplot(data_clean, aes(x = ln_MBC, y = .data[[y_var]])) +
    geom_point(color = color_point, size = 2.5, alpha = 0.7) +
    stat_smooth(method = "lm", 
                se = TRUE, color = line_color, fill = fill_color, size = 1.2) +
    geom_hline(yintercept=0,linetype="longdash",size=0.8,colour="grey65")+
    annotate("text",
             x = max(data_clean$ln_MBC, na.rm = TRUE) - 0.2,
             y = max(y_values[abs(resid_std) <= 2], na.rm = TRUE),
             label = paste0("R² = ", round(r2, 3), "\np = ", signif(pval, 3)),
             hjust = 1, vjust = 1,
             size = 5) +
    labs(
      x = "ln(MBC initial)",
      y = y_label,
      title = title_label
    ) +
    theme(
      aspect.ratio = 0.9,
      panel.background = element_rect(fill = 'white', colour = 'white'),
      plot.title = element_text(hjust = 0.5, size=20, face=2,margin = ggplot2::margin(0,0,0.5,0,'cm')),
      axis.title=element_text(colour='black', size=14,face=2),
      axis.text.x = element_text(colour = 'black', size = 12,margin = ggplot2::margin(0,0,0.3,0,'cm')),
      axis.text.y = element_text(colour = 'black', size = 12,margin = ggplot2::margin(0,0,0,0.3,'cm')),
      axis.line = element_line(colour = 'black',size=1),
      legend.position = c(0.275, 0.88),
      legend.title = element_blank(),
      legend.text  = element_text(face = "bold", size = 10)
    )
  
  return(p)
}


# --------- 分别画出两个图 ----------
# 图1：SOC响应图（蓝色点，红色线，粉色SE）
p_soc <- plot_poly_with_outlier_removal(
  dat_SOC_MBC,
  y_var = "SOC.effect.size",
  y_label = "LnRR-SOC",
  title_label = "LnRR-SOC vs ln(MBC initial)",
  color_point = "#e6550d",
  line_color = "#a63603",
  fill_color = "#fcbba1"
)

# --------- GAM函数：拟合、剔除异常值、绘图 ----------

plot_gam_with_outlier_removal <- function(data, y_var, title_label,
                                          y_label = NULL,
                                          color_point = "#1f77b4",
                                          line_color = "red",
                                          fill_color = "pink") {
  
  if (is.null(y_label)) y_label <- y_var
  
  # 构建 GAM 模型公式
  formula <- as.formula(paste(y_var, "~ s(ln_MBC)"))
  
  # 初次拟合 GAM 模型
  model <- gam(formula, data = data)
  
  # 提取残差（GAM 没有 rstandard，用原始残差近似处理）
  resid_raw <- resid(model)
  
  # 剔除残差大于 2 的异常值
  data_clean <- data[abs(resid_raw) <= 2, ]
  
  # 重新拟合 GAM 模型
  model_clean <- gam(formula, data = data_clean)
  
  # 提取 R² 和整体平滑项的 p 值
  r2 <- summary(model_clean)$r.sq
  pval <- summary(model_clean)$s.table[1, 4]
  
  y_values <- data[[y_var]]
  
  # 绘图
  p <- ggplot(data_clean, aes(x = ln_MBC, y = .data[[y_var]])) +
    geom_point(color = color_point, size = 2.5, alpha = 0.7) +
    stat_smooth(method = "gam", formula = y ~ s(x),
                se = TRUE, color = line_color, fill = fill_color, size = 1.2) +
    geom_hline(yintercept = 0, linetype = "longdash", size = 0.8, colour = "grey65") +
    annotate("text",
             x = max(data_clean$ln_MBC, na.rm = TRUE) - 0.2,
             y = max(y_values[abs(resid_raw) <= 2], na.rm = TRUE),
             label = paste0("R² = ", round(r2, 3), "\np = ", signif(pval, 3)),
             hjust = 1, vjust = 1,
             size = 5) +
    labs(
      x = "ln(MBC initial)",
      y = y_label,
      title = title_label
    ) +
    theme(
      aspect.ratio = 0.9,
      panel.background = element_rect(fill = 'white', colour = 'white'),
      plot.title = element_text(hjust = 0.5, size = 20, face = 2, margin = ggplot2::margin(0, 0, 0.5, 0, 'cm')),
      axis.title = element_text(colour = 'black', size = 14, face = 2),
      axis.text.x = element_text(colour = 'black', size = 12, margin = ggplot2::margin(0, 0, 0.3, 0, 'cm')),
      axis.text.y = element_text(colour = 'black', size = 12, margin = ggplot2::margin(0, 0, 0, 0.3, 'cm')),
      axis.line = element_line(colour = 'black', size = 1),
      legend.position = c(0.275, 0.88),
      legend.title = element_blank(),
      legend.text = element_text(face = "bold", size = 10)
    )
  
  return(p)
}


# 图2：MBC响应图（绿色点，深绿色线，浅绿SE）
p_mbc <- plot_poly_with_outlier_removal(
  dat_SOC_MBC,
  y_var = "MBC.effect.size",
  y_label = "LnRR-MBC",
  title_label = "LnRR-MBC vs ln(MBC initial)",
  color_point = "#32b671",
  line_color = "#145742",
  fill_color = "#a3d8ba")

# 并排展示

filename <- paste('Result_Regression/',"MBC_initial_2",'.jpeg', sep='')
jpeg(filename,width = 4000,height=2000,res = 300)
MBC_initial <- p_soc + p_mbc
print(MBC_initial)
graphics.off()

