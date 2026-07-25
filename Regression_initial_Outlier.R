library(ggplot2)
library(dplyr)
library(patchwork)

# 读取数据
dat_SOC_MBC <- read.csv("dat_SOC_MBC.csv", header = TRUE)

dat_SOC_MBC <- dat_SOC_MBC %>%
  rename(
    SOC.effect.size = effect.size,
    MBC.effect.size = effect.size.1
  )

# 识别高 Initial SOC 点：默认 Initial SOC 最大值为 high SOC point
high_soc_value <- max(dat_SOC_MBC$soc.peren.mean_30cm, na.rm = TRUE)

dat_SOC_MBC <- dat_SOC_MBC %>%
  mutate(
    high_SOC_point = soc.peren.mean_30cm == high_soc_value
  )

dat_no_high_SOC <- dat_SOC_MBC %>%
  filter(!high_SOC_point)

# -------------------------
# 模型拟合，用于输出表格
# -------------------------

model_soc_all <- lm(SOC.effect.size ~ soc.peren.mean_30cm, data = dat_SOC_MBC)
model_soc_no  <- lm(SOC.effect.size ~ soc.peren.mean_30cm, data = dat_no_high_SOC)

model_mbc_all <- lm(MBC.effect.size ~ soc.peren.mean_30cm, data = dat_SOC_MBC)
model_mbc_no  <- lm(MBC.effect.size ~ soc.peren.mean_30cm, data = dat_no_high_SOC)

result_table <- data.frame(
  Model = c(
    "LnRR-SOC with high SOC point",
    "LnRR-SOC without high SOC point",
    "LnRR-MBC with high SOC point",
    "LnRR-MBC without high SOC point"
  ),
  Slope = c(
    coef(model_soc_all)[2],
    coef(model_soc_no)[2],
    coef(model_mbc_all)[2],
    coef(model_mbc_no)[2]
  ),
  R2 = c(
    summary(model_soc_all)$r.squared,
    summary(model_soc_no)$r.squared,
    summary(model_mbc_all)$r.squared,
    summary(model_mbc_no)$r.squared
  ),
  p_value = c(
    summary(model_soc_all)$coefficients[2, 4],
    summary(model_soc_no)$coefficients[2, 4],
    summary(model_mbc_all)$coefficients[2, 4],
    summary(model_mbc_no)$coefficients[2, 4]
  )
)

print(result_table)

write.csv(
  result_table,
  "Initial_SOC_sensitivity_analysis_table.csv",
  row.names = FALSE
)

# -------------------------
# 统一主题
# -------------------------

my_theme <- theme_classic() +
  theme(
    aspect.ratio = 0.85,
    axis.line = element_line(colour = "black", linewidth = 0.8),
    axis.title = element_text(size = 12, colour = "black", face = "bold"),
    axis.text = element_text(size = 10, colour = "black"),
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 9),
    legend.key = element_blank()
  )

# -------------------------
# A. LnRR-SOC
# -------------------------

p1 <- ggplot(dat_SOC_MBC, aes(x = soc.peren.mean_30cm, y = SOC.effect.size)) +
  
  geom_point(
    data = dat_SOC_MBC %>% filter(!high_SOC_point),
    aes(color = "Data"),
    shape = 16,
    size = 2.2,
    alpha = 0.75
  ) +
  
  geom_point(
    data = dat_SOC_MBC %>% filter(high_SOC_point),
    aes(color = "High SOC point"),
    shape = 16,
    size = 3.2,
    alpha = 0.95
  ) +
  
  geom_smooth(
    aes(color = "With high SOC point"),
    method = "lm",
    se = TRUE,
    fill = "#f4a3b4",
    linewidth = 1.0,
    alpha = 0.35
  ) +
  
  geom_smooth(
    data = dat_no_high_SOC,
    aes(
      x = soc.peren.mean_30cm,
      y = SOC.effect.size,
      color = "Without high SOC point"
    ),
    method = "lm",
    se = FALSE,
    linetype = "dashed",
    linewidth = 1.0
  ) +
  
  geom_hline(
    yintercept = 0,
    linetype = "longdash",
    linewidth = 0.8,
    colour = "grey55"
  ) +
  
  scale_color_manual(
    values = c(
      "Data" = "#d95f5f",
      "High SOC point" = "#b2182b",
      "With high SOC point" = "#de2d26",
      "Without high SOC point" = "#a50f15"
    ),
    breaks = c(
      "Data",
      "High SOC point",
      "With high SOC point",
      "Without high SOC point"
    )
  ) +
  
  labs(
    x = expression("Initial SOC (Mg ha"^-1*")"),
    y = "Log Response Ratio of SOC",
    title = "lnRR-SOC vs Initial SOC"
  ) +
  
  my_theme

# -------------------------
# B. LnRR-MBC
# -------------------------

p2 <- ggplot(dat_SOC_MBC, aes(x = soc.peren.mean_30cm, y = MBC.effect.size)) +
  
  geom_point(
    data = dat_SOC_MBC %>% filter(!high_SOC_point),
    aes(color = "Data"),
    shape = 16,
    size = 2.2,
    alpha = 0.75
  ) +
  
  geom_point(
    data = dat_SOC_MBC %>% filter(high_SOC_point),
    aes(color = "High SOC point"),
    shape = 16,
    size = 3.2,
    alpha = 0.95
  ) +
  
  geom_smooth(
    aes(color = "With high SOC point"),
    method = "lm",
    se = TRUE,
    fill = "#c6dbef",
    linewidth = 1.0,
    alpha = 0.35
  ) +
  
  geom_smooth(
    data = dat_no_high_SOC,
    aes(
      x = soc.peren.mean_30cm,
      y = MBC.effect.size,
      color = "Without high SOC point"
    ),
    method = "lm",
    se = FALSE,
    linetype = "dashed",
    linewidth = 1.0
  ) +
  
  geom_hline(
    yintercept = 0,
    linetype = "longdash",
    linewidth = 0.8,
    colour = "grey55"
  ) +
  
  scale_color_manual(
    values = c(
      "Data" = "#6baed6",
      "High SOC point" = "#2171b5",
      "With high SOC point" = "#08519c",
      "Without high SOC point" = "#08306b"
    ),
    breaks = c(
      "Data",
      "High SOC point",
      "With high SOC point",
      "Without high SOC point"
    )
  ) +
  
  labs(
    x = expression("Initial SOC (Mg ha"^-1*")"),
    y = "Log Response Ratio of MBC",
    title = "lnRR-MBC vs Initial SOC"
  ) +
  
  my_theme

# -------------------------
# 合并图
# -------------------------

p_all <- p1 + p2 + plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.title = element_blank()
  )

p_all

# -------------------------
# 保存图片
# -------------------------

ggsave(
  filename = "Initial_SOC_sensitivity_analysis.png",
  plot = p_all,
  width = 10,
  height = 5,
  dpi = 600
)

ggsave(
  filename = "Initial_SOC_sensitivity_analysis.pdf",
  plot = p_all,
  width = 10,
  height = 5
)