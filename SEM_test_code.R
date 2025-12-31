setwd("F:/R/Meta-analysis-crop/SEM")

library(lme4)
library(lavaan)
library(semPlot)
library(car)
library(plspm)

dat_SOC <- dat_sta_fin_SOC

dat_SOC_MBC <- read.csv("dat_SOC_MBC.csv", header = TRUE)

dat_numeric <- dat_SOC_MBC %>%
  rename(SOC.effect.size = effect.size, MBC.effect.size = effect.size.1) %>%
  select(SOC.effect.size, MBC.effect.size, MAT, MAP, elevation, soil.pH, clay, bulk.density, soc.peren.mean_30cm, age.peren, leg.peren, vegetation.type) %>%
  mutate(Legume = ifelse(leg.peren == "Yes", 1, 0)) %>%
  select(-leg.peren)

# 检查变量结构
str(dat_numeric)

# 找出所有需要标准化的连续变量
continuous_vars <- c("SOC.effect.size", "MBC.effect.size", "MAT", "MAP", "soil.pH", "clay", "bulk.density", "age.peren")

# 对连续变量进行标准化处理
dat_numeric[continuous_vars] <- scale(dat_numeric[continuous_vars])

# 检查标准化后的数据
str(dat_numeric)

# 显示部分数据
head(dat_numeric)

# blocks 定义
blocks <- list(
  c("age.peren", "Legume"),                   # Plant
  c("MAT", "MAP"),                              # Climate
  c("soil.pH"),                                             # Soil chemistry
  c("clay", "bulk.density"),                                # Soil physics
  c("MBC.effect.size"),                                     # MBC
  c("SOC.effect.size")                                      # SOC
)
# 正确路径矩阵
#           Plant Soil Climate MBC SOC
path_matrix <- rbind(
  Plant        =   c(0,  0,  0,  0,  0,  0),
  Climate      =   c(0,  0,  0,  0,  0,  0),
  Soil_chem    =   c(0,  0,  0,  0,  0,  0),
  Soil_phys    =   c(0,  0,  0,  0,  0,  0),
  MBC          =   c(1,  1,  1,  1,  0,  0),
  SOC          =   c(1,  1,  1,  1,  1,  0)
)
colnames(path_matrix) <- rownames(path_matrix) <- c("Plant", "Climate", "Soil_chem", "Soil_phys", "MBC", "SOC")

modes <- c("A", "A", "A", "A", "A", "A")

# 运行模型（dat_numeric_clean 是你之前处理好的数据）
pls_model <- plspm(dat_numeric, path_matrix, blocks, modes = modes)

# 查看结果
summary(pls_model)

# 可视化路径图
plot(pls_model)


pls_model_boot <- plspm(dat_numeric, path_matrix, blocks, 
                        modes = modes, boot.val = TRUE, br = 500)
# 提取Bootstrap路径结果
boot_paths <- pls_model_boot$boot$paths

# 确认路径名存储在行名中，并转换为列
boot_paths$Path <- rownames(boot_paths)

# 计算 p 值
boot_paths$p_value <- 2 * (1 - pnorm(abs(boot_paths$Original / boot_paths$Std.Error)))

# 添加显著性星号标记
boot_paths$signif <- cut(boot_paths$p_value,
                         breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                         labels = c("***", "**", "*", ""),
                         right = FALSE)

# 最终清晰输出
results_star <- data.frame(
  Path = boot_paths$Path,
  Coefficient = round(boot_paths$Original, 3),
  p_value = round(boot_paths$p_value, 4),
  Significance = boot_paths$signif
)

# 输出结果
print(results_star)


##############SOC-MBC
##########################




blocks <- list(
  c("age.peren", "Legume"),                   # Plant
  c("MAT", "MAP"),                              # Climate
  c("soil.pH"),                                             # Soil chemistry
  c("clay", "bulk.density"),                                # Soil physics
  c("SOC.effect.size"),                                     # SOC
  c("MBC.effect.size")                                      # MBC
)
# 正确路径矩阵
#  data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAIAAAAUCAYAAACnOeyiAAAAXklEQVR4XmNgAIL///8zMYSGhjIDGYIMIiIMvECGMwMDN4M4kFEDUqIIZKwDMdSBjAsghj6Q8QPEMAAy/lOBoQekv4AYKkDGfgZeXl4RICOLQUtLiw3IUAJJMQIZ7AC2tU2tJWCy/wAAAABJRU5ErkJggg==         Plant Soil Climate MBC SOC
path_matrix <- rbind(
  Plant        =   c(0,  0,  0,  0,  0,  0),
  Climate      =   c(0,  0,  0,  0,  0,  0),
  Soil_chem    =   c(0,  0,  0,  0,  0,  0),
  Soil_phys    =   c(0,  0,  0,  0,  0,  0),
  SOC          =   c(1,  1,  1,  1,  0,  0),
  MBC          =   c(1,  1,  1,  1,  1,  0)
)
colnames(path_matrix) <- rownames(path_matrix) <- c("Plant", "Climate", "Soil_chem", "Soil_phys","SOC", "MBC")

modes <- c("A", "A", "A", "A", "A", "A")

# 运行模型（dat_numeric_clean 是你之前处理好的数据）
pls_model <- plspm(dat_numeric, path_matrix, blocks, modes = modes)

# 查看结果
summary(pls_model)

# 可视化路径图
plot(pls_model)


pls_model_boot <- plspm(dat_numeric, path_matrix, blocks, 
                        modes = modes, boot.val = TRUE, br = 500)
# 提取Bootstrap路径结果
boot_paths <- pls_model_boot$boot$paths

# 确认路径名存储在行名中，并转换为列
boot_paths$Path <- rownames(boot_paths)

# 计算 p 值
boot_paths$p_value <- 2 * (1 - pnorm(abs(boot_paths$Original / boot_paths$Std.Error)))

# 添加显著性星号标记
boot_paths$signif <- cut(boot_paths$p_value,
                         breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                         labels = c("***", "**", "*", ""),
                         right = FALSE)

# 最终清晰输出
results_star <- data.frame(
  Path = boot_paths$Path,
  Coefficient = round(boot_paths$Original, 3),
  p_value = round(boot_paths$p_value, 4),
  Significance = boot_paths$signif
)

# 输出结果
print(results_star)


