# 加载必要的库
install.packages("randomForest")
install.packages("caret")
install.packages("dplyr")
library(randomForest)
library(caret)
library(dplyr)

# 假设您已加载数据
dat_SOC <- dat_sta_fin_SOC

dat_SOC_MBC <- read.csv("dat_SOC_MBC.csv", header = TRUE)
dat_SOC_MBC <- dat_SOC_MBC %>%
  rename(SOC.effect.size = effect.size, MBC.effect.size = effect.size.1)%>%
  filter(!is.na(MBC.actual.peren.sameunit))

###########SOC
# 选择感兴趣的特征变量（连续变量和分类变量）
# 筛选连续变量和分类变量作为特征
X <- dat_SOC_MBC %>% select(soil.pH, MAP, MAT, elevation, bulk.density, clay, silt, sand,
                   crop.age, climatic.zone, PFT.peren, soc.peren.mean_30cm)

# 目标变量
y <- dat_SOC_MBC$SOC.effect.size  # 如果 'effect.size' 是目标变量

# 处理分类变量，确保它们是因子类型
X$climatic.zone <- as.factor(X$climatic.zone)
X$PFT.peren <- as.factor(X$PFT.peren)

# 分割数据为训练集和测试集
set.seed(42)  # 设置随机种子以确保结果可重复
trainIndex <- createDataPartition(y, p = 0.8, list = FALSE)  # 80%为训练集，20%为测试集
X_train <- X[trainIndex, ]
y_train <- y[trainIndex]
X_test <- X[-trainIndex, ]
y_test <- y[-trainIndex]

# 训练随机森林模型
rf_model <- randomForest(x = X_train, y = y_train, ntree = 100, importance = TRUE)

# 预测并计算误差
y_pred <- predict(rf_model, X_test)
mse <- mean((y_test - y_pred)^2)
r2 <- cor(y_test, y_pred)^2

# 输出评估结果
cat("Mean Squared Error (MSE):", mse, "\n")
cat("R²:", r2, "\n")

# 获取特征重要性
importance(rf_model)

# 输出特征重要性
varImpPlot(rf_model)

rf_importance <- read.csv("SOC_Random.csv", header = TRUE)

# 清理数据
names(rf_importance) <- c("Variable", "Importance")  # 修改列名
rf_importance$Variable <- factor(rf_importance$Variable, levels = rf_importance$Variable[order(rf_importance$Importance, decreasing = F)])
rf_importance <- rf_importance[order(rf_importance$Importance, decreasing = TRUE), ]
rf_importance$ColorGroup <- ifelse(1:nrow(rf_importance) <= 3, "Top3", "Others")

# 绘制条形图
ggplot(rf_importance, aes(x = Variable, y = Importance, fill = ColorGroup)) +
  geom_bar(stat = "identity", show.legend = FALSE, position = 'dodge', width = 0.7) +   # 绘制条形图
  scale_fill_manual(values = c("Top3" = "#c16528", "Others" = "#e2a06b"))+
  coord_flip() +                                      # 横向显示条形图
  labs(title = "Random Forest with LnRR-SOC", x = "", y = "Percentage Increase in Mean Squared Error") +
  scale_y_continuous(limits = c(0, 9),breaks = seq(0, 10, 2),expand = c(0,0))+
  theme(
    panel.background = element_blank(),  # 去掉背景色
    plot.background = element_blank(),   # 去掉图表的背景
    axis.text.x = element_text(angle = 0, colour='black', hjust = 1, size = 14),  # x轴标签不歪
    axis.text.y = element_text(size = 14, colour='black',),
    axis.title = element_text(size = 14, face = "bold"),  # 设置轴标题的字体
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),  # 设置标题
    panel.grid.major = element_blank(),  # 去掉中间的网格
    panel.grid.minor = element_blank(),  # 去掉次级网格
    axis.line = element_line(color = "black", size = 1),
    axis.ticks.length = unit(0.1, "cm"),
    plot.margin = unit(c(10, 10, 10, 10), "pt")  # 调整图边距
  )

###########MBC
# 选择感兴趣的特征变量（连续变量和分类变量）
# 筛选连续变量和分类变量作为特征
X <- dat_SOC_MBC %>% select(soil.pH, MAP, MAT, elevation, bulk.density, clay, silt, sand,
                            crop.age, climatic.zone, PFT.peren, soc.peren.mean_30cm)

# 目标变量
y <- dat_SOC_MBC$MBC.effect.size  # 如果 'effect.size' 是目标变量

# 处理分类变量，确保它们是因子类型
X$climatic.zone <- as.factor(X$climatic.zone)
X$PFT.peren <- as.factor(X$PFT.peren)

# 分割数据为训练集和测试集
set.seed(42)  # 设置随机种子以确保结果可重复
trainIndex <- createDataPartition(y, p = 0.8, list = FALSE)  # 80%为训练集，20%为测试集
X_train <- X[trainIndex, ]
y_train <- y[trainIndex]
X_test <- X[-trainIndex, ]
y_test <- y[-trainIndex]

# 训练随机森林模型
rf_model <- randomForest(x = X_train, y = y_train, ntree = 100, importance = TRUE)

# 预测并计算误差
y_pred <- predict(rf_model, X_test)
mse <- mean((y_test - y_pred)^2)
r2 <- cor(y_test, y_pred)^2

# 输出评估结果
cat("Mean Squared Error (MSE):", mse, "\n")
cat("R²:", r2, "\n")

# 获取特征重要性
importance(rf_model)

# 输出特征重要性
varImpPlot(rf_model)

rf_importance_MBC <- read.csv("MBC_Random.csv", header = TRUE)

# 清理数据
names(rf_importance_MBC) <- c("Variable", "Importance")  # 修改列名
rf_importance_MBC$Variable <- factor(rf_importance_MBC$Variable, levels = rf_importance_MBC$Variable[order(rf_importance_MBC$Importance, decreasing = F)])
rf_importance_MBC <- rf_importance_MBC[order(rf_importance_MBC$Importance, decreasing = TRUE), ]
rf_importance_MBC$ColorGroup <- ifelse(1:nrow(rf_importance_MBC) <= 3, "Top3", "Others")

# 绘制条形图
ggplot(rf_importance_MBC, aes(x = Variable, y = Importance, fill = ColorGroup)) +
  geom_bar(stat = "identity", show.legend = FALSE, position = 'dodge', width = 0.7) +   # 绘制条形图
  scale_fill_manual(values = c("Top3" = "#145742", "Others" = "#54bc57"))+
  coord_flip() +                                      # 横向显示条形图
  labs(title = "Random Forest with LnRR-MBC", x = "", y = "Percentage Increase in Mean Squared Error") +
  scale_y_continuous(limits = c(0, 9),breaks = seq(0, 10, 2),expand = c(0,0))+
  theme(
    panel.background = element_blank(),  # 去掉背景色
    plot.background = element_blank(),   # 去掉图表的背景
    axis.text.x = element_text(angle = 0, colour='black', hjust = 1, size = 14),  # x轴标签不歪
    axis.text.y = element_text(size = 14, colour='black',),
    axis.title = element_text(size = 14, face = "bold"),  # 设置轴标题的字体
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),  # 设置标题
    panel.grid.major = element_blank(),  # 去掉中间的网格
    panel.grid.minor = element_blank(),  # 去掉次级网格
    axis.line = element_line(color = "black", size = 1),
    axis.ticks.length = unit(0.1, "cm"),
    plot.margin = unit(c(10, 10, 10, 10), "pt")  # 调整图边距
  )



# 假设您的数据框叫 dat
# 加载所需包
install.packages(c("ggplot2", "reshape2", "Hmisc", "RColorBrewer"))
library(ggplot2)
library(reshape2)
library(Hmisc)
library(RColorBrewer)
library(grid)

# 选取变量（请替换为你的数据框 dat）
vars <- c("SOC.effect.size", "MBC.effect.size", "soil.pH", "MAP", "MAT", "elevation",
          "bulk.density", "clay", "silt", "sand", "soc.peren.mean_30cm")
data_sub <- dat_SOC_MBC[, vars]


# 计算 Pearson 相关性和 p 值
rc <- rcorr(as.matrix(data_sub), type = "pearson")
cor_mat <- rc$r
p_mat <- rc$P

# 提取 SOC/MBC 响应比与环境变量的相关性与 p 值
cor_df <- melt(cor_mat[c("SOC.effect.size", "MBC.effect.size"), -which(colnames(cor_mat) %in% c("SOC.effect.size", "MBC.effect.size"))])
p_df   <- melt(p_mat[c("SOC.effect.size", "MBC.effect.size"), -which(colnames(p_mat) %in% c("SOC.effect.size", "MBC.effect.size"))])

# 整理数据
colnames(cor_df) <- c("EffectType", "Variable", "Correlation")
colnames(p_df) <- c("EffectType", "Variable", "pvalue")
df_plot <- merge(cor_df, p_df)

# 显著性符号与点大小
df_plot$signif <- cut(df_plot$pvalue,
                      breaks = c(-Inf, 0.001, 0.01, 0.05, 1),
                      labels = c("***", "**", "*", ""))

df_plot$size <- cut(df_plot$pvalue,
                    breaks = c(-Inf, 0.001, 0.01, 0.05, 1),
                    labels = c(6, 4.5, 3.5, 2))

df_plot$size <- as.numeric(as.character(df_plot$size))
custom_palette <- scale_fill_gradient2(
  low = "#2166ac", mid = "white", high = "#b2182b",
  midpoint = 0, limits = c(-1, 1),
  name = "Pearson r"
)

df_plot$EffectType <- factor(df_plot$EffectType, levels = c("SOC.effect.size", "MBC.effect.size"), labels = c(
  "LnRR-SOC", "LnRR-MBC"))

df_plot$Variable <- factor(df_plot$Variable, levels = c(
  "soil.pH", "MAP", "MAT", "elevation", "bulk.density",
  "clay", "silt", "sand", "soc.peren.mean_30cm"
), labels = c(
  "Soil pH", "MAP", "MAT", "Elevation", "Bulk Density",
  "Clay(%)", "Silt(%)", "Sand(%)", "SOC(initial)"
))

# gplot 绘图（修改版：x = Variable, y = EffectType）
hot_MBC_SOC <- ggplot(df_plot, aes(x = Variable, y = EffectType)) +
  geom_tile(aes(fill = Correlation), color = "grey70", linewidth = 0.5, width = 0.9, height = 0.9) +  # 方块 + 边框
  geom_text(aes(label = signif), size = 5, color = "black", fontface = "bold") +                 # 显著性标注
  scale_fill_gradient2(
    low = "#0571b0", mid = "white", high = "#ca0020",
    midpoint = 0, limits = c(-0.5, 0.5),  # 可视需求调节范围
    name = "Pearson r"
  ) +
  scale_size(range = c(3, 8), guide = "none") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1, face = "bold", color = "black", size = 12),
    axis.text.y = element_text(face = "bold", color = "black", size = 12),
    axis.title = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 14),
    panel.grid = element_blank(),  # 不需要额外网格线
    plot.margin = unit(c(10, 10, 10, 10), "pt")
  ) +
  labs(title = "Correlations with Environmental Variables and Log Response Ratio")
print(hot_MBC_SOC)





