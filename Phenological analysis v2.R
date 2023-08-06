library(vegan)
library(dplyr)
library(viridis)
library(ggplot2)
library(pracma)
library(nls.multstart)
library(nlme)
library(pwr)
#library(drc)

options(scipen=0)
options(scipen=999)

exptype = 1
experiments = 5
replicates = 10
tstep = 100000
tsave = 5000
timesteps = seq(0, tstep, by = tsave)
tsnum = ((tstep/tsave) + 1)
start_res = 1
res_num = c(2^seq((start_res-1), (experiments-1), by=1))
trophic_connections = c(4, 7, 10, 14, 18)
wd <- getwd()

hue <- viridis(experiments, alpha = 1, begin = 0, end = 1, direction = 1, option = "viridis")

if (exptype == 1) {
  # Phen data
  phendata <- data.frame()
  auc_list <- data.frame()
  
  for (rep in 1:replicates) {
    for (ex in 1:experiments) {
      # read data from the results of the experiment
      pc <- paste(wd, "/replicate_", rep, "/experiment_", ex, "/data/phenotype_count.dat", sep = "")
      phen <- read.table(file = pc, skip = 8, sep = " ", dec = ".", header = FALSE, skipNul = TRUE, fill = TRUE)
      
      auc_shan <- trapz(phen$V1, phen$V3)
      auc_phen <- trapz(phen$V1, phen$V4)
      auc_aD <- trapz(phen$V1, phen$V6)
      newauc <- c(ex, auc_shan, auc_phen, auc_aD)
      auc_list <- rbind(auc_list, newauc)
      
      for (i in 1:tsnum) {
        newphen <- c(ex, phen$V1[i], phen$V3[i], phen$V2[i], phen$V6[i])
        phendata <- rbind(phendata, newphen)
      }
    }
  }
  
  colnames(phendata) <- c("Exp", "ts", "Shan", "phen", "aD")
  colnames(auc_list) <- c("Exp", "auc_Shan", "auc_phen", "auc_aD")
  
  # Combining data from replicates and computing the average and standard deviation
  phen_combined <- phendata %>%
    group_by(Exp, ts) %>%
    summarize(mean_Shan = mean(Shan),
              std_Shan = sd(Shan),
              mean_phen = mean(phen),
              std_phen = sd(phen),
              mean_aD = mean(aD),
              std_aD = sd(aD))
  
  auc_combined <- auc_list %>%
    group_by(Exp) %>%
    summarize(mean_phenAUC = mean(auc_phen),
              std_phenAUC = sd(auc_phen))
  
  # plot Shannon's Index with error bars
  p1 <- ggplot(data = phen_combined, aes(x = ts, y = mean_Shan, color = factor(Exp))) +
    geom_line(size = 1) +
    geom_errorbar(aes(ymin = mean_Shan - std_Shan, ymax = mean_Shan + std_Shan), width = 0.2) +
    scale_color_manual(values = hue, labels = res_num) +
    labs(x = "Timesteps", y = "Shannon's Index") +
    ylim(0, 4.5) +
    theme_classic() +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 12)
    ) +
    labs(color = "Number of initial resources") +
    guides(
      color = guide_legend(
        title = "Number of initial resources",
        title.position = "top",
        label.position = "bottom",
        label.hjust = 0.5
      )
    )
  print(p1)
  # plot Unique Phenotypes with error bars
  p2 <- ggplot(data = phen_combined, aes(x = ts, y = mean_phen, color = factor(Exp))) +
    geom_line(size = 1) +
    geom_errorbar(aes(ymin = mean_phen - std_phen, ymax = mean_phen + std_phen), width = 0.2) +
    scale_color_manual(values = hue, labels = res_num) +
    labs(x = "Timesteps", y = "Unique Phenotypes") +
    ylim(0, 75) +
    theme_classic() +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 12)
    ) +
    labs(color = "Number of initial resources") +
    guides(
      color = guide_legend(
        title = "Number of initial resources",
        title.position = "top",
        label.position = "bottom",
        label.hjust = 0.5
      )
    )
  print(p2)
  # plot Average Task Diversity with error bars
  p3 <- ggplot(data = phen_combined, aes(x = ts, y = mean_aD, color = factor(Exp))) +
    geom_line(size = 1) +
    geom_errorbar(aes(ymin = mean_aD - std_aD, ymax = mean_aD + std_aD), width = 0.2) +
    scale_color_manual(values = hue, labels = res_num) +
    labs(x = "Timesteps", y = "Average Task Diversity") +
    ylim(0, 4) +
    theme_classic() +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 12)
    ) +
    labs(color = "Number of initial resources") +
    guides(
      color = guide_legend(
        title = "Number of initial resources",
        title.position = "top",
        label.position = "bottom",
        label.hjust = 0.5
      )
    )
  print(p3)
  # AUC Plots
  p4 <- ggplot(data = auc_list, aes(x = factor(Exp), y = auc_Shan, fill = factor(Exp))) +
    geom_boxplot() +
    scale_fill_manual(values = hue, labels = res_num) +
    labs(x = "Experiment", y = "Shannon's Index AUC") +
    theme_classic() +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 12)
    ) +
    guides(
      fill = guide_legend(
        title = "Food Webs",
        title.position = "top",
        label.position = "bottom",
        label.hjust = 0.5
      )
    )
  print(p4)
  # Final species richness
  final_species <- phendata %>%
    group_by(Exp) %>%
    filter(ts == max(ts)) %>%
    summarize(
      mean_Shan = mean(Shan),
      std_Shan = sd(Shan),
      mean_phen = mean(phen),
      std_phen = sd(phen),
      mean_aD = mean(aD),
      std_aD = sd(aD)
    )
  
  final_species$resources <- res_num
  
  # Compute the average AUC for each experiment and number of initial resources
  auc_by_res <- auc_list %>%
    group_by(Exp) %>%
    mutate(Res = res_num[Exp]) %>%
    group_by(Res, add = TRUE) %>%
    summarize(avg_auc_phen = mean(auc_phen))
  
  # Plot Resource Diversity vs. AUC
  p5 <- ggplot(data = auc_by_res, aes(x = Res, y = avg_auc_phen)) +
    geom_point(size = 3, aes(color = factor(Exp))) +
    scale_color_manual(values = hue, labels = res_num) +
    labs(x = "Resource Diversity", y = "Phenotype AUC") +
    theme_classic() +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 12)
    ) +
    labs(color = "Experiment") +
    guides(
      color = guide_legend(
        title = "Experiment",
        title.position = "top",
        label.position = "bottom",
        label.hjust = 0.5
      )
    )
  print(p5)
  # Compute the final species richness for each experiment and number of initial resources
  sp_richness <- phen_combined %>%
    group_by(Exp) %>%
    mutate(Res = res_num[Exp]) %>%
    slice_tail(n = 1)
  
  # Plot Resource Diversity vs. final species richness
  p6 <- ggplot(data = sp_richness, aes(x = Res, y = mean_phen)) +
    geom_point(size = 3, aes(color = factor(Exp))) +
    scale_color_manual(values = hue, labels = res_num) +
    labs(x = "Resource Diversity", y = "Final Species Richness") +
    theme_classic() +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 12)
    ) +
    labs(color = "Experiment") +
    guides(
      color = guide_legend(
        title = "Experiment",
        title.position = "top",
        label.position = "bottom",
        label.hjust = 0.5
      )
    )
  print(p6)
  
  # Fit a GLM model

  glm_model <- glm(sp_richness$mean_Shan ~ sp_richness$Res, data = sp_richness)
  summary(glm_model)
  
  residuals <- resid(glm_model)
  fitted.values <- fitted(glm_model)
  plot(fitted.values, residuals)
  
  qqnorm(residuals)
  qqline(residuals)
  
  cooks.distance <- cooks.distance(glm_model)
  plot(cooks.distance, main="Cook's Distance")
  abline(h = 1, col="red")
  

  
  # Logistic curve
  model <- drm(sp_richness$mean_Shan ~ sp_richness$Res, fct = L.3(), data = sp_richness)
  summary(model)
  p7 <- plot(model, log = "", main = "Logistic function", xlab = "Number of resources", ylab = "Mean Shannon's")
  print(p7)
  residuals <- resid(model)
  plot(predict(model), residuals)
  abline(h = 0, lty = 2)
  
  qqnorm(residuals)
  qqline(residuals)
  
  SSE <- sum(residuals^2)
  SST <- sum((sp_richness$mean_Shan - mean(sp_richness$mean_Shan))^2)
  pseudoR2 <- 1 - SSE/SST
  pseudoR2
  
  # Compute Shannon's index for each replicate and experiment
  shan_index <- phendata %>%
    group_by(Exp, ts) %>%
    summarize(mean_Shan = mean(Shan), std_Shan = sd(Shan))
  
  # Perform ANOVA to compare Shannon's index between experiments
  anova_shan <- aov(mean_Shan ~ factor(Exp), data = shan_index)
  # Extract ANOVA results
  anova_summary_shan <- summary(anova_shan)
  # Print ANOVA results
  cat("ANOVA Results for Shannon's Index:\n")
  cat("-----------------------------------------------------\n")
  print(anova_summary_shan)
}

if (exptype == 2) {
  tstep = 200000
  # Phen data
  phendata <- data.frame()
  auc_list <- data.frame()
  
  for (rep in 1:replicates) {
    for (ex in 1:experiments) {
      # read data from the results of the experiment
      pc <- paste(wd, "/replicate_", rep, "/experiment_", ex, "/data/phenotype_count.dat", sep = "")
      phen <- read.table(file = pc, skip = 8, sep = " ", dec = ".", header = FALSE, skipNul = TRUE, fill = TRUE)
      
      auc_shan <- trapz(phen$V1, phen$V3)
      auc_phen <- trapz(phen$V1, phen$V4)
      auc_aD <- trapz(phen$V1, phen$V6)
      newauc <- c(ex, auc_shan, auc_phen, auc_aD)
      auc_list <- rbind(auc_list, newauc)
      
      for (i in 1:tsnum) {
        newphen <- c(ex, phen$V1[i], phen$V3[i], phen$V2[i], phen$V6[i])
        phendata <- rbind(phendata, newphen)
      }
    }
  }
  
  colnames(phendata) <- c("Exp", "ts", "Shan", "phen", "aD")
  colnames(auc_list) <- c("Exp", "auc_Shan", "auc_phen", "auc_aD")
  
  # Combining data from replicates and computing the average and standard deviation
  phen_combined <- phendata %>%
    group_by(Exp, ts) %>%
    summarize(mean_Shan = mean(Shan),
              std_Shan = sd(Shan),
              mean_phen = mean(phen),
              std_phen = sd(phen),
              mean_aD = mean(aD),
              std_aD = sd(aD))
  
  auc_combined <- auc_list %>%
    group_by(Exp) %>%
    summarize(mean_phenAUC = mean(auc_phen),
              std_phenAUC = sd(auc_phen))
  
  # plot Shannon's Index with error bars
  p1 <- ggplot(data = phen_combined, aes(x = ts, y = mean_Shan, color = factor(Exp))) +
    geom_line(size = 1) +
    geom_errorbar(aes(ymin = mean_Shan - std_Shan, ymax = mean_Shan + std_Shan), width = 0.2) +
    scale_color_manual(values = hue, labels = trophic_connections) +
    labs(x = "Timesteps", y = "Shannon's Index") +
    ylim(0, 4.5) +
    theme_classic() +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 12)
    ) +
    labs(color = "Number of trophic connection") +
    guides(
      color = guide_legend(
        title = "Number of trophic connection",
        title.position = "top",
        label.position = "bottom",
        label.hjust = 0.5
      )
    )
  print(p1)
  # plot Unique Phenotypes with error bars
  p2 <- ggplot(data = phen_combined, aes(x = ts, y = mean_phen, color = factor(Exp))) +
    geom_line(size = 1) +
    geom_errorbar(aes(ymin = mean_phen - std_phen, ymax = mean_phen + std_phen), width = 0.2) +
    scale_color_manual(values = hue, labels = trophic_connections) +
    labs(x = "Timesteps", y = "Unique Phenotypes") +
    ylim(0, 75) +
    theme_classic() +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 12)
    ) +
    labs(color = "Number of trophic connection") +
    guides(
      color = guide_legend(
        title = "Number of trophic connection",
        title.position = "top",
        label.position = "bottom",
        label.hjust = 0.5
      )
    )
  print(p2)
  # plot Average Task Diversity with error bars
  p3 <- ggplot(data = phen_combined, aes(x = ts, y = mean_aD, color = factor(Exp))) +
    geom_line(size = 1) +
    geom_errorbar(aes(ymin = mean_aD - std_aD, ymax = mean_aD + std_aD), width = 0.2) +
    scale_color_manual(values = hue, labels = trophic_connections) +
    labs(x = "Timesteps", y = "Average Task Diversity") +
    ylim(0, 4) +
    theme_classic() +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 12)
    ) +
    labs(color = "Number of trophic connection") +
    guides(
      color = guide_legend(
        title = "Number of trophic connection",
        title.position = "top",
        label.position = "bottom",
        label.hjust = 0.5
      )
    )
  print(p3)
  # AUC Plots
  p4 <- ggplot(data = auc_list, aes(x = factor(Exp), y = auc_Shan, fill = factor(Exp))) +
    geom_boxplot() +
    scale_fill_manual(values = hue, labels = trophic_connections) +
    labs(x = "Experiment", y = "Shannon's Index AUC") +
    theme_classic() +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 12)
    ) +
    guides(
      fill = guide_legend(
        title = "Number of trophic connection",
        title.position = "top",
        label.position = "bottom",
        label.hjust = 0.5
      )
    )
  print(p4)
  # Final species richness
  final_species <- phendata %>%
    group_by(Exp) %>%
    filter(ts == max(ts)) %>%
    summarize(
      mean_Shan = mean(Shan),
      std_Shan = sd(Shan),
      mean_phen = mean(phen),
      std_phen = sd(phen),
      mean_aD = mean(aD),
      std_aD = sd(aD)
    )
  
  final_species$resources <- res_num
  
  # Compute the average AUC for each experiment and number of initial resources
  auc_by_res <- auc_list %>%
    group_by(Exp) %>%
    mutate(trophic = trophic_connections[Exp]) %>%
    group_by(trophic, add = TRUE) %>%
    summarize(avg_auc_phen = mean(auc_phen))
  
  # Plot Resource Diversity vs. AUC
  p5 <- ggplot(data = auc_by_res, aes(x = trophic, y = avg_auc_phen)) +
    geom_point(size = 3, aes(color = factor(Exp))) +
    scale_color_manual(values = hue, labels = trophic_connections) +
    labs(x = "Trophic diversity", y = "Phenotype AUC") +
    theme_classic() +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 12)
    ) +
    labs(color = "Experiment") +
    guides(
      color = guide_legend(
        title = "Experiment",
        title.position = "top",
        label.position = "bottom",
        label.hjust = 0.5
      )
    )
  print(p5)
  # Compute the final species richness for each experiment and number of initial resources
  sp_richness <- phen_combined %>%
    group_by(Exp) %>%
    mutate(trophic = trophic_connections[Exp]) %>%
    slice_tail(n = 1)
  
  # Plot Resource Diversity vs. final species richness
  p6 <- ggplot(data = sp_richness, aes(x = trophic, y = mean_phen)) +
    geom_point(size = 3, aes(color = factor(Exp))) +
    scale_color_manual(values = hue, labels = trophic_connections) +
    labs(x = "Trophic Diversity", y = "Final Species Richness") +
    theme_classic() +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 12)
    ) +
    labs(color = "Experiment") +
    guides(
      color = guide_legend(
        title = "Trophic connections",
        title.position = "top",
        label.position = "bottom",
        label.hjust = 0.5
      )
    )
  print(p6)  
  

  # Logistic curve
  model <- drm(sp_richness$mean_Shan ~ sp_richness$trophic, fct = L.3(), data = sp_richness)
  summary(model)
  p7 <- plot(model, log = "", main = "Logistic function", xlab = "Number of resources", ylab = "Mean Shannon's")
  residuals <- resid(model)
  plot(predict(model), residuals)
  abline(h = 0, lty = 2)
  
  qqnorm(residuals)
  qqline(residuals)
  
  SSE <- sum(residuals^2)
  SST <- sum((sp_richness$mean_Shan - mean(sp_richness$mean_Shan))^2)
  pseudoR2 <- 1 - SSE/SST
  pseudoR2
  

  
  
  # Compute Shannon's index for each replicate and experiment
  shan_index <- phendata %>%
    group_by(Exp, ts) %>%
    summarize(mean_Shan = mean(Shan), std_Shan = sd(Shan))
  
  # Perform ANOVA to compare Shannon's index between experiments
  anova_shan <- aov(mean_Shan ~ factor(Exp), data = shan_index)
  # Extract ANOVA results
  anova_summary_shan <- summary(anova_shan)
  # Print ANOVA results
  cat("ANOVA Results for Shannon's Index:\n")
  cat("-----------------------------------------------------\n")
  print(anova_summary_shan)
  
}
