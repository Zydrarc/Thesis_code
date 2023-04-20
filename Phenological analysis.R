library(vegan)
library("dplyr")
library(viridis)
library(ggplot2)
library(pracma)
library(nls.multstart)
library(nlme)
library(pwr)
library(drc)


options(scipen=0) # scientific notation on
options(scipen=999) # scientific notation off (this is needed as the loop would not be calling the right file with out this)

experiments = 5 # number of experiments run
replicates = 10
tstep= 100000
tsave= 5000 # every how many timesteps the experiment is saved 
timesteps = seq(0,tstep, by = tsave)
tsnum = ((tstep/tsave)+1)
start_res = 1 # how much resources does the experiment start out with
res_num = c(2^seq((start_res-1),(experiments-1),by=1))
#res_num = c(2,4,8,16,32)
#res_num = c(1,2,3,4,5)
wd<-getwd()

hue <- viridis(experiments, alpha = 1, begin = 0, end = 1, direction = 1,option = "viridis")

#Phen data

phendata <- data.frame()
auc_list <- data.frame()

for (rep in 1:replicates) {
  for (ex in 1:experiments){
    #read data from the results of experiment
    pc <- paste(wd, "/replicate_", rep, "/experiment_", ex, "/data/phenotype_count.dat", sep = "")
    phen <- read.table(file = pc, skip = 8, sep = " ", dec = ".", header = FALSE, skipNul = TRUE, fill = TRUE )
    
    auc_shan <- trapz(phen$V1, phen$V3)
    auc_phen <- trapz(phen$V1, phen$V4)
    auc_aD <- trapz(phen$V1, phen$V6)
    newauc <- c(ex, auc_shan, auc_phen, auc_aD)
    auc_list <- rbind(auc_list, newauc)
    
    for (i in 1:tsnum){
      newphen <- c(ex, phen$V1[i], phen$V3[i], phen$V4[i], phen$V6[i])
      phendata <- rbind(phendata, newphen)
    }
  }
}

colnames(phendata) <- c("Exp", "ts", "Shan", "phen", "aD")
colnames(auc_list) <- c("Exp", "auc_Shan", "auc_phen", "auc_aD")

colnames(phendata)=c("Exp","ts","Shan","phen","aD")
colnames(auc_list)=c("Exp","auc_Shan","auc_phen","auc_aD")
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
ggplot(data = phen_combined, aes(x = ts, y = mean_Shan, color = factor(Exp))) +
  geom_line(size = 1) +
  geom_errorbar(aes(ymin = mean_Shan - std_Shan, ymax = mean_Shan + std_Shan), width = 0.2) +
  scale_color_manual(values = hue,labels = res_num) +
  labs(x = "Timesteps", y = "Shannon's Index") +
  ylim(0, 4.5) +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 12)) +  # increase font size of legend title
  labs(color = "Number of initial resources") +
  guides(color = guide_legend(title = "Number of initial resources", title.position = "top", label.position = "bottom", label.hjust = 0.5)) +
  theme(legend.text=element_text(size=12)) 
  

# plot Unique Phenotypes with error bars
ggplot(data = phen_combined, aes(x = ts, y = mean_phen, color = factor(Exp))) +
  geom_line(size = 1) +
  geom_errorbar(aes(ymin = mean_phen - std_phen, ymax = mean_phen + std_phen), width = 0.2) +
  scale_color_manual(values = hue,labels = res_num) +
  labs(x = "Timesteps", y = "Unique Phenotypes") +
  ylim(0, 1000) +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 12)) +  # increase font size of legend title
  labs(color = "Number of initial resources") +
  guides(color = guide_legend(title = "Number of initial resources", title.position = "top", label.position = "bottom", label.hjust = 0.5)) +
  theme(legend.text=element_text(size=12))   

# plot Average Task Diversity with error bars
ggplot(data = phen_combined, aes(x = ts, y = mean_aD, color = factor(Exp))) +
  geom_line(size = 1) +
  geom_errorbar(aes(ymin = mean_aD - std_aD, ymax = mean_aD + std_aD), width = 0.2) +
  scale_color_manual(values = hue,labels = res_num) +
  labs(x = "Timesteps", y = "Average Task Diversity") +
  ylim(0, 4) +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 12)) +  # increase font size of legend title
  labs(color = "Number of initial resources") +
  guides(color = guide_legend(title = "Number of initial resources", title.position = "top", label.position = "bottom", label.hjust = 0.5)) +
  theme(legend.text=element_text(size=12))   


##AUC Plots

ggplot(data = auc_list, aes(x = factor(Exp), y = auc_Shan, fill = factor(Exp))) +
  geom_boxplot() +
  scale_fill_manual(values = hue, labels = res_num) +
  labs(x = "Experiment", y = "Shannon's Index AUC") +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 12)) +  # increase font size of legend title
  guides(fill = guide_legend(title = "Number of initial resources", title.position = "top", label.position = "bottom", label.hjust = 0.5)) +
  theme(legend.text = element_text(size = 12))

###
final_species <- phendata %>%
  group_by(Exp) %>%
  filter(ts == max(ts))%>%
summarize(mean_Shan = mean(Shan),
          std_Shan = sd(Shan),
          mean_phen = mean(phen),
          std_phen = sd(phen),
          mean_aD = mean(aD),
          std_aD = sd(aD))

# compute the average AUC for each experiment and number of initial resources
auc_by_res <- auc_list %>%
  group_by(Exp) %>%
  mutate(Res = res_num[Exp]) %>%
  group_by(Res, add = TRUE) %>%
  summarize(avg_auc_phen = mean(auc_phen))

# plot Resource Diversity vs. AUC
ggplot(data = auc_by_res, aes(x = Res, y = avg_auc_phen)) +
  geom_point(size = 3, aes(color = factor(Exp))) +
  scale_color_manual(values = hue,labels = res_num) +
  labs(x = "Resource Diversity", y = "Phenotype AUC") +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 12)) +  # increase font size of legend title
  labs(color = "Experiment") +
  guides(color = guide_legend(title = "Experiment", title.position = "top", label.position = "bottom", label.hjust = 0.5)) +
  theme(legend.text=element_text(size=12))

# compute the final species richness for each experiment and number of initial resources
sp_richness <- phen_combined %>%
  group_by(Exp) %>%
  mutate(Res = res_num[Exp]) %>%
  slice_tail(n = 1)


# plot Resource Diversity vs. final species richness
ggplot(data = sp_richness, aes(x = Res, y = mean_phen)) +
  geom_point(size = 3, aes(color = factor(Exp))) +
  scale_color_manual(values = hue,labels = res_num) +
  labs(x = "Resource Diversity", y = "Final Species Richness") +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 12),
        legend.text=element_text(size=12)) +
  labs(color = "Experiment") +
  guides(color = guide_legend(title = "Experiment", title.position = "top", label.position = "bottom", label.hjust = 0.5)) 


##logistic curve
model <- drm(sp_richness$mean_Shan~sp_richness$Res, fct = L.3(), data = sp_richness)
summary(model)
plot(model, log="", main = "Logistic function", xlab = "Number of resources", ylab = "Mean Shannon's")



