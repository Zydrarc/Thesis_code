library(vegan)
library("dplyr")
library(viridis)
library(ggplot2)
library(pracma)
library(nls.multstart)
library(nlme)
library(pwr)
library(drc)
library(tidyr)
library(stringr)

options(scipen=0) # scientific notation on
options(scipen=999) # scientific notation off

experiments = 5 # number of experiments run
replicates = 20
tstep= 200000
tsave= 5000 # every how many timesteps the experiment is saved 
timesteps = seq(0,tstep, by = tsave)
tsnum = ((tstep/tsave)+1)
start_res = 1 # how much resources does the experiment start out with
res_num = c(2^seq((start_res-1),(experiments-1),by=1))

wd<-getwd()

hue <- viridis(experiments, alpha = 1, begin = 0, end = 1, direction = 1,option = "viridis")

trophic_data <- data.frame()

for (rep in 1:replicates) {
  for (ex in 1:experiments){
    pc <- paste(wd, "/replicate_", rep, "/experiment_", ex, "/data/tasks.dat", sep = "")
    troph <- read.table(file = pc, skip = 11, sep = " ", dec = ".", header = FALSE, skipNul = TRUE, fill = TRUE )
    
    if(ex == 1){
      for (i in 1:tsnum){
        newtroph <- c(ex, troph$V1[i],(troph$V2[i] +troph$V3[i]), troph$V4[i],troph$V5[i],troph$V6[i])
        trophic_data <- rbind(trophic_data, newtroph)
      }
    }
    if(ex == 2){
      for (i in 1:tsnum){
        newtroph <- c(ex, troph$V1[i],(troph$V2[i] +troph$V3[i]), (troph$V4[i]+troph$V5[i]),troph$V6[i],troph$V7[i])
        trophic_data <- rbind(trophic_data, newtroph)
      }
    }
    if(ex == 3){
      for (i in 1:tsnum){
        newtroph <- c(ex, troph$V1[i],(troph$V2[i] +troph$V3[i]), (troph$V4[i]+troph$V5[i]),(troph$V6[i]+troph$V7[i]),troph$V8[i])
        trophic_data <- rbind(trophic_data, newtroph)
      }
    }
    if(ex == 4){
      for (i in 1:tsnum){
        newtroph <- c(ex, troph$V1[i],(troph$V2[i] +troph$V3[i]), (troph$V4[i]+troph$V5[i]+troph$V6[i]),(troph$V7[i]+troph$V8[i]),troph$V9[i])
        trophic_data <- rbind(trophic_data, newtroph)
      }   
    }
    if(ex == 5){
      for (i in 1:tsnum){
        newtroph <- c(ex, troph$V1[i],(troph$V2[i] +troph$V3[i]), (troph$V4[i]+troph$V5[i]+troph$V6[i]),(troph$V7[i]+troph$V8[i]+troph$V9[i]),troph$V10[i])
        trophic_data <- rbind(trophic_data, newtroph)
      }   
    }
  }
}
colnames(trophic_data) <- c("Experiment", "Timestep", "Tl_1", "Tl_2", "Tl_3", "Tl_4")


# Calculate sum
grouped_trophic_data_sum <- trophic_data %>% 
  group_by(Experiment, Timestep) %>%
  summarise(Tl_1 = sum(Tl_1),
            Tl_2 = sum(Tl_2),
            Tl_3 = sum(Tl_3),
            Tl_4 = sum(Tl_4),
            .groups = "drop")

# Calculate mean and standard deviation
grouped_trophic_data_mean_sd <- trophic_data %>% 
  group_by(Experiment, Timestep) %>%
  summarise_each(funs(mean, sd), matches("^Tl_"))

print(grouped_trophic_data_sum)
print(grouped_trophic_data_mean_sd)

# percentage of each trophic level relative to the total
grouped_trophic_data_sum <- grouped_trophic_data_sum %>%
  mutate(Tl_1_pct = Tl_1 / (Tl_1 + Tl_2 + Tl_3 + Tl_4) * 100,
         Tl_2_pct = Tl_2 / (Tl_1 + Tl_2 + Tl_3 + Tl_4) * 100,
         Tl_3_pct = Tl_3 / (Tl_1 + Tl_2 + Tl_3 + Tl_4) * 100,
         Tl_4_pct = Tl_4 / (Tl_1 + Tl_2 + Tl_3 + Tl_4) * 100)

# ADD THESE LINES
# Convert wide data to long
# Convert wide data to long
long_data_pct <- grouped_trophic_data_sum %>% 
  gather(key = "Trophic_Level", value = "Percentage", Tl_1_pct:Tl_4_pct)

color_palette <- viridis(length(unique(long_data_pct$Trophic_Level)))
long_data_pct$Trophic_Level <- case_when(
  long_data_pct$Trophic_Level == "Tl_1_pct" ~ "1st trophic level",
  long_data_pct$Trophic_Level == "Tl_2_pct" ~ "2nd trophic level",
  long_data_pct$Trophic_Level == "Tl_3_pct" ~ "3rd trophic level",
  long_data_pct$Trophic_Level == "Tl_4_pct" ~ "4th trophic level"
)


# Plot
pct_plot <- ggplot(long_data_pct, aes(x = Timestep, y = Percentage, fill = Trophic_Level)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Experiment, scales = "free_x") +
  labs(x = "Timestep", y = "Relative Percentage (%)", fill = "Trophic Level", title = "Graph of relative percentage of trophic levels with increasing complexity") +
  scale_fill_manual(values = color_palette) +  # Set the color palette
  #scale_y_reverse() +  # Reverse the order of y-axis
  theme_minimal()
print(pct_plot)

# Plot for actual data
long_data_no_pct <- long_data %>%
  filter(!str_detect(Trophic_Level, "_pct")) %>%
  mutate(Percentage = Percentage * 100) # Convert to percentage

ggplot(long_data_no_pct, aes(x = Timestep, y = Percentage, fill = Trophic_Level)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Experiment) +
  labs(x = "Timestep", y = "Value", fill = "Trophic Level") +
  theme_minimal()

# Reshape data to long format for boxplot
trophic_data_long <- reshape2::melt(grouped_trophic_data_sum, id.vars = c("Experiment", "Timestep"),
                                    variable.name = "Trophic_Level", value.name = "Value")

# Filter the data for the trophic levels of interest
trophic_data_filtered <- trophic_data_long %>%
  filter(Trophic_Level %in% c("Tl_1_pct", "Tl_2_pct", "Tl_3_pct", "Tl_4_pct"))

# Update the Trophic_Level labels
trophic_data_filtered$Trophic_Level <- case_when(
  trophic_data_filtered$Trophic_Level == "Tl_1_pct" ~ "1st trophic level",
  trophic_data_filtered$Trophic_Level == "Tl_2_pct" ~ "2nd trophic level",
  trophic_data_filtered$Trophic_Level == "Tl_3_pct" ~ "3rd trophic level",
  trophic_data_filtered$Trophic_Level == "Tl_4_pct" ~ "4th trophic level",
  TRUE ~ trophic_data_filtered$Trophic_Level
)

# Create a boxplot with experiments on x-axis
ggplot(trophic_data_filtered, aes(x = factor(Experiment), y = Value, fill = Trophic_Level)) +
  geom_boxplot() +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "Experiment", y = "Trophic Prevalence (%)") +
  facet_wrap(~Trophic_Level, scales = "free")
