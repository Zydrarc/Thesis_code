library("dplyr")
library(stats)
library(ggplot2)

experiments = 5 # number of experiments run
replicates = 20
wd<-getwd()
## Evo of EQU

taskdata <- data.frame()


results <- data.frame(Experiment = integer(), `Presence of EQU` = integer()) 

for (rep in 1:replicates) {
  for (ex in 1:experiments){
    #read data from the results of experiment
    pw <- paste(wd, "/replicate_", rep, "/experiment_", ex, "/data/tasks.dat", sep = "")
    taskdata <- read.table(file = pw, skip = 13, sep = " ", dec = ".", header = FALSE, skipNul = TRUE, fill = TRUE )
    
    last_row <- nrow(taskdata)
    last_column <- (length(taskdata)-1)
    
    # Get the value in the last row and last column of the table
    val <- taskdata[last_row, last_column]
    if (!is.na(val)) {
      if (val > 0) {
        results[nrow(results) + 1,] <- list(ex, 1)
      } else {
        results[nrow(results) + 1,] <- list(ex, 0)
      }
    }
  }
}

# results is your new dataframe
presence_of_EQU <- results
presence_of_EQU <- presence_of_EQU %>%
  arrange(Experiment)

presence_of_EQU <- presence_of_EQU %>%
  mutate(Connections = case_when(
    Experiment == 1 ~ 4,
    Experiment == 2 ~ 7,
    Experiment == 3 ~ 10,
    Experiment == 4 ~ 14,
    Experiment == 5 ~ 18,
    TRUE ~ NA_real_
  ))

presence_of_EQU <- dplyr::select(presence_of_EQU, Experiment, Connections, `Presence.of.EQU`)


EQU_avg <- presence_of_EQU %>%
  group_by(Experiment, Connections) %>%
  summarise(Average = mean(`Presence.of.EQU`, na.rm = TRUE),
            SD = sd(`Presence.of.EQU`, na.rm = TRUE))

# Group data by Experiment and count number of times EQU emerged
counts_exp <- presence_of_EQU %>%
  group_by(Experiment) %>%
  summarise(count = sum(`Presence.of.EQU`, na.rm = TRUE))

# Fit the glm
model_exp <- glm(count ~ Experiment, data = counts_exp, family = poisson)

# Print summary of the model
summary(model_exp)

# Group data by Connections and count number of times EQU emerged
counts_con <- presence_of_EQU %>%
  group_by(Connections) %>%
  summarise(count = sum(`Presence.of.EQU`, na.rm = TRUE))

# Fit the glm
model_con <- glm(count ~ Connections, data = counts_con, family = poisson)

# Print summary of the model
summary(model_con)

model <- drm(count ~ Connections, fct = L.3(), data = counts_con)
summary(model)
p7 <- plot(model, log = "", main = "Logistic function", xlab = "Number of trophic connections", ylab = "Number of time EQU evolved")
print(p7)
residuals <- resid(model)
plot(predict(model), residuals)

abline(h = 0, lty = 2)



##model validation


##checking residuals
# Residual plots for model_con
par(mfrow = c(2, 1))

plot(predict(model_con, type = "response"), residuals(model_con, type = "deviance"),
     xlab = "Fitted values", ylab = "Deviance Residuals", main = "Residual plot for model_con")


plot(counts_con$Connections, residuals(model_con, type = "deviance"),
     xlab = "Connections", ylab = "Deviance Residuals", main = "Residual plot against Connections")

par(mfrow = c(1, 1))

##Check for overdispersion 

dispersion_con <- sum(residuals(model_con, type = "pearson")^2) / model_con$df.residual
print(paste("Dispersion parameter for model_con: ", dispersion_con))

# Check the goodness-of-fit with Chi-square test:

chisq_con <- sum(residuals(model_con, type = "pearson")^2)
p_value_con <- pchisq(chisq_con, df = model_con$df.residual, lower.tail = FALSE)
print(paste("Chi-square test p-value for model_con: ", p_value_con))

# Plotting the model
ggplot(counts_con, aes(x = Connections, y = count)) +
  geom_point() +
  geom_smooth(method = "glm", method.args = list(family = poisson), se = FALSE) +
  labs(x = "Connections", y = "Number of times EQU emerged") +
  ggtitle("Model: Number of times EQU emerged by Connections")


# Residual plots for model_exp
par(mfrow = c(2, 1))


plot(predict(model_exp, type = "response"), residuals(model_exp, type = "deviance"),
     xlab = "Fitted values", ylab = "Deviance Residuals", main = "Residual plot for model_exp")


plot(counts_exp$Experiment, residuals(model_exp, type = "deviance"),
     xlab = "Experiment", ylab = "Deviance Residuals", main = "Residual plot against Experiment")

par(mfrow = c(1, 1))


dispersion_exp <- sum(residuals(model_exp, type = "pearson")^2) / model_exp$df.residual
print(paste("Dispersion parameter for model_exp: ", dispersion_exp))

chisq_exp <- sum(residuals(model_exp, type = "pearson")^2)
p_value_exp <- pchisq(chisq_exp, df = model_exp$df.residual, lower.tail = FALSE)
print(paste("Chi-square test p-value for model_exp: ", p_value_exp))

# Plotting the model
ggplot(counts_exp, aes(x = Experiment, y = count)) +
  geom_point() +
  geom_smooth(method = "glm", method.args = list(family = poisson), se = FALSE) +
  labs(x = "Experiment", y = "Number of times EQU emerged") +
  ggtitle("Model: Number of times EQU emerged by Experiment")

