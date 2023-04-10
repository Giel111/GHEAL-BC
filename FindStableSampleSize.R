# Find Stable Sample Size

library(dplyr)
library(ggplot2)
library(stats)
library(scales)
setwd('C:/Users/gielv/OneDrive/UT/IEM/001 AFSTUDEREN/RProjects/ThesisModel') #<-home

hour <- 1
day=hour*24;
year = day*365;


# df <- read.csv('0_1_Output_DFs/1_20230317_3_Baseline1M_1080000EndDatDf.csv') #<- any enddat dataframe
df <- df %>%
  dplyr::select(TotalUtility) %>%
  mutate(TotalUtility = as.numeric(TotalUtility) / year)

df['TotalUtility'][is.na(df['TotalUtility'])] <- 786867.1/year

ShufVec <- df[sample(1:nrow(df)),]


# create a data frame with the vector and the index
df <- data.frame(ShufVec = ShufVec, Index = 1:length(ShufVec))

# calculate the running mean
df$Avg <- cumsum(df$ShufVec) / 1:nrow(df)


# calculate the confidence intervals
n <- nrow(df)
df$CI99 <- qt(0.995, n - 1) * sqrt(cumsum(df$ShufVec^2) / 1:n - df$Avg^2)
df$CI975 <- qt(0.9875, n - 1) * sqrt(cumsum(df$ShufVec^2) / 1:n - df$Avg^2)
df$CI95 <- qt(0.975, n - 1) * sqrt(cumsum(df$ShufVec^2) / 1:n - df$Avg^2)

# plot the running mean
ggplot(df, aes(x = Index)) +
  geom_line(aes(y = Avg), color = "black") +
  geom_line(aes(y = Avg + CI99), color = "red") +
  geom_line(aes(y = Avg - CI99), color = "red") +
  geom_line(aes(y = Avg + CI975), color = "orange") +
  geom_line(aes(y = Avg - CI975), color = "orange") +
  geom_line(aes(y = Avg + CI95), color = "yellow") +
  geom_line(aes(y = Avg - CI95), color = "yellow") +
  ylab("Running Mean")+
  labs(title = "Average of all patient utilities up to that index")

# Sample size
n <- length(ShufVec)

# Degrees of freedom
df <- n - 1

# Mean of the data
xbar <- mean(ShufVec)

# Standard deviation of the data
s <- sd(ShufVec)

# Calculate the standard error
se <- s / sqrt(n)

# t-value for 95% CI
t95 <- qt(0.975, df)

# t-value for 97.5% CI
t97.5 <- qt(0.9875, df)

# t-value for 99% CI
t99 <- qt(0.99995, df)

# Calculate the confidence intervals
ci95 <- c(xbar - t95 * se, xbar + t95 * se)
ci97.5 <- c(xbar - t97.5 * se, xbar + t97.5 * se)
ci99 <- c(xbar - t99 * se, xbar + t99 * se)

# Sample size
n <- round(length(ShufVec)/2)

# Degrees of freedom
df <- n - 1

# Mean of the data
xbar <- mean(ShufVec)

# Standard deviation of the data
s <- sd(ShufVec)

# Create an empty data frame
df_ci <- data.frame(index = integer(),
                    lower95 = double(),
                    upper95 = double(),
                    lower975 = double(),
                    upper975 = double(),
                    lower99 = double(),
                    upper99 = double())

# Loop through the data and calculate the confidence intervals
for (i in 1:n) {
  if (i%%5000==0){
    print(i)
  }
  # Subset the data
  x <- ShufVec[1:i]
  
  # Mean of the subset
  xbar_i <- mean(x)
  
  # Standard deviation of the subset
  s_i <- sd(x)
  
  # Calculate the standard error
  se_i <- s_i / sqrt(i)
  
  # t-value for 95% CI
  t95 <- qt(0.975, i - 1)
  
  # t-value for 97.5% CI
  t975 <- qt(0.9875, i - 1)
  
  # t-value for 99% CI
  t99 <- qt(0.995, i - 1)
  
  # Calculate the confidence intervals
  ci95 <- c(xbar - t95 * se_i, xbar + t95 * se_i)
  ci975 <- c(xbar - t975 * se_i, xbar + t975 * se_i)
  ci99 <- c(xbar - t99 * se_i, xbar + t99 * se_i)
  
  # Add the confidence intervals to the data frame
  df_ci[i, "index"] <- i
  df_ci[i, "lower95"] <- ci95[1]
  df_ci[i, "upper95"] <- ci95[2]
  df_ci[i, "lower975"] <- ci975[1]
  df_ci[i, "upper975"] <- ci975[2]
  df_ci[i, "lower99"] <- ci99[1]
  df_ci[i, "upper99"] <- ci99[2]
}

# create a data frame with the vector and the index
df <- data.frame(ShufVec = ShufVec, Index = 1:length(ShufVec))

# calculate the running mean
df$Avg <- cumsum(df$ShufVec) / 1:nrow(df)

# Add actual mean as well:
df_ci['actual_mean'] <- df$Avg[1:nrow(df_ci)]

df_ci <- read.csv('0_ExperimentResults/ConfidenceIntervalOfTotalQALY.csv')

# Plot the results
ggplot(df_ci, aes(x = index, y = actual_mean)) +
  geom_line(color = "black") +
  # ylim(c(71,74))+
  geom_ribbon(aes(ymin = lower95, ymax = upper95), fill = "Red", alpha = 0.8) +
  geom_ribbon(aes(ymin = lower975, ymax = upper975), fill = "Yellow", alpha = 0.5) +
  geom_ribbon(aes(ymin = lower99, ymax = upper99), fill = "Green", alpha = 0.1) +
  ggtitle(paste0("Running Mean and Confidence Intervals\n", "95%, 97.5%, and 99% for the average QALY")) +
  labs(x = "Number of patients",
       y = "Average QALY") +
  scale_color_manual(values = c("Red", "Yellow", "Green"), name = "Confidence Interval") +
  ylim(c(71.5,73))+
  theme(legend.position = c(0.85, 0.85))+
  theme(plot.title = element_text(hjust = 0.5))

