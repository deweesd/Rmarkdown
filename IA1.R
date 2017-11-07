setwd("/Users/Deweesd/Desktop/Advanced_Stats_with_R/IA1")
getwd


library(tidyverse)
library(magrittr)
library(pander)



Zebra_fish_diet_data <- read.table("zfish_diet_IA.tsv", sep = '\t', header = TRUE)

#View(Zebra_fish_diet_data)

# deleting outliers
Zebra_fish_diet_data_update <- Zebra_fish_diet_data[-27,]
Zebra_fish_diet_data_final <- Zebra_fish_diet_data_update[-140,]


#View(Zebra_fish_diet_data_final)



#Protein is independent variable

#Weight/SL is dependent variable 



#plot for entire dataset 
total_dataset_plot <- ggplot(data = Zebra_fish_diet_data_final, aes(y = Weight, x = SL, fill = Diet)) + geom_dotplot(binaxis = "x", stackdir = "center") + scale_fill_manual(values=c("slateblue1", "olivedrab1")) + theme(axis.text.x = element_text(angle = 40, vjust = 0.3, size = 6)) +
  labs(x=expression(bold("SL")), y=expression(bold("Weight")), title=expression(bold("Protein/Weight_Distribution")))



total_dataset_plot

#plot(weight_diet_plot)




# weight plot
Weight_Diet_plot <- ggplot(data = Zebra_fish_diet_data_final, aes(y = Weight, x = Diet)) + geom_boxplot(aes(fill = Diet)) + scale_fill_manual(values=c("slateblue1", "olivedrab1")) + theme(axis.text.x = element_text(angle = 40, vjust = 0.3, size = 6)) +
  labs(x=expression(bold("Diet")), y=expression(bold("Weight")), title=expression(bold("Diet to Weight Distribution")))

SL_Diet_plot <- ggplot(data = Zebra_fish_diet_data_final, aes(y = SL, x = Diet)) + geom_boxplot(aes(fill = Diet)) + scale_fill_manual(values=c("slateblue1", "olivedrab1")) + theme(axis.text.x = element_text(angle = 40, vjust = 0.3, size = 6)) +
  labs(x=expression(bold("Diet")), y=expression(bold("SL")), title=expression(bold("Diet to SL Distribution")))

require(gridExtra)
grid.arrange(Weight_Diet_plot, SL_Diet_plot, ncol=2)







#group needs to be only categorical variables

Data_summary_SL <- group_by(Zebra_fish_diet_data_final, Diet) %>%
  summarise_at("SL", .funs = c(Variance = "var", Mean="mean", SD="sd", median="median")) %>%
  write.csv('SL_Summary_Stats.csv', row.names = F)


Data_summary_Weight <- group_by(Zebra_fish_diet_data_final, Diet) %>%
  summarise_at("Weight", .funs = c(Variance = "var", Mean="mean", SD="sd", median="median")) %>%
  write.csv('Weight_Summary_Stats.csv', row.names = F)


# stats


Weight_summary_table <- read.table("Weight_Summary_Stats.csv", sep = ",", header = TRUE)
SL_summary_table <- read.table("SL_Summary_Stats.csv", sep = ",", header = TRUE)

Summary_table <- rbind(Weight_summary_table, SL_summary_table)

Summary_table$Variable <- c("SL", "SL", "Weight", "Weight")


# flipping columns 
Summary_table <- Summary_table[c("Variable", "Diet", "Variance", "Mean", "SD", "median")]

Summary_table



# part 3 -- test for both continuous variables 

# creating subset DFs for both control and enriched groups
subset_data_control <- subset(Zebra_fish_diet_data_final, Diet=="Control")
subset_data_enriched <- subset(Zebra_fish_diet_data_final, Diet=="Enriched")

#View(subset_data_control)

SL_ttest <- t.test(subset_data_control$SL, subset_data_enriched$SL)
Weight_ttest <- t.test(subset_data_control$Weight, subset_data_enriched$Weight)


print(SL_ttest)
print(Weight_ttest)












# Part 4

# Bootstrap approach addon - function
bootstrap_function <- function(value){
  samples = 1000
  sample_size = 100
  
  # creating a matrix of the 1000 trials among the 20 samples
  boot_sample_matrix <- matrix(sample(value, size = samples * sample_size, replace = TRUE), samples, sample_size)
  
  # calculating the mean for each 1000 samples with 1 row
  boot_stats <- apply(boot_sample_matrix, 1, mean)
  
  # only print out boot_stats values (martix, 1 value, and mean)
  return(boot_stats)
}

# set an open dataframe so I can add the subset columns from initial dataframe for both control and enriched group
bootstrap_DF <- NULL





#View(subset_data_control)

# calling bootstrap function made above for both variables SL & weight to the corresponding subset DFs
bootstrap_DF$bootstrap_control_SL <- bootstrap_function(subset_data_control$SL)
bootstrap_DF$bootstrap_enriched_SL <- bootstrap_function(subset_data_enriched$SL)
bootstrap_DF$bootstrap_control_Weight <- bootstrap_function(subset_data_control$Weight)
bootstrap_DF$bootstrap_enriched_Weight <- bootstrap_function(subset_data_enriched$Weight)



# set as dataframe so I can call it during ggplot
bootstrap_DF <- as.data.frame(bootstrap_DF)


#View(bootstrap_DF)

# histogram of the bootstrap DFs and using column [,1] for aesetics 

Bootstrap_control_weight_Hist <- ggplot(data = bootstrap_DF, aes(bootstrap_DF$bootstrap_control_Weight)) +
  geom_histogram(aes(y = ..density..), alpha = 0.5, colour="black", fill="white", bins = 15) +
  geom_density(alpha=.3, fill = "darkseagreen2") +
  labs(x=expression(bold("mean_dist_weight")), y=expression(bold("Frequency")), title=expression(bold("distribution of weight_control")))

Bootstrap_enriched_weight_Hist <- ggplot(data = bootstrap_DF, aes(bootstrap_DF$bootstrap_enriched_Weight)) +
  geom_histogram(aes(y = ..density..), alpha = 0.5,colour="black", fill="white", bins = 15) +
  geom_density(alpha=.3, fill = "steel blue") +
  labs(x=expression(bold("mean_dist_weight")), y=expression(bold("Frequency")), title=expression(bold("distribution of weight_enriched")))

Bootstrap_control_SL_Hist <- ggplot(data = bootstrap_DF, aes(bootstrap_DF$bootstrap_control_SL)) +
  geom_histogram(aes(y = ..density..), alpha = 0.5,colour="black", fill="white", bins = 30) +
  geom_density(alpha=.3, fill = "steel blue") +
  labs(x=expression(bold("mean_dist_weight")), y=expression(bold("Frequency")), title=expression(bold("distribution of SL_control")))

Bootstrap_enriched_SL_Hist <- ggplot(data = bootstrap_DF, aes(bootstrap_DF$bootstrap_enriched_SL)) +
  geom_histogram(aes(y = ..density..), alpha = 0.5,colour="black", fill="white", bins = 30) +
  geom_density(alpha=.3, fill = "darkseagreen2") +
  labs(x=expression(bold("mean_dist_weight")), y=expression(bold("Frequency")), title=expression(bold("distribution of SL_enriched")))



require(gridExtra)
grid.arrange(Bootstrap_control_weight_Hist, Bootstrap_enriched_weight_Hist, Bootstrap_control_SL_Hist, Bootstrap_enriched_SL_Hist, ncol=2, nrow=2)
  


# basic stats for each column value (4 total) -- these will be used when I plot my resampled values

mean_SL_E <- mean(bootstrap_DF$bootstrap_enriched_SL)
mean_SL_C <- mean(bootstrap_DF$bootstrap_control_SL)


mean_Weight_C <- mean(bootstrap_DF$bootstrap_control_Weight)
mean_Weight_E <- mean(bootstrap_DF$bootstrap_enriched_Weight)


#View(mean_Weight_C)
#View(mean_Weight_E)

#View(mean_SL_C)
#View(mean_SL_E)


SD_SL_E <- sd(bootstrap_DF$bootstrap_enriched_SL)
SD_SL_C <- sd(bootstrap_DF$bootstrap_control_SL)
SD_Weight_C <- sd(bootstrap_DF$bootstrap_control_Weight)
SD_Weight_E <- sd(bootstrap_DF$bootstrap_enriched_Weight)


#View(SD_SL_E)

# confidence interval function
CI_function <- function(value){
  Gene_Standard_err <- sd(value)
  # equation for variable that establishes marginal error (0.25% & 1 - sample size) times value of GSE
  Gene_Margin_err <- qt(.975,999)*Gene_Standard_err
  # Factoring in lower bound of distribution model in this variable (0.25%)
  lower_CI <- mean(value) - Gene_Margin_err
  # factoring in upper bound of distribution model in this variable (0.25%)
  upper_CI <- mean(value) + Gene_Margin_err
  # setting a list so it prints more than one variable when calling my function
  stats_list <- c(lower_CI, upper_CI)
  return(stats_list)
}

# SL CI
CI_SL_E <- CI_function(bootstrap_DF$bootstrap_enriched_SL)
CI_SL_C <- CI_function(bootstrap_DF$bootstrap_control_SL)

# Weight CI
CI_W_E <- CI_function(bootstrap_DF$bootstrap_enriched_Weight)
CI_W_C <- CI_function(bootstrap_DF$bootstrap_control_Weight)

#View(CI_W_E)

# resample data
SL_dataframe <- NULL
Weight_dataframe <- NULL


# Appending the two datasets (SL & Weight) into one column 
SL_dataframe$combined_SL_Data <- c(bootstrap_DF[,1], bootstrap_DF[,2])
Weight_dataframe$combined_weight_Data <- c(bootstrap_DF[,3], bootstrap_DF[,4])  

#View(SL_dataframe)
#View(Weight_dataframe)
  
SL_dataframe <- as.data.frame(SL_dataframe)
Weight_dataframe <- as.data.frame(Weight_dataframe) 

#View(Weight_dataframe)

# part 4 figures with resampled data and stats overlapping hist 

combined_hist_data_SL <- ggplot(data = SL_dataframe, aes(SL_dataframe$combined_SL_Data)) +
  geom_histogram(aes(y = ..density..), alpha = 0.8, colour="black", fill="peachpuff", bins = 30) +
  geom_vline(aes(xintercept= 3.98253, colour = "mean_Control"), linetype="solid", size=1, show.legend=TRUE) +
  geom_vline(aes(xintercept= 4.474164, colour= "mean_Enriched"), linetype="solid", size=1, show.legend = TRUE) +
  geom_vline(aes(xintercept= 4.432707, colour= "C_interval_En_lower"), linetype="dashed", size=0.9, show.legend = TRUE) +
  geom_vline(aes(xintercept= 4.513602, colour= "C_interval_En_upper"), linetype="dashed", size=0.9, show.legend = TRUE) +
  geom_vline(aes(xintercept= 3.939077, colour= "C_interval_Ctrl_lower"), linetype="dashed", size=0.9, show.legend = TRUE) +
  geom_vline(aes(xintercept= 4.025749, colour= "CI_interval_Ctrl_upper"), linetype="dashed", size=0.9, show.legend = TRUE) +
  guides(colour = guide_legend(title=expression(bold("Statistics")))) +
  geom_density(alpha=.3, fill = "darkseagreen2") +
labs(x=expression(bold("mean_distribution_Stand_Length")), y=expression(bold("Frequency")), title=expression(bold("Distribution of SL_combined")))


# ggplot for Weight DF (combined data)

combined_hist_data_Weight <- ggplot(data = Weight_dataframe, aes(Weight_dataframe$combined_weight_Data)) +
  geom_histogram(aes(y = ..density..), alpha = 0.8, colour="black", fill="peachpuff", bins = 30) +
  geom_vline(aes(xintercept= 0.6302481, colour = "mean_Control"), linetype="solid", size=1, show.legend=TRUE) +
  geom_vline(aes(xintercept= 0.6806097, colour= "mean_Enriched"), linetype="solid", size=1, show.legend = TRUE) +
  geom_vline(aes(xintercept= 0.6753320, colour= "C_interval_En_lower"), linetype="dashed", size=0.9, show.legend = TRUE) +
  geom_vline(aes(xintercept= 0.6858875, colour= "C_interval_En_upper"), linetype="dashed", size=0.9, show.legend = TRUE) +
  geom_vline(aes(xintercept= 0.6258413, colour= "C_interval_Ctrl_lower"), linetype="dashed", size=0.9, show.legend = TRUE) +
  geom_vline(aes(xintercept= 0.6346550, colour= "CI_interval_Ctrl_upper"), linetype="dashed", size=0.9, show.legend = TRUE) +
  guides(colour = guide_legend(title=expression(bold("Statistics")))) +
  geom_density(alpha=.3, fill = "seagreen2") +
  labs(x=expression(bold("mean_distribution_Weight")), y=expression(bold("Frequency")), title=expression(bold("Distribution of Weight_combined")))




require(gridExtra)
grid.arrange(combined_hist_data_Weight, combined_hist_data_SL, ncol=2, nrow=1)



# part 5



x <- Zebra_fish_diet_data_final$SL
y <- Zebra_fish_diet_data_final$Weight

#View(x)

# model 1 (four plots supplemental) assumptions we did for model 1
zfish_lm <- lm(y ~ x)

summary(zfish_lm)

plot(y ~ x, col = "blue")
abline(zfish_lm, col = "green")

plot(zfish_lm)




#influence.measures(zfish_lm)

#model 2 

install.packages("lmodel2")
library(lmodel2)

size_lm_modII <- lmodel2(y ~ x)

size_lm_modII


#SMAs is when there is error in y/x with different units in scales but there is a correlation between the two

# y = mx + b (b is equal to a and m is equal to b)

plot(size_lm_modII, xlab = "SL", ylab = "Weight")
abline(a=0.1907835, b=0.1098969, lty=2, col = "blue") # model 2
abline(a=0.216241, b=0.103876, lty=3, lwd = 3, col = "seagreen3") #model 1

#can see that model 1 is more optimal than model 2 due to the line of best fit overlap

# confidence interval for model 1
confint(zfish_lm)




In order to accuretly represent the output, we must identify the level of significance in which we compare our p-value from our test statistic to be at the traditional level of p < 0.05. Furthermore, we will compute a 95% confidence interval for the random variable mean difference μC, to discover approximately how much of a difference there is on average in weight and SL, with or without a controlled diet.


