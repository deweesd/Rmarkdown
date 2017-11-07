setwd("/Users/Deweesd/Desktop/Advanced_Stats_with_R/PS1")
getwd

# Read in dataset
RNA_data <- read.table("HW_RNAseq.csv", header = TRUE, sep = ",")



# setting categorical (independent) vairables 
Cat_variable <- data.frame(RNA_data[,1:4])


# setting response variables 
response_variable <- data.frame(RNA_data[,5:14])



# creating a new orientation of data set and compiles rows to colums
converted_RNA_seq_data <- reshape2::melt(response_variable)

# plotting each gene column via histagram

# setting value as an representation of gene expression levels
ggplot(data = converted_RNA_seq_data, aes(value)) + 
  # Setting variable to be the header title "gene" and its corresponding number
  facet_wrap(~variable, scales = 'free') +
  # setting width to 12 for each histogram
  geom_histogram(binwidth = 12) +
  # adjusting the theme of each figure by adjusting the text size, colors, etc.
  theme_gray(strip.text.x = element_text(size = 8, angle = 0),
        strip.text.y = element_text(size = 10, face = "bold"),
        strip.background = element_rect(colour = "black", fill = "red")) +

  # adding labels to both X and Y variables 
labs(x=expression(bold("Gene_Exp_Level")), y=expression(bold("Frequency")))


# Transform all original continuous variables to z-scores and plot the distribution for just two genes

response_variable$gene_1_z_score <- (response_variable$Gene01-mean(response_variable$Gene01)) / sd(response_variable$Gene01)
response_variable$gene_2_z_score <- (response_variable$Gene02-mean(response_variable$Gene02)) / sd(response_variable$Gene02)
response_variable$gene_3_z_score <- (response_variable$Gene03-mean(response_variable$Gene03)) / sd(response_variable$Gene03)
response_variable$gene_4_z_score <- (response_variable$Gene04-mean(response_variable$Gene04)) / sd(response_variable$Gene04)
response_variable$gene_5_z_score <- (response_variable$Gene05-mean(response_variable$Gene05)) / sd(response_variable$Gene05)
response_variable$gene_6_z_score <- (response_variable$Gene06-mean(response_variable$Gene06)) / sd(response_variable$Gene06)
response_variable$gene_7_z_score <- (response_variable$Gene07-mean(response_variable$Gene07)) / sd(response_variable$Gene07)
response_variable$gene_8_z_score <- (response_variable$Gene08-mean(response_variable$Gene08)) / sd(response_variable$Gene08)
response_variable$gene_9_z_score <- (response_variable$Gene09-mean(response_variable$Gene09)) / sd(response_variable$Gene09)
response_variable$gene_10_z_score <- (response_variable$Gene10-mean(response_variable$Gene10)) / sd(response_variable$Gene10)

Gene_1_Z_Hist <- ggplot(data = response_variable, aes(gene_1_z_score)) +
  geom_histogram(aes(y=..density..), colour="black", fill="white", bins = 8) +
  geom_density(alpha=.3, fill = "darkseagreen2") +
                 labs(x=expression(bold("Z-Score")), y=expression(bold("Frequency")), title=expression(bold("Gene-01-Z_score")))

Gene_2_Z_Hist <- ggplot(data = response_variable, aes(gene_2_z_score)) +
  geom_histogram(aes(y=..density..), colour="black", fill="white", bins = 8) +
  geom_density(alpha=.3, fill = "darkseagreen2") +
  labs(x=expression(bold("Z-Score")), y=expression(bold("Frequency")), title=expression(bold("Gene-02-Z_score")))

require(gridExtra)
grid.arrange(Gene_1_Z_Hist, Gene_2_Z_Hist, ncol=2)



RNA_data$gene1_zscore <- response_variable$gene_1_z_score

RNA_data$gene2_zscore <- response_variable$gene_2_z_score

# using interaction to categorize and filter through all possible combonation of sex and treatment levels 
RNA_data$sex_treamtment_combination <- interaction(RNA_data$Sex, RNA_data$Treatment)

Gene_1_boxplot_Z <- ggplot(data = RNA_data, aes(y = gene1_zscore, x = sex_treamtment_combination, fill=Population)) + geom_boxplot() + scale_fill_manual(values=c("slateblue1", "olivedrab1")) + theme(axis.text.x = element_text(angle = 40, vjust = 0.3, size = 6)) +
  labs(x=expression(bold("Sex/Treatment")), y=expression(bold("Gene1_Z-Score")), title=expression(bold("Gene1_Z-score_to_Sex/Treatment")))

Gene_2_boxplot_Z <- ggplot(data = RNA_data, aes(y = gene2_zscore, x = sex_treamtment_combination, fill=Population)) + geom_boxplot() + scale_fill_manual(values=c("yellow1", "red1")) + theme(axis.text.x = element_text(angle = 40, vjust = 0.3, size = 6)) +
  labs(x=expression(bold("Sex/Treatment")), y=expression(bold("Gene2_Z-Score")), title=expression(bold("Gene2_Z-score_to_Sex/Treatment")))

Gene_1_boxplot <- ggplot(data = RNA_data, aes(y = Gene01, x = sex_treamtment_combination, fill=Population)) + geom_boxplot() + scale_fill_manual(values=c("slateblue1", "olivedrab1")) + theme(axis.text.x = element_text(angle = 40, vjust = 0.3, size = 6)) +
  labs(x=expression(bold("Sex/Treatment")), y=expression(bold("Gene1_Z-Score")), title=expression(bold("Gene1_Z-score_to_Sex/Treatment")))

Gene_2_boxplot <- ggplot(data = RNA_data, aes(y = Gene02, x = sex_treamtment_combination, fill=Population)) + geom_boxplot() + scale_fill_manual(values=c("yellow1", "red1")) + theme(axis.text.x = element_text(angle = 40, vjust = 0.3, size = 6)) +
  labs(x=expression(bold("Sex/Treatment")), y=expression(bold("Gene2_Z-Score")), title=expression(bold("Gene2_to_Sex/Treatment")))


require(gridExtra)
grid.arrange(Gene_1_boxplot, Gene_1_boxplot_Z, Gene_2_boxplot, Gene_2_boxplot_Z)


# setting table for mean, variance, and SD for the categorical vairables
Gene1_summary <- group_by(RNA_data, Sex, Treatment, Population) %>%
  summarise_at("Gene01", .funs = c(Variance = "var", Mean="mean", SD="sd")) %>%
  write.csv('Gene_1_Summary_Stats.csv', row.names = F)

# setting table for mean, variance, and SD for the categorical vairables
Gene2_summary <- group_by(RNA_data, Sex, Treatment, Population) %>%
  summarise_at("Gene02", .funs = c(Variance = "var", Mean="mean", SD="sd")) %>%
  write.csv('Gene_2_Summary_Stats.csv', row.names = F)

# setting table for mean, variance, and SD for the categorical vairables
Gene3_summary <- group_by(RNA_data, Sex, Treatment, Population) %>%
  summarise_at("Gene03", .funs = c(Variance = "var", Mean="mean", SD="sd")) %>%
  write.csv('Gene_3_Summary_Stats.csv', row.names = F)

# setting table for mean, variance, and SD for the categorical vairables
Gene4_summary <- group_by(RNA_data, Sex, Treatment, Population) %>%
  summarise_at("Gene04", .funs = c(Variance = "var", Mean="mean", SD="sd")) %>%
  write.csv('Gene_4_Summary_Stats.csv', row.names = F)

# setting table for mean, variance, and SD for the categorical vairables
Gene5_summary <- group_by(RNA_data, Sex, Treatment, Population) %>%
  summarise_at("Gene05", .funs = c(Variance = "var", Mean="mean", SD="sd")) %>%
  write.csv('Gene_5_Summary_Stats.csv', row.names = F)

# setting table for mean, variance, and SD for the categorical vairables
Gene6_summary <- group_by(RNA_data, Sex, Treatment, Population) %>%
  summarise_at("Gene06", .funs = c(Variance = "var", Mean="mean", SD="sd")) %>%
  write.csv('Gene_6_Summary_Stats.csv', row.names = F)

# setting table for mean, variance, and SD for the categorical vairables
Gene7_summary <- group_by(RNA_data, Sex, Treatment, Population) %>%
  summarise_at("Gene07", .funs = c(Variance = "var", Mean="mean", SD="sd")) %>%
  write.csv('Gene_7_Summary_Stats.csv', row.names = F)

# setting table for mean, variance, and SD for the categorical vairables
Gene8_summary <- group_by(RNA_data, Sex, Treatment, Population) %>%
  summarise_at("Gene08", .funs = c(Variance = "var", Mean="mean", SD="sd")) %>%
  write.csv('Gene_8_Summary_Stats.csv', row.names = F)

# setting table for mean, variance, and SD for the categorical vairables
Gene9_summary <- group_by(RNA_data, Sex, Treatment, Population) %>%
  summarise_at("Gene09", .funs = c(Variance = "var", Mean="mean", SD="sd")) %>%
  write.csv('Gene_9_Summary_Stats.csv', row.names = F)

# setting table for mean, variance, and SD for the categorical vairables
Gene10_summary <- group_by(RNA_data, Sex, Treatment, Population) %>%
  summarise_at("Gene10", .funs = c(Variance = "var", Mean="mean", SD="sd")) %>%
  write.csv('Gene_10_Summary_Stats.csv', row.names = F)

# setting parametic function for SE and 95% conf. interval (assuming distribution is normal = PA)
parametric_approach_func <- function(value){
  # equation for variable that establistes standard error (value represents gene integers in table)
  Gene_Standard_err <- sd(value)/sqrt(20)
  # equation for variable that establishes marginal error (0.25% & 1 - sample size) times value of GSE
  Gene_Margin_err <- qt(.025,19)*Gene_Standard_err
  # Factoring in lower bound of distribution model in this variable (0.25%)
  lower_bound <- mean(value) - Gene_Margin_err
  # factoring in upper bound of distribution model in this variable (0.25%)
  upper_bound <- mean(value) + Gene_Margin_err
  # setting a list so it prints more than one variable when calling my function
  stats_list <- c(Gene_Standard_err, lower_bound, upper_bound)
  return(stats_list)
}

# This is where I calculate the parametic value function onto the 10 established gene columns set in RNA datatable

PA_Gene1_mean <- parametric_approach_func(RNA_data$Gene01)
PA_Gene2_mean <- parametric_approach_func(RNA_data$Gene02)
PA_Gene3_mean <- parametric_approach_func(RNA_data$Gene03)
PA_Gene4_mean <- parametric_approach_func(RNA_data$Gene04)
PA_Gene5_mean <- parametric_approach_func(RNA_data$Gene05)
PA_Gene6_mean <- parametric_approach_func(RNA_data$Gene06)
PA_Gene7_mean <- parametric_approach_func(RNA_data$Gene07)
PA_Gene8_mean <- parametric_approach_func(RNA_data$Gene08)
PA_Gene9_mean <- parametric_approach_func(RNA_data$Gene09)
PA_Gene10_mean <- parametric_approach_func(RNA_data$Gene10)


parametric_approach_table <- rbind(PA_Gene1_mean, PA_Gene2_mean, PA_Gene3_mean, PA_Gene4_mean, PA_Gene5_mean, PA_Gene6_mean, PA_Gene7_mean, PA_Gene8_mean, PA_Gene9_mean, PA_Gene10_mean)

#View(parametric_approach_table)



# Bootstrap approach addon

bootstrap_function <- function(value){
  samples = 1000
  sample_size = 20
  
  # creating a matrix of the 1000 trials among the 20 samples
  boot_sample_matrix <- matrix(sample(value, size = samples * sample_size, replace = TRUE),samples, sample_size)
  
  # calculating the mean for each 1000 samples with 1 row
  boot_stats <- apply(boot_sample_matrix, 1, mean)
  
  Gene_Standard_err <- sd(boot_stats)
  # equation for variable that establishes marginal error (0.25% & 1 - sample size) times value of GSE
  Gene_Margin_err <- qt(.025,19)*Gene_Standard_err
  # Factoring in lower bound of distribution model in this variable (0.25%)
  lower_bound <- mean(value) - Gene_Margin_err
  # factoring in upper bound of distribution model in this variable (0.25%)
  upper_bound <- mean(value) + Gene_Margin_err
  # setting a list so it prints more than one variable when calling my function
  stats_list <- c(Gene_Standard_err, lower_bound, upper_bound)
  return(stats_list)
}


Boot_Gene1_mean <- bootstrap_function(RNA_data$Gene01)
Boot_Gene2_mean <- bootstrap_function(RNA_data$Gene02)
Boot_Gene3_mean <- bootstrap_function(RNA_data$Gene03)
Boot_Gene4_mean <- bootstrap_function(RNA_data$Gene04)
Boot_Gene5_mean <- bootstrap_function(RNA_data$Gene05)
Boot_Gene6_mean <- bootstrap_function(RNA_data$Gene06)
Boot_Gene7_mean <- bootstrap_function(RNA_data$Gene07)
Boot_Gene8_mean <- bootstrap_function(RNA_data$Gene08)
Boot_Gene9_mean <- bootstrap_function(RNA_data$Gene09)
Boot_Gene10_mean <- bootstrap_function(RNA_data$Gene10)

Boot_strap_table <- rbind(parametric_approach_table, Boot_Gene1_mean, Boot_Gene2_mean, Boot_Gene3_mean, Boot_Gene4_mean, Boot_Gene5_mean, Boot_Gene6_mean, Boot_Gene7_mean, Boot_Gene8_mean, Boot_Gene9_mean, Boot_Gene10_mean)


bootstrap_var_function <- function(value){
  samples = 1000
  sample_size = 20
  
  # creating a matrix of the 1000 trials among the 20 samples
  boot_sample_matrix <- matrix(sample(value, size = samples * sample_size, replace = TRUE),samples, sample_size)
  
  # calculating the mean for each 1000 samples with 1 row
  boot_stats <- apply(boot_sample_matrix, 1, var)
  
  Gene_Standard_err <- sd(boot_stats)
  # equation for variable that establishes marginal error (0.25% & 1 - sample size) times value of GSE
  Gene_Margin_err <- qt(.025,19)*Gene_Standard_err
  # Factoring in lower bound of distribution model in this variable (0.25%)
  lower_bound <- mean(value) - Gene_Margin_err
  # factoring in upper bound of distribution model in this variable (0.25%)
  upper_bound <- mean(value) + Gene_Margin_err
  # setting a list so it prints more than one variable when calling my function
  stats_list <- c(Gene_Standard_err, lower_bound, upper_bound)
  return(stats_list)
}

Boot_Gene1_var <- bootstrap_var_function(RNA_data$Gene01)
Boot_Gene2_var <- bootstrap_var_function(RNA_data$Gene02)
Boot_Gene3_var <- bootstrap_var_function(RNA_data$Gene03)
Boot_Gene4_var <- bootstrap_var_function(RNA_data$Gene04)
Boot_Gene5_var <- bootstrap_var_function(RNA_data$Gene05)
Boot_Gene6_var <- bootstrap_var_function(RNA_data$Gene06)
Boot_Gene7_var <- bootstrap_var_function(RNA_data$Gene07)
Boot_Gene8_var <- bootstrap_var_function(RNA_data$Gene08)
Boot_Gene9_var <- bootstrap_var_function(RNA_data$Gene09)
Boot_Gene10_var <- bootstrap_var_function(RNA_data$Gene10)

boot_var_table <- rbind(Boot_strap_table, Boot_Gene1_var, Boot_Gene2_var, Boot_Gene3_var, Boot_Gene4_var, Boot_Gene5_var, Boot_Gene6_var, Boot_Gene7_var, Boot_Gene8_var, Boot_Gene9_var, Boot_Gene10_var)

colnames(boot_var_table) <- c("Standard Error", "Mean", "Variance")

pander(boot_var_table)


#Part 3

Num_of_trees_obs_in_plot <- c(0,1,2,3,4,5,6,7,8,9,10,11)
Num_of_plots_with_this_many_trees <- c(74,149,228,181,169,84,49,24,19,12,9,4)

#making data frame of said vectors generated above 
tree_plot_dataframe <- data.frame(Num_of_trees_obs_in_plot, Num_of_plots_with_this_many_trees)

View(tree_plot_dataframe)

tree_plot_dataframe$trees <- tree_plot_dataframe$Num_of_trees_obs_in_plot*tree_plot_dataframe$Num_of_plots_with_this_many_trees


tree_plot_dataframe$sum <- sum(tree_plot_dataframe$trees)/1000

View(tree_plot_dataframe)

mean(tree_plot_dataframe$sum)
var(tree_plot_dataframe$Num_of_trees_obs_in_plot)

ggplot(data = tree_plot_dataframe) + 
  geom_col(aes(x = tree_plot_dataframe[,1], y = tree_plot_dataframe[,2])) +
  geom_area(aes(x = tree_plot_dataframe[,1], y = tree_plot_dataframe[,2]), alpha = 0.3, fill = "yellow") +
  geom_smooth(aes(x = tree_plot_dataframe[,1], y = tree_plot_dataframe[,2]), color = "black", fill = "steel blue") +
labs(x=expression(bold("Tree_obs_in_plot")), y=expression(bold("Num_of_plot/tree_freq")), title=expression(bold("Plot in Tree Distribution")))




mle_list <- rep.int(tree_plot_dataframe$Num_of_trees_obs_in_plot, tree_plot_dataframe$Num_of_plots_with_this_many_trees)

mean_of_trees <- mean(mle_list)
var_of_trees <- var(mle_list)

# set of parameters
lambda <- seq(1,10,by=0.001)

ln_poisson_probs <- log((exp(1)^-lambda)*(lambda^mean_of_trees)/factorial(mean_of_trees))


print(ln_poisson_probs)


plot(lambda, ln_poisson_probs, type = "l")
LNL_function_dataframe <- data.frame(lambda, ln_poisson_probs)

Confidence_interval <- subset(_func, ln_poisson_probs>=max(ln_poisson_probs)-0.00192)

print(ci_interval)
