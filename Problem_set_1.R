#RNAseq_Data <- read.table('GacuRNAseq_subset.csv', header=T,sep=',')
#print(RNAseq_Data)
#head(RNAseq_Data)
#dim(RNAseq_Data) #20 186

#RNAseq_Data_txt <- read.table('GacuRNAseq_Subset.txt', header=T,sep=',')
#print(RNAseq_Data_txt)
#head(RNAseq_Data_txt)
#dim(RNAseq_Data_txt) #20 1

setwd ("/Users/Deweesd/Desktop/Bi623_assignments")
getwd()

data_csv <- read.table('GacuRNAseq_subset.csv', header=T, row.names=1, sep=',')
data_txt <- read.table('GacuRNAseq_Subset.txt', header=T, row.names=1, sep='\t')

#head(data_csv)
dim(data_csv) #20 185

#head(data_txt)
dim(data_txt) #20 185

class(data_txt) #data.frame

#explanatory variables: population, treatment, and sex (dependent variables are genes)

unique(data_txt[,1])
#Levels Boot RabbitSlough

unique(data_txt[,2])
#conventional monoAssoc

unique(data_txt[,3])
#female male 

class(data_txt[,1])
#factor class type 

class(data_txt[,2])
#factor class type

class(data_txt[,3])
#factor class type

#182 response variables (genes) 185 total variables (-3 to exclude explanatory variables)

#use class() in a for loop to iterate over all response variable columns
for (i in data_txt[,4:185]){
  print(class(i))
}
#all columns are numberic (182)

mean_col_1 <- data_txt[,4]
mean(mean_col_1)
#mean is 20.0366

mean_col_2 <- data_txt[,5]
mean(mean_col_2)
#mean is 5.014469

#help(lapply)
lapply(data_txt[,4:5], mean)

#$ENSGACG00000000003
#[1] 20.0366

#$ENSGACG00000000004
#[1] 5.014469
expr_means_100 <- lapply(data_txt[,4:103], mean)

mean(expr_means_100)
class(expr_means_100)
#argument is not numeric or logical --> need to turn lists into numerical values

#6
expr_means_100_vector <- as.numeric(expr_means_100) #convert list to numerical vector
class(expr_means_100_vector)

mean(expr_means_100_vector)
#mean is 314.8202
class(expr_means_100_vector)


hist(expr_means_100_vector, main="mean expression level", xlab="mean expression level")

pdf("hist_mean_exp.pdf")

#7
expr_means_100_vector_log_10 <- log10(expr_means_100_vector)

hist(expr_means_100_vector_log_10, main="log 10 mean expression level for 100 genes", xlab="log 10 mean expression level")

pdf("hist_logmean_exp.pdf")

#8
sub_less_than_500 <- subset(expr_means_100_vector, expr_means_100_vector<500) #take values (means) less than 500
hist(sub_less_than_500, main="mean<500", xlab="mean expression level", col = "steel blue")

pdf("hist_subset_exp.pdf")

#9

?tapply
#do quantitative list first and explanatory variable
tapply(data_txt[,5], data_txt[,1], mean)#apply function to each cell of a ragged array

# Boot RabbitSlough 
#9.8480121    0.1809256

tapply(data_txt[,5], data_txt[,1], sd)
#SD Boot RabbitSlough 
#5.6041829    0.3815023

mean(data_txt[,5])
#5.014469

sd(data_txt[,5])
#6.287977 total 

boot <- rnorm(1000, 9.848012, 5.6041829) #creates normal distribution
RabbitSlough <- rnorm(1000, 0.1809256, 0.3815023)
Rab_Boots <- rnorm(1000, 5.014469, 6.287977)

par(mfrow = c(3,1)) #able to put three plots in one variable 

hist(boot, main="Boot Lake Distribution", xlab="boot expression level", xlim = c(-15,30), ylim = c(0,300), col = "steel blue")
hist(RabbitSlough, main="Rabbit Lake Distribution", xlab="Rabbit expression level", xlim = c(-15,30), ylim = c(0,300),  col = "red")
hist(Rab_Boots, main="Total Lake Distribution", xlab="Total expression level", xlim = c(-15,30), ylim = c(0,300), col = "green")

pdf("Gene2_Pop_NormSamp.pdf")

#10
three_genes <- data.frame(data_txt[,1:6])

write.table(three_genes, "three_genes.csv", quote=F, row.names=T, sep=",")

#file checked out when ran on terminal

#11
coefvar <- function(v){
  vec <- c() #open set for vec
  for(x in v) {
    m <- mean(x)
    std <- sd(x)
    covar <- (std/m)
    vec <- c(vec, covar) #concatonate the values from vec and covar
  }
  vec 
}


co_variance <- coefvar(data_txt[,4:185])

par(mfrow = c(1,1)) #set just to plot 1 fifure
hist(co_variance, main="Covariance Distribution", xlab="Covariance", col = "steel blue")

boxplot(co_variance, main="Covariance Distribution", xlab="Covariance", col = "red")
#dev.off()

#12

subset_Boot_pop <- subset(data_txt, Population == "Boot") #generating table of just Boot and gene values
subset_Rabbit_pop <- subset(data_txt, Population == "RabbitSlough")

cov_values_Boot <- coefvar(subset_Boot_pop[,4:185]) #setting new vector that applies my function "coefvar" to the subset vector I made earlier with all gene values
cov_values_Rabbit <- coefvar(subset_Rabbit_pop[,4:185])

#histo and boxplots for both cov_values of rabbitslough and boot

par(mfrow = c(2,2))
hist(cov_values_Boot, main="Histogram CV Boot Distribution", xlab="Covariance", col = "blue")
hist(cov_values_Rabbit, main="Histogram CV RabbitSlough Distribution", xlab="Covariance", col = "red" )
boxplot(cov_values_Boot, main="Covariance Distribution of Boot pop", xlab="Covariance", col = "steel blue")
boxplot(cov_values_Rabbit, main="Covariance Distribution of RabbitSlough Pop", xlab="Covariance", col = "steel blue")



#13

#histogram for rabbitslough has a wider spread of cov compared to boot as well as a larger frequency at the 0.5 position of the x-axis




#14

#CV is a better representaion over std because it has a larger spread over the mean values compared to std. You can compare coefficient 
#varient to do different means. It normalizes the standard deviations





#cv_vector_Boot = c() #open vector to append the CV values at Boot
#cv_vector_Rabbit = c()
#response_variable <- colnames(data_txt[,4:185]) #variable of all the gene values

#for (gene_header in response_variable)
  #if (is.numeric(subset_Boot_pop[[gene_header]])) {
   #

#}

#View(subset_Boot_pop)


#15
#create 4 sample test and use rnorm to show sample size(n), mean, and std
sample_1 <- rnorm(10, 20, 0.3)
sample_2 <- rnorm(10, 21, 0.3)

sample_3 <- rnorm(1000, 20, 0.5)
sample_4 <- rnorm(1000, 21, 0.5)

#1 row, 2 columns
par(mfrow = c(1,2))

#print out boxplot figures and use cbind to allow two plots on one plane.
boxplot(cbind(sample_1, sample_2), ylab ="Mass (grams)", xlab = "N = 10       N = 10", ylim = c(18,24))
boxplot(cbind(sample_3, sample_4), ylab ="Mass (grams)", xlab = "N = 10       N = 10", ylim = c(18,24))


#generate ttest value (p<) for each sample
t.test(sample_1)
t.test(sample_2)
t.test(sample_3)
t.test(sample_4)

dev.off()
