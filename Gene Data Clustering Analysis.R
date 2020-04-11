#Applied Machine Learning for Health Data 
#Gene Data Clustering Analysis
#Name: Ivy Fong 
#Date: December 20, 2018

#load packages to be used
library(cluster)

#set working directory 
setwd("C:/Users/ivyfo/Dropbox/Master of Public Health/Master of Public Health - Courses/Fall 2018 - Courses/CHL7001 - Machine Learning/CHL7001 - Assignments/CHL7001 A3")

#set seed
set.seed(123)

#load data and create dataset d without missing values
d <- read.csv("Data_Cortex_Nuclear.csv", header=T, na.strings="?") #read csv data into R, specify variable names in header, return ? for missing values
d <- na.omit(d) #only keep observations with complete information
summary(d) #print summary of d dataset

#subset data
d.data <- d[,2:78] #create subset d.data gene expression variables
d.labels <- d[,82] #create subset d.labels with class variable
summary(d.data) #print summary of d.data dataset
summary(d.labels) #print summary of d.labels dataset

#calculate mean and variance of the gene expression variables for each observation
apply(d.data, 2, mean) #print mean of gene variables
apply(d.data, 2, var) #print variance of gene variables
#scale of features looks similar - will not scale data


#PAM
#run PAM algorithm for different number of clusters and compare the associated silhouette widths
sil_width <- c() #create sil_width vector

for(i in 1:9){
  pam_fit <- pam(d.data, k=i+1) #run PAM algorithm for 2-10 clusters
  sil_width[i] <- pam_fit$silinfo$avg.width #extract average silhouette width from silinfo$avg.width to build vector
}

plot(2:10, sil_width, #plot sihouette width (higher is better)
     xlab = "Number of clusters",
     ylab = "Silhouette Width")
lines(2:10, sil_width) #highest silhouette width with 3 clusters 

#apply PAM to partition gene data into 3 clusters around medoids - a more robust version of K-means
pam <- pam(d.data, 3) #apply PAM to gene data with k = 3 clusters
pam$clusinfo #print PAM cluster information including cluster size = 169 220 163
pam$silinfo$clus.avg.widths #print cluster average silhouette widths = 0.1851879 0.3212409 0.2780063
pam$silinfo$avg.width #print average silhouette width = 0.2668203
  
#create tables to present the distribution of outcomes across clusters
table(pam$cluster,d.labels) #create table of clusters by class
table(pam$cluster,d$Genotype) #create table of clusters by genotype
table(pam$cluster,d$Treatment) #create table of clusters by treatment
table(pam$cluster,d$Behavior) #create table of clusters by behaviour

#create silhouette plot
si <- silhouette(pam) #use silhouette function to compute silhouette information for PAM clustering
plot(si) #generate silhouette plot


#K-means
#run K-means algorithm for different number of clusters and compare the associated silhouette widths
km2 <- kmeans(d.data, centers=2, nstart=20)
si.km2 <- silhouette(km2$cluster, dist(d.data))
summary(si.km2) #k = 2, silhouette width = 0.2752

km3 <- kmeans(d.data, centers=3, nstart=20)
si.km3 <- silhouette(km3$cluster, dist(d.data))
summary(si.km3) #k = 3, silhouette width = 0.27490

km4 <- kmeans(d.data, centers=4, nstart=20)
si.km4 <- silhouette(km4$cluster, dist(d.data))
summary(si.km4) #k = 4, silhouette width = 0.21944

km5 <- kmeans(d.data, centers=5, nstart=20)
si.km5 <- silhouette(km5$cluster, dist(d.data))
summary(si.km5) #k = 5, silhouette width = 0.20543

km6 <- kmeans(d.data, centers=6, nstart=20)
si.km6 <- silhouette(km6$cluster, dist(d.data))
summary(si.km6) #k = 6, silhouette width = 0.21206

km7 <- kmeans(d.data, centers=7, nstart=20)
si.km7 <- silhouette(km7$cluster, dist(d.data))
summary(si.km7) #k = 7, silhouette width = 0.22344

km8 <- kmeans(d.data, centers=8, nstart=20)
si.km8 <- silhouette(km8$cluster, dist(d.data))
summary(si.km8) #k = 8, silhouette width = 0.2123

km9 <- kmeans(d.data, centers=9, nstart=20)
si.km9 <- silhouette(km9$cluster, dist(d.data))
summary(si.km9) #k = 9, silhouette width = 0.21825

km10 <- kmeans(d.data, centers=10, nstart=20)
si.km10 <- silhouette(km10$cluster, dist(d.data))
summary(si.km10) #k = 10, silhouette width = 0.2085

#create sil_width2 vector with silhouette widths for K-means clustering with 2 to 10 clusters
sil_width2 <- c(0.2752, 0.27490, 0.21944, 0.20543, 0.21206, 0.22344, 0.2123, 0.21825, 0.2085)

#plot sihouette width (higher is better)
plot(2:10, sil_width2, 
     xlab = "Number of clusters",
     ylab = "Silhouette Width")
lines(2:10, sil_width2)  #highest silhouette width with 2 clusters

#apply K-means to partition gene data into 2 clusters around centroids
km <- kmeans(d.data, centers=2, nstart=20) #perform K-means clustering 20 times, create 2 clusters, pick the best solution
km$iter #print number of iterations needed to find solution = 1

#create tables to present the distribution of outcomes across clusters
table(km$cluster,d.labels) #create table of clusters by class
table(km$cluster,d$Genotype) #create table of clusters by genotype
table(km$cluster,d$Treatment) #create table of clusters by treatment
table(km$cluster,d$Behavior) #create table of clusters by behaviour
#no distinct clusters by genotype, behaviour, treatment, or class

#create silhouette plot
si2 <- silhouette(km$cluster, dist(d.data)) #compute silhouette information, dist matrix is needed - dist outputs the distance between 1st and 2nd point and so on
summary(si2) #print silhouette information summary
#cluster size = 265 287 
#cluster average silhouette widths = 0.2094640 0.3359831 
#average silhouette width = 0.2752
plot(si2) #generate silhouette plot


#PCA
#perform PCA
pr.out <- prcomp(d.data, scale=F) #perform PCA on gene data, don't scale the data
pr.out$rotation[,1] #print pc1 loadings - except NR2A_N, ERK_N, pCAMKII_N, Bcatenin_N, loadings > 0.2, all other pc1 loadings are close to 0
pr.out$rotation[,2] #print pc2 loadings - except NR2A_N, ERK_N, pCAMKII_N loadings > 0.2, all other pc2 loadings are close to 0

#create biplot
biplot(pr.out, scale=0) #generate biplot, plotting together the points and the features based on the first 2 pc's

#calculate pve = proportion of variance explained by each component
pr.var <- pr.out$sdev^2
pve <- pr.var/sum(pr.var) 
pve #print proportion of variance explained by each component

#plot proportion of variance explained
plot(pve, xlab="Principal Component", ylab="Proportion of
     Variance Explained ", ylim=c(0,1) ,type="b")

#plot cumulative proportion of variance explained
plot(cumsum(pve), xlab="Principal Component", ylab="
     Cumulative Proportion of Variance Explained ", ylim=c(0,1) ,
     type="b")

#cumsum function gives cumulative sums of previous units - proportion of variance explained
cumsum(pve) #print cumulative proportion of variance explained - first 6 pc's explain about 92% of the data variation

#colour-code PAM clusters on plot of pc1 vs. pc2
plot(pr.out$x[,1:2], col=4-as.numeric(pam$cluster)) #generate plot of pc2 vs. pc2 with colour-coded PAM clusters
legend("bottomright", legend=levels(as.factor(pam$cluster)), text.col=4-(1:3), y.intersp=0.8) #add legend with cluster ID
