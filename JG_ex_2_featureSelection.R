# Exercise 2: Feature selection
# Javier Galindos
rm(list=ls())
library("rstudioapi")  
set.seed(123)
# Hopkins statistics from scratch

# Set working directory
setwd(dirname(getActiveDocumentContext()$path)) 
getwd()  

# Load distance function
source("JG_ex_1_distance.R")

# loading the data from file
load(file="Data/JGdata.RData")
#load(file="Data/SmallData.RData")

# Paramters for Hopkins Statistics
sample_size <- 150
distFun= "Euclidean"

# Function
JG_hopkinsStats = function(dataset, sample_size){
  # sample of the dataset
  ind <- sample(nrow(dataset), size = sample_size) # nrow() not working for one dimesion
  datasetSampled = dataset[ind,]
  
  distFun= "Euclidean"
  
  #1 Generating synthetic dataset
  uniformGenerator = function(element){
    return( runif(length(element),min(element),max(element)))
  }
  
  syntheticData <- matrix(,nrow=nrow(dataset), ncol=ncol(dataset))
  
  for (i in seq(along=1:ncol(dataset))){
    syntheticData[,i] = uniformGenerator(dataset[,i])
  }
  
  #2 Compute distance to each nearest neighbor
  distNN = function(element, dataset) {
    dist <- matrix(,nrow = nrow(dataset))
    for (i in seq(along=1:nrow(dataset))){
      dist[i] <- JG_dist(element,dataset[i,],distFun)
    }
    dist <- dist[dist>0]
    return(min(dist))
  }
  
  # Sampled dataset
  distSample <- matrix(,nrow = nrow(datasetSampled))
  for (i in seq(along=1:nrow(datasetSampled))){
    distSample[i] <- distNN(datasetSampled[i,],dataset)
  }
  
  # Synthetic dataset
  distSynth <- matrix(,nrow = nrow(datasetSampled))
  for (i in seq(along=1:nrow(datasetSampled))){
    distSynth[i] <- distNN(syntheticData[i,],dataset)
  }
  
  # Compute Hopkins stats
  H = (sum(distSynth))/(sum(distSynth) + sum(distSample))
}

H = JG_hopkinsStats(x,sample_size)
H 

# R function
library(clustertend)
H.R = hopkins(x, sample_size)
H.R

# Entropy from scrach
dataset <- x
JG_entropy <- function(dataset) {
  n_regions <- 10
  
  axis_grid <- matrix(,nrow = ncol(dataset), ncol= (n_regions+1))
  #Create the regions
  for (i in seq(along=1:ncol(dataset))){
    axis_grid[i,] <- seq(from = min(dataset[,i]), to = max(dataset[,i]), by = (abs((max(dataset[,i]) - min(dataset[,i]))/ n_regions)))
  }
  
  grid_ranges <- n_regions**ncol(dataset)
  # Check number of points in each region (2D)
  p <- array(0, dim = grid_ranges)
  horz <- 1
  vert <- 1 
  for (i in seq(along=1:grid_ranges)){
    x_min <- axis_grid[1,horz]
    x_max <- axis_grid[1,horz+1]
    y_min <- axis_grid[2,vert]
    y_max <- axis_grid[2,vert+1]
    # Check if each datapoint belong to this region
    for (j in seq(along=1:nrow(dataset))){
      x <- dataset[j,1]
      y <- dataset[j,2]
      if (x >= x_min && x <= x_max){ # > prev pointx && < next pointx
        if (y >= y_min && y <= y_max){ # >prev point y && < next pointy
          p[i] <- p[i] + 1
        }
      }
    }
    
    # Update x and y
    if (i %in% seq((n_regions),(grid_ranges - n_regions), n_regions)){
      vert <- vert + 1 #Increment y
      horz <- 0
    }
    horz <- horz + 1 #Increment x
  }
  
  p <- p / nrow(dataset)
  p <- p[p>0]
  
  #compute entropy
  E <- -sum(p * log(p) + (1-p) * log(1-p))
  return (E)
}

E = JG_entropy(x)

## Feature selection

# Load silhouette coefficient function
source("JG_ex_5_silhouette.R")

# Dataset
load(file="Data/3gaussiandata.RData")
library(rgl)
plot3d(x[,1], x[,2], x[,3], type="p",col=(x[,4]+1)) 


# Check Hopkins and Entropy in features 1 and 2
E_12= JG_entropy(x[,1:2])
Hopkins_12=JG_hopkinsStats(x[,1:2],sample_size)
Sil_12 = JG_silhouette(x[,1:2],x[,4])

# Check Hopkins and Entropy in features 2 and 3
E_23= JG_entropy(x[,2:3])
Hopkins_23=JG_hopkinsStats(x[,2:3],sample_size)
Sil_23 = JG_silhouette(x[,2:3],x[,4])

# Check Hopkins and Entropy in features 1 and 3
E_13= JG_entropy(x[,c(1,3)])
Hopkins_13=JG_hopkinsStats(x[,c(1,3)],sample_size)
Sil_13 = JG_silhouette(x[,c(1,3)],x[,4])


# 1 & 2 is the best option

# K-means
library(shotGroups)

# Features 1 & 2
train = x[,1:2]
results <- kmeans(train,2)
idx <- results$cluster
plot(x[,1],x[,2], col= (idx+1),main='K-means: Feature 1 and 2')
Sil_12_kmeans = JG_silhouette(train, idx)

# Features 2 & 3
train = x[,2:3]
results <- kmeans(train,2)
idx <- results$cluster
plot(x[,2],x[,3], col= (idx+1),main='K-means: Feature 2 and 3')
Sil_23_kmeans = JG_silhouette(train, idx)

# Features 1 & 3
train = x[,c(1,3)]
results <- kmeans(train,2)
idx <- results$cluster
plot(x[,1],x[,3], col= (idx+1),main='K-means: Feature 1 and 3')
Sil_13_kmeans = JG_silhouette(train, idx)

# All features
train = x[,1:3]
results <- kmeans(train,2)
idx <- results$cluster
Sil_all_kmeans = JG_silhouette(train, idx)

# EM Mixture of gaussians
library(mixtools)

# Features 1 & 2
train = x[,1:2]
gm<-mvnormalmixEM(train,k=2)  #multivariate normal distribution EM
pred<-apply(gm$posterior, 1, function(row) which.max(row))
plot(gm, which=2, main2='EM: Features 1 and 2')
Sil_12_EM = JG_silhouette(train, pred)

# Features 2 & 3
train = x[,2:3]
gm<-mvnormalmixEM(train,k=2)  #multivariate normal distribution EM
pred<-apply(gm$posterior, 1, function(row) which.max(row))
plot(gm, which=2,main2='EM: Features 2 and 3',xlim=c(-4,18),ylim=c(-4,13))
Sil_23_EM = JG_silhouette(train, pred)

# Features 1 & 3
train = x[,c(1,3)]
gm<-mvnormalmixEM(train,k=2)  #multivariate normal distribution EM
pred<-apply(gm$posterior, 1, function(row) which.max(row))
plot(gm, which=2,main2='EM: Features 1 and 3',xlim=c(-6.5,9.5),ylim=c(-2,12))
Sil_13_EM = JG_silhouette(train, pred)

# All features
train = x[,1:3]
gm<-mvnormalmixEM(train,k=2)  #multivariate normal distribution EM
pred<-apply(gm$posterior, 1, function(row) which.max(row))
#plot(gm, which=2)
Sil_all_EM = JG_silhouette(train, pred)





