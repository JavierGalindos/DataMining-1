# Exercise 6: Outlier analysis
# Javier Galindos

# Local outlier factor (LOF)

rm(list=ls())

# Libraries
library(DDoutlier) # LOF
# Set working directory
setwd(dirname(getActiveDocumentContext()$path)) 
getwd()  

# loading the data from file
load(file="Data/JGdata.RData")

# Load distance function
source("JG_ex_1_distance.R")

# Distance metric
distFun= "Euclidean"

# Reachability distance
reachability_distance <- function (dataset, X, Y, v_k){
  s_x_y <- JG_dist(dataset[X,], dataset[Y,], distFun)
  return (max(s_x_y,v_k))
}

# This function return the LOF score for each data point
JG_LOF <- function(dataset, k){
  # Get the distance to k nearest neighbor of every point (knn_distance) 
  # and index of k-neighbors (neighbors)
  knn_distance <- matrix(,nrow = nrow(dataset)) # Create the array for the distance
  neighbors <- matrix(,nrow = nrow(dataset), ncol = k) # Array for index of neighbors
  for (i in seq(along=1:nrow(dataset))){
    temp_dist <- matrix(,nrow = nrow(dataset))
    # Compute distance to each point
    for (j in seq(along=1:nrow(dataset))){
      temp_dist[j] = JG_dist(dataset[i,], dataset[j,], distFun)
    }
    # Take the k lowest distance
    temp_dist <- sort(temp_dist,index.return=TRUE)
    knn_distance[i] <- temp_dist$x[k+1] # First element is the same point (k+1) to skip the same point
    # Save the index of the neighbors of distance less than knn_distance
    # For better understanding of the code, only k point are taken.
    neighbors[i,] <- temp_dist$ix[2:(k+1)]
  }
  
  # Average reachability distance
  AR_x <- matrix(0,nrow = nrow(dataset))
  for (i in seq(along=1:nrow(dataset))){
    # average of its reachability distances to all objects in its neighborhood.
    for (j in seq(along=1:k)){
      AR_x[i] <- AR_x[i] + reachability_distance(dataset, i, neighbors[i,j], knn_distance[neighbors[i,j]])
    }
    AR_x[i] <- AR_x[i] / ncol(neighbors)
  }
  
  # LOF
  LOF <- matrix(0,nrow = nrow(dataset))
  for (i in seq(along=1:nrow(dataset))){
    for (j in seq(along=1:k)){
      LOF[i] <- LOF[i] + AR_x[i] / AR_x[neighbors[i,j]]
    }
    LOF[i] <- LOF[i] / ncol(neighbors)
  }
  
  return(LOF)
}

LOF <- JG_LOF(x,3)
maxColorValue <- 50
palette <- colorRampPalette(c("white","black"))(maxColorValue)
par(bg = 'white')
plot(x[,1],x[,2],col=palette[cut(LOF, maxColorValue)],pch=20,main = 'LOF Output')


LOF=sort(LOF,decreasing = TRUE)

# R implementation
lof <- LOF(x, k = 3)

lof=sort(lof,decreasing = TRUE)
