# Exercise 5: Silohuette coeffiecent
# Javier Galindos
# Ref[1]: https://gist.github.com/AlexandreAbraham/5544803

#rm(list=ls())
library("rstudioapi")  
set.seed(123)

# Set working directory
setwd(dirname(getActiveDocumentContext()$path)) 
getwd()  

# Load distance function
source("JG_ex_1_distance.R")

# loading the data from file
load(file="Data/JGdata.RData")

distFun= "Euclidean"


# The Silhouette Coefficient is calculated using the mean intra-cluster
# distance (a) and the mean nearest-cluster distance (b) for each sample.
# The Silhouette Coefficient for a sample is (b - a) / max(a, b)

intra_cluster_distance_a <- function(dataset, labels,i){
  distance <- array(NA,dim = nrow(dataset))
  n_points <- 0 # Save num of point within the cluster
  # Get the average distance of point xi to points within the cluster it belong to
  for (j in seq(nrow(dataset))){
    if (labels[j] == labels[i]){
      if (i != j){  
        distance[j] <- JG_dist(dataset[i,],dataset[j,],distFun)
        n_points <- n_points + 1
      }
    }
  }
  # Average of distances
  #a <- sum(distance) / n_points
  a <- mean(distance,na.rm=TRUE)
  return (a)
}

nearest_cluster_distance_b <- function(dataset, labels,i){
  n_clusters <- length(unique(labels))
  dist_cluster <- vector()
  clusters <- unique(labels[labels[i] != labels]) # Take out the labels of i
  # Get the average distance to each different cluster 
  for (cluster in clusters){
    distance <- array(NA,dim = nrow(dataset))
    n_points <- 0 # Save num of point within the cluster
    for (j in seq(nrow(dataset))){
        if (labels[j] == cluster){ # If the point belong to the corresponding cluster
            distance[j] <- JG_dist(dataset[i,],dataset[j,],distFun)
            n_points <- n_points + 1
      }
    }
    # Append Average of distances
    dist_cluster <- c(dist_cluster, mean(distance,na.rm=TRUE))
    #dist_cluster[cluster] <- mean(distance,na.rm=TRUE)
  }
  # Return the minimum distance 
  return (min(dist_cluster))
}

JG_silhouette <- function(dataset,labels){
  # Compute the average silhouette coefficient for each cluster
  clusters <- unique(labels)
  # One column for each cluster
  silouette_score <- matrix(,nrow=nrow(dataset), ncol=length(clusters))
  silouette_score_cluster <- vector()
  
  for (cluster in seq_along(clusters)){
    j <- 1
    for (i in seq(nrow(dataset))){
      # If belongs to this cluster
      if (labels[i] == clusters[cluster]){
        # intra-cluster distance (a)
        a <- intra_cluster_distance_a(dataset, labels,i)
        # mean nearest-cluster distance (b)
        b <- nearest_cluster_distance_b(dataset, labels,i)
        silouette_score[j,cluster] <- (b - a) / (max(a,b))
        j <- j + 1 
      }
    }
  }
  silouette_score_cluster <- colMeans(silouette_score,na.rm = TRUE)
  return (silouette_score_cluster)
}

# Test it

#dataset = x
#test <- JG_silhouette(dataset,labels)
