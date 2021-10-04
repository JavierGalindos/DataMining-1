# Exercise 3: K-means
# Javier Galindos

# K-means from scrach with stopping criteria
# this file is to implement my own clustering 
rm(list=ls())
graphics.off()
library("rstudioapi")  
set.seed(123)

# Set working directory
setwd(dirname(getActiveDocumentContext()$path)) 
getwd()  

# loading the data from file
load(file="Data/JGdata.RData")
#load(file="Data/2gaussian.RData")

# Load distance function
source("JG_ex_1_distance.R")

# function to find the label
get_label <- function(dataset, centroids){
  for (i in seq(along=1:nrow(dataset))){
    temp_dist <- matrix(,nrow=nrow(centroids))
    for (j in seq(along=1:nrow(centroids))){
      temp_dist[j,] = JG_dist(centroids[j,], dataset[i,1:2], "Euclidean")
    }
    dataset[i, 3] = which.min(temp_dist)
  }
  return(dataset)
}

# update centroids on each iteration
update_centroids <- function(dataset){
  last_column = ncol(dataset)
  K = max(dataset[,last_column])
  centroids <- matrix(, nrow=K, ncol=last_column-1)
  for (k in seq(along=1:K)){
    cluster <- dataset[dataset[, last_column] == k, 1:2]
    
    centroids[k,] = colMeans(cluster)
    
  }
  cat("___")
  cat(centroids)
  cat("___ \n")
  return(centroids)
}


# do we need the following step?
#split the data in proportion 70/30 for training and validation purposes.
sample_size <- floor(0.7 * nrow(x))
data_dimensionality = ncol(x)
#set.seed(123) # seed is necessary to make it reproducible
train_ind <- sample(seq_len(nrow(x)), size = sample_size)

train_set <- x[train_ind, ]
test <- x[-train_ind, ]
train <- train_set[,1:2] # just to assure that we have two columns

# initialize
# Assume that K = 3
K <- 3 
centroids <- matrix(, nrow=K,ncol=data_dimensionality)


for (k in seq(along=1:K)){
  centroids[k,] = runif(2)
}
cat("Initial centroids are: ", centroids)

# this example does not contain any convergence criteria
# students are asked to design and add it.

# add one column (to be used for labels)
aux_column <- matrix(0, nrow(train_set))
train_set = cbind(train_set,aux_column)

# initialize 
#centroids <- matrix(, nrow=K,ncol=data_dimensionality)
step = 1

# Add stopping criteria: 1. No variation in centroids or 2. More than 100 iterations
while (step < 100){
  cat("Step: ", step)
  train_set <- get_label(train_set, centroids)
  centroids_prev <- centroids
  centroids <- update_centroids(train_set)
  cat(" \n Centroids have been updated: ", centroids)
  #plot(centroids, col="red")
  #par(new=TRUE)
  if (identical(centroids_prev, centroids)) {
    break
  }
  step = step + 1
}

plot(x[,1],x[,2], main = 'JG Dataset', pch= 20)
plot(train_set[, 1], train_set[, 2], type='p',col=(train_set[,3]+1),main='K-means output')
par(new=TRUE)
points(centroids[,1], centroids[,2], col="black", type='p',pch=20)


