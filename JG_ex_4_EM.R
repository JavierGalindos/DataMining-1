# Exercise 4: EM 3 dimensional
# Javier Galindos

# Implementation of EM algorithm (3 dimensional)

# clear everything and load required libraries/codes
rm(list=ls())
#library(mvtnorm)
library(mclust) # To compute dmvnorm
library(car)
library(scales)
library("rstudioapi")  
library(rgl)

set.seed(7)

mySum <- function(x) {
  #this function returns sum of the elements of x, whereas infinitly large elements are ignored
  sum(x[is.finite(x)])
}

# Set working directory
 
setwd(dirname(getActiveDocumentContext()$path)) 
getwd()      

# Load test dataset
load("Data/3gaussiandata.RData")
plot3d(x[,1], x[,2], x[,3], type="p") 


# Plotting data
# 2D
plot(x[,1],x[,2], col=scales::alpha(x[,4],0.3), pch=20, xlim=c(-7,10),ylim=c(-10,20)) 
# 3D
plot3d(x[,1], x[,2], x[,3], type="p") 

# data preparation
classes <- x[,4] #labels vector
x <-x[,1:3] #removing labels from data
hist(x, col="blue")
plot(density(x))


#Step 1: initialization
pi1<-0.5  
pi2<-0.5

#initial values for means
mu1<-matrix(c(sample(-10:0, 1),sample(-10:0, 1),sample(1:9,1)), nrow=3) #random uniformly selected starting values
mu2<-matrix(c(sample(0:10, 1),sample(0:10, 1),sample(1:9, 1)),nrow=3) #random uniformly selected starting values

#initial values for covariance matrix for each cluster
#NB!! covariance matrix of sigmas, must be symmetric!!

sym12 = sample(-3:3, 1)
sym13 = sample(-3:3, 1)
sym23 = sample(-3:3, 1)

sigma1<-matrix(c(sample(1:5, 1),sym12,sym13,sym12,sample(1:5, 1),sym23,sym13,sym23,sample(1:5, 1)),nrow=3)

# Sample again to get random numbers
sym12 = sample(-3:3, 1)
sym13 = sample(-3:3, 1)
sym23 = sample(-3:3, 1)
sigma2<-matrix(c(sample(1:5, 1),sym12,sym13,sym12,sample(1:5, 1),sym23,sym13,sym23,sample(1:5, 1)),nrow=3)

#Plotting initizations in the data set
# 2D representation
plot(x[,1],x[,2], col=scales::alpha(classes,0.3), pch=20, xlim=c(-7,10),ylim=c(-10,20)) #plotting data
par(new=TRUE) #to include the previous plot on the previous = combine plots

points(mu1[1], mu1[2], pch=18, cex=1, col="steelblue")
points(mu2[1], mu2[2], pch=18, cex=1, col="yellowgreen")

car::ellipse(center=as.vector(mu1[1:2,]),shape=sigma1[1:2,1:2], radius=3, col="steelblue")
car::ellipse(center=as.vector(mu2[1:2,]),shape=sigma2[1:2,1:2], radius=3, col="yellowgreen")

#initial conditions for stopping the algo
loglik<- rep(NA, 2000) #log likelihoods storage
loglik[1]<-0 #initial log likelihood value
loglik[2]<-mySum(pi1*(log(pi1)+log(matrix(dmvnorm(x,mu1,sigma1),nrow=3,ncol=1000))))+mySum(pi2*(log(pi2)+log(matrix(dmvnorm(x,mu2,sigma2),nrow=3,ncol=1000))))
k<-2

#main loop with step 2, 3 - EM
while(abs(loglik[k]-loglik[k-1]) >= 0.00001) {  #if no significant improvement, finish
  
  # Step 2 -> E-step: Expectation - Calculating the "Soft Labels" of Each Data Point
  tau1<-pi1*matrix(dmvnorm(x,mu1,sigma1))
  tau2<-pi2*matrix(dmvnorm(x,mu2,sigma2))
  
  normalizer<-tau1 + tau2
  
  
  tau1<-tau1/normalizer
  tau2<-tau2/normalizer
  
  # Step 3 -> M step: Maximization - Re-estimate the Component Parameters
  n<-dim(x)[1] #number of datapoints
  
  pi1<-mySum(tau1)/n #recomputing responsabilities
  pi2<-mySum(tau2)/n
  
  # Implementation with for loop
  # Updating mean values
  for (i in seq(along=1:ncol(x))){
    mu1[i]<-(t(tau1)%*%x[,i])/mySum(tau1)
    mu2[i]<-(t(tau2)%*%x[,i])/mySum(tau2)
  }
  # Other implementation
  #mu1[1]<-(t(tau1)%*%x[,1])/mySum(tau1) #recalculating mean values
  #mu1[2]<-(t(tau1)%*%x[,2])/mySum(tau1)  #t(tau) to perform matrix multiplication, row by vector == scalar
  #mu1[3]<-(t(tau1)%*%x[,3])/mySum(tau1)
  
  #mu2[1]<-(t(tau2)%*%x[,1])/mySum(tau2)
  #mu2[2]<-(t(tau2)%*%x[,2])/mySum(tau2)
  #mu2[3]<-(t(tau2)%*%x[,3])/mySum(tau2)
  
  # Updating covariance matrix
  for (i in seq(along=1:ncol(x))){
    for (j in seq(along=1:ncol(x))){
      sigma1[i,j]<-t(tau1)%*%((x[,i]-mu1[i])*(x[,j]-mu1[j]))/(mySum(tau1))
      sigma2[i,j]<-t(tau2)%*%((x[,i]-mu2[i])*(x[,j]-mu2[j]))/(mySum(tau2))
    }
  }
  
  # Other implementation
  # sigma1[1,1]<-t(tau1)%*%((x[,1]-mu1[1])*(x[,1]-mu1[1]))/(mySum(tau1)) #recalculating covariance matrix
  # sigma1[1,2]<-t(tau1)%*%((x[,1]-mu1[1])*(x[,2]-mu1[2]))/(mySum(tau1))
  # sigma1[2,1]<-sigma1[1,2]
  # sigma1[2,2]<-t(tau1)%*%((x[,2]-mu1[2])*(x[,2]-mu1[2]))/(mySum(tau1))
  # 
  # sigma2[1,1]<-t(tau2)%*%((x[,1]-mu2[1])*(x[,1]-mu2[1]))/(mySum(tau2))
  # sigma2[1,2]<-t(tau2)%*%((x[,1]-mu2[1])*(x[,2]-mu2[2]))/(mySum(tau2))
  # sigma2[2,1]<-sigma2[1,2]
  # sigma2[2,2]<-t(tau2)%*%((x[,2]-mu2[2])*(x[,2]-mu2[2]))/(mySum(tau2))
  
  dev.off() #new graph plot
  plot(x[,1],x[,2], col=scales::alpha(classes,0.3), pch=20, xlim=c(-7,10),ylim=c(-10,20)) #plotting data
  points(mu1[1], mu1[2], pch=18, cex=1, col="steelblue")
  points(mu2[1], mu2[2], pch=18, cex=1, col="yellowgreen")
  
  car::ellipse(center=as.vector(mu1[1:2,]),shape=sigma1[1:2,1:2], radius=3, col="steelblue")
  car::ellipse(center=as.vector(mu2[1:2,]),shape=sigma2[1:2,1:2], radius=3, col="yellowgreen")
  
  #new loglik calculation
  loglik[k+1]<-mySum(pi1*(log(pi1)+log(matrix(dmvnorm(x,mu1,sigma1),nrow=3,ncol=1000))))+mySum(pi2*(log(pi2)+log(matrix(dmvnorm(x,mu2,sigma2),nrow=3,ncol=1000))))
  
  k<-k+1
}

# 3D visualization
plot3d(x[,1], x[,2], x[,3], type="p") 
points3d(mu1[1], mu1[2], mu1[3], size=14, col="firebrick1")
points3d(mu2[1], mu2[2], mu2[3], size=14,  col="yellowgreen")
plot3d( ellipse3d(sigma1, centre = mu1), col = "firebrick1", alpha = 0.5, add = TRUE)
plot3d( ellipse3d(sigma2, centre = mu2), col = "yellowgreen", alpha = 0.5, add = TRUE)


# Comparison with library implementation of EM
# library(mixtools)
# gm<-mvnormalmixEM(x,k=2,epsilon=1e-04)  #multivariate normal distribution EM
#                                         #normalmixEM is only for univariate normal
# gm$lambda
# gm$mu
# gm$sigma
# gm$loglik
# plot(gm, which=2)
# 
# (head(gm$posterior))
# pred<-apply(gm$posterior, 1, function(row) which.max(row))
# (confusionMatrix[1,1]+confusionMatrix[2,2])/1000 # accuracy
# (confusionMatrix[1,2]+confusionMatrix[2,1])/1000 # error rate
