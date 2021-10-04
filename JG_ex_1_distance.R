# Exercise 1: Disrance funtions
# Javier Galindos

# Function to compute distances
# this function returns the distace between the element1 and element2
# according to the metricf
# Euclidean, Manhattan, Minkowski, Mahalanobis, Cosine and Canberra distances
# are implemented 

JG_dist<- function(element1,element2,metricf,p,dataset){
  dimensions=length(element1)
  sqd<-matrix(,dimensions,1)
  
  # Euclidean
  if (metricf=="Euclidean"){
    for(i in seq(along=element1)){
      sqd[i]<-(element1[i]-element2[i])^2
    }
    dist<-sqrt(colSums(sqd))
  }
  
  # Manhattan
  if (metricf=="Manhattan"){
    for(i in seq(along=element1)){
      sqd[i]<- abs(element1[i]-element2[i])
    }
    dist<-colSums(sqd)
  }
  
  # Minkowski for any p value
  if (metricf=="Minkowski"){
    for(i in seq(along=element1)){
      sqd[i]<- abs(element1[i]-element2[i])^p
    }
    dist<-colSums(sqd)^(1/p)
  }
  
  # Mahalanobis
  # Covariance should be of the all dataset
  if (metricf=="Mahalanobis"){
    C= cov(dataset)
    dist <- sqrt( t(element1-element2) %*% solve(C) %*% (element1-element2) )
  }
  
  # Cosine distance
  if (metricf=="Cosine"){
    dist <- sum(element1*element2)/sqrt(sum(element1^2)*sum(element2^2))
  }
  
  # Canberra distance
  if (metricf=="Canberra"){
    for(i in seq(along=element1)){
      sqd[i]<-(abs(element1[i]-element2[i])) / (abs(element1[i]) + abs(element2[i]))
    }
    dist<-colSums(sqd)
  }
  
  return(dist)
}
