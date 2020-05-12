library(tseries)
library(urca)
library(vars)
library(dynlm)
library(fpp2)
library(quantmod)
library(dplyr)
library(matlib)
library(lattice)
library(shrink)

setwd("~/MA 810")
studentdata <- read.csv("Student_Data.csv", header=TRUE)
studentdataWOnames <- read.csv("Student_Data_woNames.csv", header=TRUE)
pastdata <- subset(studentdataWOnames, select= -c(Alg2, Geom, Alg2Proj))
alg2score <- studentdataWOnames[, "Alg2"]
alg2projREAL <- studentdataWOnames[,"Alg2Proj"]


alg2proj <- numeric(nrow(pastdata))
alg1score <- studentdataWOnames[, "Alg1"]
pd_noalg1 <- subset(pastdata, select= -c(Alg1))
error <- numeric(length(alg2proj))
percenterror <- numeric(length(alg2proj))
score17 <- studentdataWOnames[, "y2017"]
score15 <- studentdataWOnames[, "y2015"]
score14 <- studentdataWOnames[, "y2014"]
score13 <- studentdataWOnames[, "y2013"]

#betas <- matrix(nrow =nrow(dummydata), ncol = ncol(dummydata))
#difference <- matrix(nrow =nrow(dummydata), ncol = ncol(dummydata))

betas <- matrix(nrow =nrow(pastdata), ncol = ncol(pastdata))
difference <- matrix(nrow =nrow(pastdata), ncol = ncol(pastdata))

for (i in 1:nrow(pastdata)){
  dummydata <- pastdata
  for (j in 1:ncol(pastdata)){
    #Remove missing data
    if (is.na(pastdata[i,j])){
      dummydata <- subset(dummydata, select=-c(j))
      }
  }
  #Get Betas for a student
  beta <- numeric(ncol(dummydata))
  covxx <- cov(dummydata, use = "complete.obs")
  covxy <- cov(dummydata, alg2score, use = "complete.obs")
  covxx_inv <- solve(covxx)
  beta <- covxx_inv %*% covxy
  sumbeta <- sum(beta)
  for (m in 1:length(beta)){
    beta[m] <- beta[m]/sumbeta
    betas[i,m] <- beta[m]
  }
  #Start process of finding student estimates
  total = studentdata[1,"Alg2"]
  for (k in 1:ncol(dummydata)){
    #create individual entries and sum them up
    if (is.na(dummydata[i,k])){
      difference[i,k] = 0
    } else{
      difference[i,k] = beta[k]*(dummydata[i,k] - dummydata[1,k])
    }
    
    total = total + difference[i,k]
  }
  alg2proj[i] <- total
  error[i] <- alg2proj[i] - studentdata[i, "Alg2Proj"]
  percenterror[i] <- abs(alg2proj[i] - studentdata[i, "Alg2Proj"])/studentdata[i,"Alg2Proj"]
}

meanerror <- mean(error)
meanpe <- mean(percenterror)

##### IGNORE BELOW FOR NOW
#Weights done a second time with Rough Algebra II projections as response variable.

alg2projNEW <- numeric(nrow(pastdata))

for (i in 1:nrow(pastdata)){
  dummydata <- pastdata
  for (j in 1:ncol(pastdata)){
    #Remove missing data
    if (is.na(pastdata[i,j])){
      dummydata <- subset(dummydata, select=-c(i))
    }
  }
  #Get Betas for a student
  beta <- numeric(ncol(dummydata))
  covxx <- cov(dummydata, use = "complete.obs")
  covxy <- cov(dummydata, alg2proj, use = "complete.obs")
  covxx_inv <- solve(covxx)
  beta <- covxx_inv %*% covxy
  sumbeta <- sum(beta)
  for (m in 1:length(beta)){
    #beta[m] <- beta[m]/sumbeta
    betas[i,m] <- beta[m]
  }
  #Start processof finding student estimates
  total <- studentdata[1,"Alg2"]
  for (k in 1:ncol(dummydata)){
    #create individual entries and sum them up
    if (is.na(dummydata[i,k])){
      difference[i,k] = 0
    } else{
      difference[i,k] = beta[k]*(dummydata[i,k] - dummydata[1,k])
    }
    
    total = total + difference[i,k]
  }
  alg2projNEW[i] <- total
  error[i] <- alg2projNEW[i] - studentdata[i, "Alg2Proj"]
  percenterror[i] <- abs(alg2projNEW[i] - studentdata[i, "Alg2Proj"])/studentdata[i,"Alg2Proj"]
}

meanerror <- mean(error)
meanpe <- mean(percenterror)
########################################################


#Teacher Model
#18 should have beat prediction
studentscore <- alg2score[-1]
#studentprojREAL <- alg2projREAL[-1]
studentproj <- alg2proj[-1]
lm(studentscore ~ studentproj) #not sure here...
modelpred <- 0.7997*studentproj + 70.4354
growth <- studentscore - modelpred
sumgrowth <- sum(growth)
stderror <- sd(growth)/sqrt(length(growth))
index <- sum(growth)/stderror


easydiff <- studentscore - studentproj
easyerror <- sd(easydiff)/sqrt(length(easydiff))
easyindex <- sum(easydiff)/easyerror


#Linear Regression with data?
summary(aov(alg2score ~ alg1score + score17 + score15 + score14 + score13))
summary(lm(alg2score ~ alg1score + score17 + score15 + score14 + score13))


#Use Ancova Model with data that already exists?