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
library(ggplot2)

setwd("~/MA 810")

studentdata <- read.csv("Student_Data.csv", header=TRUE)
studentdataWOnames <- read.csv("Student_Data_woNames.csv", header=TRUE)
pastdata <- subset(studentdataWOnames, select= -c(Alg2, Geom, Alg2Proj))
pastdata2 <- subset(studentdataWOnames, select=-c(Alg2, Alg2Proj))
alg2score <- studentdataWOnames[, "Alg2"]
alg2projREAL <- studentdataWOnames[,"Alg2Proj"]


alg2proj <- numeric(nrow(pastdata2))
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

betas <- matrix(nrow =nrow(pd_noalg1), ncol = ncol(pd_noalg1))
difference <- matrix(nrow =nrow(pd_noalg1), ncol = ncol(pd_noalg1))
stu_z <- numeric(ncol(pd_noalg1))
alg2maybe <- numeric(ncol(pd_noalg1))

for (i in 1:nrow(pd_noalg1)){
  dummydata <- pd_noalg1
  for (j in 1:ncol(pd_noalg1)){
    #Remove missing data
    if (is.na(pd_noalg1[i,j])){
      dummydata <- subset(dummydata, select=-c(j))
    }
  }
  #Get Betas for a student
  beta <- numeric(ncol(dummydata))
  covxx <- cov(dummydata, use = "complete.obs")
  covxy <- cov(dummydata, alg1score, use = "complete.obs")
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
  stu_z[i] <- (total - studentdata[1, "Alg1"])/21.063
  alg2maybe[i] <- stu_z[i]*21.063 + studentdata[1,"Alg2"]
}

meanerror <- mean(error)
meanpe <- mean(percenterror)

##### IGNORE BELOW FOR NOW
#Weights done a second time with Rough Algebra II projections as response variable.

alg2projNEW <- numeric(nrow(pastdata))

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
  covxy <- cov(dummydata, alg2proj, use = "complete.obs")
  covxx_inv <- solve(covxx)
  beta <- covxx_inv %*% covxy
  sumbeta <- sum(beta)
  for (m in 1:length(beta)){
    beta[m] <- beta[m]/sumbeta
    betas[i,m] <- beta[m]
  }
  #Start process of finding student estimates
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

finaldata <- cbind(studentdataWOnames, error, percenterror)

#Teacher Model
#18 should have beat prediction
studentscore <- alg2score[-1]
#studentprojREAL <- alg2projREAL[-1]
studentproj <- alg2projREAL[-1]
lm(studentscore ~ studentproj) #not sure here...
modelpred <- 1.138*studentproj -38.559
growth <- studentscore - modelpred
avggrowth <- sum(growth)/length(growth)
stderror <- sd(growth)/sqrt(length(growth))
index <- avggrowth/stderror

easydiff <- studentscore - studentproj
easyerror <- sd(easydiff)/sqrt(length(easydiff))
easyindex <- sum(easydiff)/easyerror


#Linear Regression with data?
summary(aov(alg2score ~ alg1score + score17 + score15 + score14 + score13))
summary(lm(alg2score ~ alg1score + score17 + score15 + score14 + score13))
summary(lm(alg2score[-1] ~ alg1score[-1] + score17[-1] + score15[-1] + score14[-1] + score13[-1]))
summary(aov(lm(alg2score ~ alg1score + score17 + score15 + score14 + score13), projections=TRUE))

modelproj <- numeric(nrow(pastdata))
growthNEW <- numeric(nrow(pastdata))
#MODEL WITHOUT THE MEAN
for (i in 2:nrow(pastdata)){
  modelproj[i] <- alg1score[i]*0.516473 + score17[i]*-.03514 + score15[i]*.007344 + score14[i]*0.243914 + score13[i]*.001644 + 141.829526
  if (!is.na(modelproj[i])){
    growthNEW[i] <- studentscore[i-1] - modelproj[i]
  }
}

#MODEL WITH THE MEAN
for (i in 2:nrow(pastdata)){
  modelproj[i] <- alg1score[i]*0.5238 + score17[i]*-.005947 + score15[i]*.001335 + score14[i]*0.2615 + score13[i]*.0005849 + 136.4
  if (!is.na(modelproj[i])){
    growthNEW[i] <- studentscore[i-1] - modelproj[i]
  }
}
sumgrowth <- sum(growthNEW)
stderrorNEW <- sd(growthNEW)/sqrt(length(growthNEW))
index <- sum(growthNEW)/stderrorNEW

#Use Ancova Model with data that already exists?

graphx <- numeric(1000)
graphy <- numeric(1000)

for (i in 0:length(graphx)){
  pointx = rnorm(1, mean = 305, sd=21.063)
  pointy = rnorm(1, mean = pointx+2, sd=21.063)
  graphx[i] = pointx
  graphy[i] = pointy
}

#graphx <- round(rnorm(1000, mean=305, sd=21.063))
#graphy <- round(rnorm(1000, mean=305, sd=21.063))

fakelist <- data.frame(graphx, graphy)
fakelist.new <- data.frame(append(fakelist, c(x1='random_data')))
reallist <- data.frame(graphx=alg2projREAL, graphy=alg2score)
reallist.new <- data.frame(append(reallist, c(x1='real_data')))

totallist.new <- rbind(fakelist.new, reallist.new)

ggplot(reallist, aes(alg2projREAL, alg2score))+
  geom_point(color = "red", size=3) +
  geom_smooth(method='lm', color="red", se=FALSE, size=2)+
  geom_abline(slope=1, color="black", size=2)+
  xlab("Projected Score")+
  ylab("Actual Score")+
  ggtitle("Projected Versus Actual Score")

ggplot(totallist.new, aes(graphx, graphy))+
  geom_point(aes(graphx, graphy, color=x1), size=3) +
  scale_color_manual(values=c("plum","red"))+
  geom_smooth(method='lm', se=FALSE, size=2)+
  geom_abline(slope=1, color="black", size=2, alpha=0.2)+
  xlab("Projected Score")+
  ylab("Actual Score")+
  ggtitle("Projected Versus Actual Score")


#Regression based on random data
fit <- lm(totallist.new[,"graphy"]~totallist.new[,"graphx"])

predictedfit <- predict(fit)

growth <- alg2score - predictedfit[(length(predictedfit)-22):(length(predictedfit))]
avggrowth <- sum(growth)/length(growth)
stderror <- sd(growth)/sqrt(length(growth))
index <- avggrowth/stderror

#Average teacher score
(204*5+98*4+212*3+90*2+210*1)/(204+98+212+90+210)
