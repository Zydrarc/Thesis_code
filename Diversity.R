library(vegan)
library("dplyr")
library(viridis)

options(scipen=0) # scientific notation on
options(scipen=999) # scientific notation off (this is needed as the loop would not be calling the right file with out this)


experiments = 10 # number of experiments run
tstep= 100000
tsave= 5000 # every how many timesteps the experiment is saved 
timesteps = seq(0,tstep, by = tsave)
tsnum = ((tstep/tsave)+1)

wd<-getwd()
divres<-data.frame()


for (x in 1:experiments)
{
  for (i in timesteps)
  {
    #read data from the results of experiment
    loc<- paste(wd,"/experiment_",x,"/data/detail-",i,".spop",sep = "")
    data<- read.table(file = loc,skip =25, sep = " ", dec = ".", header = FALSE,skipNul = TRUE, fill = TRUE )
    
    #Calculate Shannon's diversity index
  num <-data$V5
  D= diversity(num,"invsimpson")
  shan=diversity(num,"shannon")
  ind = sum(data$V5)
  
  newres<-c(x,i,D,shan,ind)
  divres<- rbind(divres,newres)
  print(paste("experiment",x,"timestep",i,"complete"))
  }
  print(paste("experiment",x,"complete"))
}


colnames(divres)=c("Res","ts","D","Shan","indidivuals")

hue <- viridis(experiments, alpha = 1, begin = 0, end = 1, direction = 1,option = "viridis")

plot(divres$ts[divres$Res==1],divres$D[divres$Res==1],type="l",xlab = "Timesteps", ylab = "Simpson's Reciprocal Index", ylim = c(0,3500),col = hue[1],lwd=2.0)
for (n in 2:experiments){
  lines(divres$ts[divres$Res==n],divres$D[divres$Res==n],type="l",col= hue[n],lwd=2.0)
}
legend("bottomright", legend=c(1:10),col=hue[1:10], lty=1, lwd = 2.0)

plot(divres$ts[divres$Res==1],divres$Shan[divres$Res==1],type="l",xlab = "Timesteps", ylab = "Shannon's Index",col= hue[n],lwd=2.0,ylim = c(0,10))
for (n in 2:experiments){
  lines(divres$ts[divres$Res==n],divres$Shan[divres$Res==n],type="l",col= hue[n],lwd=2.0)
}
legend("bottomright", legend=c(1:10),col=hue[1:10], lty=1, lwd = 2.0)

## plot for number of individuals

plot(divres$ts[divres$Res==1],divres$indidivuals[divres$Res==1],type="l",xlab = "Timesteps", ylab = "Number of Individuals",xlim = c(0,100000),col = hue[1],lwd=2.0)
for (n in 2:experiments){
  lines(divres$ts[divres$Res==n],divres$indidivuals[divres$Res==n],type="l",col= hue[n],lwd=2.0)
}
legend("bottomright", legend=c(1:10),col=hue[1:10], lty=1, lwd = 2.0)



#using the phenotype data

phendata<-data.frame()

for (ex in 1:experiments){
  
  #read data from the results of experiment
  pc<- paste(wd,"/experiment_",ex,"/data/phenotype_count.dat",sep = "")
  phen<- read.table(file = pc,skip =8, sep = " ", dec = ".", header = FALSE,skipNul = TRUE, fill = TRUE )
  
for (i in 1:tsnum){
  newphen <- c(ex,phen$V1[i],phen$V3[i],phen$V4[i],phen$V5[i])
  phendata<-rbind(phendata,newphen)
}
  
}

colnames(phendata)=c("Exp","ts","Shan","phen","aD")

plot(phendata$ts[phendata$Exp==1],phendata$Shan[phendata$Exp==1],type="l",xlab = "Timesteps", ylab = "Shanon's Index", ylim = c(0,5),col = hue[1],lwd=2.0)
for (n in 2:experiments){
  lines(phendata$ts[phendata$Exp==n],phendata$Shan[phendata$Exp==n],type="l",col= hue[n],lwd=2.0)
}
legend("bottomright", legend=c(1:10),col=hue[1:10], lty=1, lwd = 2.0)

plot(phendata$ts[phendata$Exp==1],phendata$phen[phendata$Exp==1],type="l",xlab = "Timesteps", ylab = "Unique Phenotypes", ylim = c(0,2000),col = hue[1],lwd=2.0)
for (n in 2:experiments){
  lines(phendata$ts[phendata$Exp==n],phendata$phen[phendata$Exp==n],type="l",col= hue[n],lwd=2.0)
}
legend("topleft", legend=c(1:10),col=hue[1:10], lty=1, lwd = 2.0)

plot(phendata$ts[phendata$Exp==1],phendata$aD[phendata$Exp==1],type="l",xlab = "Timesteps", ylab = "Average Task Diversity", ylim = c(0,5),col = hue[1],lwd=2.0)
for (n in 2:experiments){
  lines(phendata$ts[phendata$Exp==n],phendata$aD[phendata$Exp==n],type="l",col= hue[n],lwd=2.0)
}
legend("topleft", legend=c(1:10),col=hue[1:10], lty=1, lwd = 2.0)



#regression model

phenmodel<-data.frame()

for (ex in 2:experiments){
  
  #read data from the results of experiment
  pc<- paste(wd,"/experiment_",ex,"/data/phenotype_count.dat",sep = "")
  phen<- read.table(file = pc,skip =28, sep = " ", dec = ".", header = FALSE,skipNul = TRUE, fill = TRUE )
  
  for (i in 1){
    newphen <- c(ex,phen$V1[i],phen$V3[i],phen$V4[i],phen$V5[i])
    phenmodel<-rbind(phenmodel,newphen)
  }
  
}
colnames(phenmodel)=c("Exp","ts","Shan","phen","aD")

model <-lm(phenmodel$Shan ~ phenmodel$Exp, data = phenmodel)
summary(model)
plot(model$residuals)

model2 <-lm(phenmodel$Shan ~ I(phenmodel$Exp^2), data = phenmodel)
summary(model2)
plot(model2$residuals)
