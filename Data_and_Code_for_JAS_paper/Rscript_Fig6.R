library(depmixS4)
library(dplyr)
library(wavethresh)
# Read data----------
data = read.csv('SIINFEKL data reshaped.csv',check.name = F)
data = as.matrix(data)
p_active = matrix(0,nrow(data),60)


# Fit model (Figure A----------------
for(CellID in 1:nrow(data)){
  #CellID = 3 or 8
  goodData = (data[CellID,]>0)
  Y = data[CellID,goodData]
  time = (1:ncol(data))[goodData]
  data_N = length(Y)
  
 
  results = Bflux(Y,time)
  #plotting#########
  Act = unlist(results$pos_Act)
  Act_p = table(Act)/iter
  Act_n = as.integer(names(Act_p))
  Act_n = Act_n[Act_n<=60&Act_n>0]
  p_active[CellID,Act_n]=Act_p[as.character(Act_n)]
  
  par(mfrow=c(2,2))
  plot(time,Y,ylim=c(0,max(Y)+1),main="Data",xlab="Time",ylab="Ca2+")
  plot(time,Y,ylim=c(0,max(Y)+1),main="Fitted Data",xlab="Time",ylab="Ca2+")
  fitData=apply(results$pos_Y, 2,mean)
  points(1:length(fitData),fitData,col="red")
  #hist(Act,xlim=c(0,time[data_N]),main="Activation Time Points",freq=T,xlab="Time",breaks = time[data_N])
  plot(Act_p, xlim=c(0,time[data_N]),main="Prob of being activated",xlab="Time",ylab='probability')
  ylim = c((min(results$pos_NAct)-1),(max(results$pos_NAct)+1))
  plot(results$pos_NAct,main="Number of Activation",xlab="Iteration", ylab="",ylim = ylim ,yaxt = "n")
  axis(2, at = ylim[1]:ylim[2])
}

write.csv(p_active,'POA.csv')

# plot inital activation-------------
prob = read.csv('POA.csv')[-1]
max_prob=apply(prob,1,which.max)

initial = function(x,a){
  for(i in 1:length(x)){
    if(x[i]>=a){break}
  }
  return(i)
}

j=1
first = data.frame(IT=rep(NA,200),cutoff=rep(NA,200),cellID=rep(NA,200))
for(cutoff in seq(0.1,0.9,0.1)){
  for(k in 1:nrow(prob)){
    first[j,"IT"] = initial(prob[k,],cutoff)
    first[j,'cutoff']=cutoff
    first[j,'cellID']=k
    j = j+1
  }
}

first = first[-which.max(first[,"IT"]),]
first$IT=first$IT/4
boxplot(first[,"IT"]~first[,'cutoff'],xlab="Cutoff", ylab="Estimated Initial Activation Time",ylim=c(1,7))
abline(15/4,0,col='black',lty="dashed",lwd=2)
sd(first[first[,'cutoff']==0.5,'IT'],na.rm=T)

#smooth with wavelet and inference with HMM----------
library("wavethresh")
p_active=NULL
for(CellID in 1:nrow(data)){
  
  Y = data[CellID,]
  Y[is.na(Y)]=0
  w=which(Y>0)
  for(i in 1:length(Y)){
    if(Y[i]<=0){
      near_id= w[order(abs(w-i))[1:2]]
      Y[i]=mean(Y[near_id])
    }
  }
  time = (1:ncol(data))
  #plot(time, Y,ylim=c(0,2))
  
  Y2=c(Y, rep(0,2^(ceiling(log2(length(Y))))-length(Y))  )
  wy <- wd(Y2)
  thresh <- threshold(wy)
  yr <- wr(thresh)[1:length(Y)]
  #lines(time, yr)
  
  mod <- depmix(Y ~ 1, family = gaussian(), nstates = 3, data = data.frame("Y"=yr))
  set.seed(1)
  fm2 <- fit(mod, verbose = FALSE)
  probs <- posterior(fm2) 
  
  plot(time,yr,col=probs[,1])
  t1=data.frame("Y"=yr,"state"=probs$state)
  t1=t1%>%group_by(state)%>%summarise(mean(Y))
  #t1=which.max(t1$`mean(Y)`)
  t1= order(t1$`mean(Y)`)[2]
  p_active=rbind(p_active,probs[,(t1+1)])
}
write.csv(p_active,"POA_HMM.csv")

# plot inital activation estimated by HMM and wavelet-------------
prob = read.csv('POA_HMM.csv')[-1]
max_prob=apply(prob,1,which.max)

initial = function(x,a){
  for(i in 1:length(x)){
    if(x[i]>=a){break}
  }
  return(i)
}

j=1
first = data.frame(IT=rep(NA,200),cutoff=rep(NA,200),cellID=rep(NA,200))
for(cutoff in seq(0.1,0.9,0.1)){
  for(k in 1:nrow(prob)){
    first[j,"IT"] = initial(prob[k,],cutoff)
    first[j,'cutoff']=cutoff
    first[j,'cellID']=k
    j = j+1
  }
}

first = first[-which.max(first[,"IT"]),]
first$IT=first$IT/4
boxplot(first[,"IT"]~first[,'cutoff'],xlab="Cutoff", ylab="Estimated Initial Activation Time",ylim=c(1,7))
abline(15/4,0,col='black',lty="dashed",lwd=2)
sd(first[first[,'cutoff']==0.5,'IT'],na.rm=T)


