### Figure 4 A-C
S = fluxSimu(jump=1,sdY = 0.5) # simulate data with SY/b = 0.5
results = Bflux(S$Ca2,S$Time) # fit the data with the bayesian model
p = fluxPlot(results,S$Ca2,S$Time) # plot the result from Bflu, returns the probability of activation 

### Figure 4D
library(ROCR) #requires package (ROCR), available in CRAN
ID_a = which(S$Expect==max(S$Expect))
End = max(S$Time);Begin = min(S$Time)
true_ST = rep(0,End)
true_ST[ID_a]=1

pred = prediction(as.vector(p[Begin:End]), true_ST)
roc = performance(pred, measure = 'tpr', x.measure = 'fpr') #Generating the ROC curve
plot(roc,lwd=3,main='ROC curve')


### Figure 4E (takes time to run)
auc_result = data.frame(sdY=rep(0,750),auc=rep(0,750))
NSR_seq = seq(0.2,3,0.2)
index = 1
for(sdY in NSR_seq){
  print(paste('simulating with NSR=',sdY))
  for(i in 1:50){
    S = fluxSimu(jump=1,sdY = sdY)
    results = Bflux(S$Ca2,S$Time)
    p = fluxPlot(results,S$Ca2,S$Time,plot=F)
    
    ID_a = which(S$Expect==max(S$Expect))
    End = max(S$Time);Begin = min(S$Time)
    true_ST = rep(0,End)
    true_ST[ID_a]=1
    pred = prediction(as.vector(p[Begin:End]), true_ST)
    auc = performance(pred, measure='auc') 
    
    auc_result[index,2]= auc@y.values[[1]]
    auc_result[index,1]= sdY
    index = index+1
  }
  
}

write.csv(auc_result, "acu_result.csv") # store AUC result, can be used in the future.
auc_result2=subset(auc_result, auc_result[,1]>0)
boxplot(auc_result2[,2]~auc_result2[,1],xlab="Noise to Signal Ratio",ylab="AUC of the ROC Curve")

