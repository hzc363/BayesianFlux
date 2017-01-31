# Read the data (OVA- data: WT_indo.csv, OVA+ data: OVAhi_indo.csv)
data = read.csv('WT_indo.csv',check.name = F)
data = as.matrix(data)

# Fit the data and calculate the probability of activation 
P_act = matrix(NA,dim(data)[1],dim(data)[2])
ID_seq = 1:nrow(data)
for(CellID in ID_seq){
  goodData = (data[CellID,]>0 & is.na(data[CellID,])==F )
  Y = data[CellID,goodData]
  time = (1:ncol(data))[goodData]
  End = max(time)
  Begin = min(time)
  time_c = Begin:End
  
  results = Bflux(Y,time,BI=10000,iter=5000)
  p = fluxPlot(results,Y,time,plot=F)
  P_act[CellID,time_c]=p[as.character(time_c)]
  print(paste(CellID,'of', dim(data)[1]))
}

write.csv(P_act,'POA.csv')


# Derive descriptive statistics 
poa= read.csv('POA.csv')[,-1]
poa= as.matrix(poa)
data = read.csv('WT_indo.csv')
par(mfrow=c(2,2))

NT_active = sum(poa[is.finite(poa)]>0.5,na.rm=T) #total number of activation time points
NT_total = sum(is.finite(poa),na.rm=T) #total number of time points
PC_active = NT_active/NT_total*100 # percent of time in activation

active_count = 0 #total number of activation events
active_count_per_cell = rep(0, nrow(poa)) #number of activation of each cell
length_active = c() #length of each activation event in min
n=1

for(i in 1:nrow(poa)){
  
  k=0
  good = !is.na(poa[i,])
  if(prod(poa[i,good]>0.5)==1){next} 
  NT_active = NT_active + sum(poa[i,good]>0.5)
  for(j in 2:(ncol(poa)-1)){
    
    if(is.na(poa[i,j])){next}
    if(is.na(poa[i,j-1])){next}
    if(is.na(poa[i,j+1])){next}
    if(poa[i,j-1]<0.5 & poa[i,j]>0.5){
      active_count = active_count+1
      active_count_per_cell[i]=active_count_per_cell[i]+1
      k = j
    }
    if(poa[i,j]>0.5 & poa[i,j+1]<0.5&k!=0){length_active[n]=j-k+1; n = n+1;k=0}
  }
}

freq_active = active_count/nrow(data)*4 #number of activation per hour

hist(active_count_per_cell+0.1,breaks = 0:4,freq=F,ylim = c(0,1)) # histagram of number of activation of each cell
hist(length_active/4,freq=F,ylim = c(0,1),breaks = 0:7) # histagram of the length of each flux

