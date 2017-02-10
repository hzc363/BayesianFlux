########Functions for Simulating data########
## Returns a list contaning simulated calcium concentration (Ca2), mean of calcium concentration (Expect), and time (Time)
fluxSimu=function(Nt = 60, base = 1.5,jump=1,slope=0.5,sdY = 0.5,
                  ST1 = 0, transition = c(0,20,30,50,60)){
  state = rep(c(0,1),10)
  state = ST1 * (1-state) + (1-ST1)*state
  state = state[1:(length(transition)-1)]
  
  Ca2= c()
  Expect = c()
  
  for(x in 1:Nt){
    expect = base
    for(i in 1:(length(transition)-1)){
      inc = jump* (x>transition[i])*(x<=transition[i+1]) * ((state[i]==1) + (state[i]==0)*exp(-slope*(x-transition[i])))
      expect = expect + inc
    }
    Ca2[x] = rnorm(1, expect, sdY)
    Expect[x] = expect
  }
  results = list(Ca2=Ca2,Expect =Expect,Time=1:Nt )
}

##### Functions for data fitting #########
#calculate Ndt from dt
dt2Ndt_func = function(dt, End){
  sum = 0
  for(i in 2:length(dt)){
    sum = sum + dt[i]
    if(sum>=End){break}
  }
  return(i)
}

#convert from dt to Tt
dt2Tt_func = function(dt,Ndt){
  Tt = rep(0,Ndt)
  Tt[1] = -dt[1]
  Tt[-1]=cumsum(dt[2:Ndt])
  return(Tt)
}

#calculate Id (Id is a middle step in calculating the expectation)
Id_func = function(ST,Tt,c,Ndt,data_N,time){
  Id = rep(0,length(time))
  j=1
  for(t in time){
    NL = which.max( Tt[ Tt<t] )
    Id[j] = ((ST[NL]==1) + (ST[NL]==0)*exp(-c*(t-Tt[NL])))
    j=j+1
  }
  return(Id)
}

#Calculate the expectation of Y given parameters
EY_func = function(a,b,c,ST,Tt,Ndt,data_N,time){
  Id = Id_func(ST,Tt,c,Ndt,data_N,time)
  EY = a + b*Id
  return(EY)
}

#calculate the likelyhood given Y and parameters
LLH_func = function(a,b,c,ST,Tt,SY,Ndt,data_N,time,Y){
  EY = EY_func(a,b,c,ST,Tt,Ndt,data_N,time)
  LLH = sum(dnorm(Y,EY,SY,log=T))
  return(LLH)
}

#identify times in which cells are activated from Tt, ST
Tt2Act=function(Tt, ST,Ndt){
  Act = c()
  for(i in 1:(Ndt-1)){
    if(ST[i]==1){Act = c(Act, round(Tt[i]):round(Tt[i+1]))}
  }
  return(Act)
}

## The main function that runs the MCMC, returns a list containing samples from the posterior distribution of 
# parameters: c (pos_c), a (pos_a), b (pos_b), Y (pos_Y), Number of activation states (pos_NAct) , time point in which activation is assigned(pos_Act)
# SY (pos_SY)
Bflux = function(Y,time,
                 #Settings
                 iter = 800,BI=200,Sw = 0.1,
                 #hyper parameters
                 Ea = 0, Sa = 1000,Eb = 0, Sb = 1000,
                 ac = 0.0001, bc = 0.0009,ld = 30,
                 aS = 0.0001, bS = 0.0009, 
                 #Initial values
                 a = 2, b = 2, c = 2,ST=rep(c(0,1),10),
                 dt = rep(15, 20),SY = 1){
  
  data_N = length(Y)
  End = max(time)
  Ndt = dt2Ndt_func(dt, End) #number of dt that matters
  Tt = dt2Tt_func(dt,Ndt)
  K = length(dt)
  seqEY = seq(0,End,1)
  
  #record posterior sample
  pos_c = rep(NA,iter)
  pos_a = rep(NA,iter)
  pos_b = rep(NA,iter)
  pos_Y = matrix(NA,iter,length(seqEY))
  pos_NAct = rep(NA,iter)
  pos_Act = vector("list", iter)
  pos_SY = rep(NA,iter)
  
  #Run MCMC######
  for(i in 1:(iter+BI)){
    
    #Sample c 
    c_new = rnorm(1,c,Sw)
    lr = LLH_func(a,b,c_new,ST,Tt,SY,Ndt,data_N,time,Y)+dgamma(c_new,ac,bc,log = T)-LLH_func(a,b,c,ST,Tt,SY,Ndt,data_N,time,Y)-dgamma(c,ac,bc,log = T)
    r = exp(lr)
    if(runif(1)<r){c = c_new}
    if(i>BI){pos_c[i-BI]= c}
    
    #Sample ST 
    lp_new = LLH_func(a,b,c,(1-ST),Tt,SY,Ndt,data_N,time,Y)
    lp_old = LLH_func(a,b,c,ST,Tt,SY,Ndt,data_N,time,Y)
    rp = 1 + exp(lp_old - lp_new)
    p = 1/rp
    if(runif(1)<p){ST = 1-ST}
    
    #Sample dt[k] 
    dtj_new = rexp(1,ld)
    dt_new = dt; dt_new[1] = dtj_new; 
    Tt_new = dt2Tt_func(dt_new,Ndt)
    lr1= LLH_func(a,b,c,ST,Tt_new,SY,Ndt,data_N,time,Y)+dexp(dtj_new,ld,log = T)
    lr2 = LLH_func(a,b,c,ST,Tt,SY,Ndt,data_N,time,Y)+dexp(dt[1],ld,log=T)
    lr = lr1-lr2
    r = exp(lr)
    if(runif(1)<r){dt= dt_new; Tt = Tt_new}
    
    
    for(j in 2:Ndt){
      dtj_new = runif(1,0,(dt[j]+dt[j+1]))
      dtj1_new = dt[j]+dt[j+1] - dtj_new
      dt_new = dt; dt_new[j] = dtj_new; dt_new[j+1] = dtj1_new
      Ndt_new = dt2Ndt_func(dt_new, End)
      Tt_new = dt2Tt_func(dt_new,Ndt_new)
      lr1= LLH_func(a,b,c,ST,Tt_new,SY,Ndt_new,data_N,time,Y)+dexp(dtj_new,ld,log = T)+dexp(dtj1_new,ld,log = T)
      lr2 = LLH_func(a,b,c,ST,Tt,SY,Ndt,data_N,time,Y)+dexp(dt[j],ld,log=T)+dexp(dt[j+1],ld,log=T)
      lr = lr1-lr2
      r = exp(lr)
      if(runif(1)<r){dt= dt_new; Tt = Tt_new; Ndt=Ndt_new}
    }
    if(i>BI){pos_Act[[i-BI]]=Tt2Act(Tt,ST,Ndt)}
    if(i>BI){pos_NAct[i-BI] = sum(ST[1:(Ndt-1)])}
    
    #Sample a
    Pa = data_N/SY^2 + 1/Sa^2
    Id = Id_func(ST,Tt,c,Ndt,data_N,time)
    SUMa = sum(Y)-b*sum(Id)
    a_mean = (SUMa/SY^2+Ea/Sa^2)/Pa
    a_sd = Pa^(-1/2)
    a_new = rnorm(1,a_mean,a_sd)
    if(a_new >0){a = a_new}
    if(i>BI){pos_a[i-BI]= a}
    
    #Sample b 
    Pb = sum(Id^2)/SY^2+1/Sb^2
    SUMb = sum( (Y-a)*Id  )
    b_mean = (SUMb/SY^2+Eb/Sb^2)/Pb
    b_sd = Pb^(-1/2)
    b_new = rnorm(1,b_mean, b_sd)
    if(b_new >0){b = b_new} 
    if(i>BI){pos_b[i-BI]= b}
    
    #Sample SY 
    EY = EY_func(a,b,c,ST,Tt,Ndt,data_N,time)
    SSE = t(Y-EY)%*%(Y-EY)
    PS = rgamma(1,(data_N/2+aS), (SSE/2+bS) )
    SY = PS^(-1/2)
    if(i>BI){pos_SY[i-BI] = SY}
    
    EY = EY_func(a,b,c,ST,Tt,Ndt,data_N,seqEY)
    #calculate expected Y
    if(i>BI){pos_Y[(i-BI),] = EY}
  }
  
  results = list(pos_c=pos_c, pos_a=pos_a, pos_b=pos_b, 
                 pos_Y=pos_Y, pos_NAct=pos_NAct,pos_Act=pos_Act,pos_SY=pos_SY)
  return(results)
}



#Make plots. take the list returned by the Bflux function as the first argument
#returns the probability of activatio at each time point

fluxPlot = function(results,Y,time,plot=T){
  
  iter = length(results$pos_b)
  End = max(time)
  Begin = min(time)
  Act = c(unlist(results$pos_Act), Begin:End)
  Act_p = (table(Act)-1)/iter
  seqEY = seq(0,End,1)
  time_c = Begin:End
  
  #png(file = paste(CellID,".png"))
  if(plot == F){return(Act_p)}
  par(mfrow=c(2,2))
  plot(time,Y,ylim=c(min(Y)*0.8,max(Y)*1.2),main="Data",xlab="Time",ylab="Ca2+")
  plot(time,Y,ylim=c(min(Y)*0.8,max(Y)*1.2),main="Fitted Data",xlab="Time",ylab="Ca2+")
  points(seqEY,apply(results$pos_Y, 2,mean),col="black",pch=19)
  barplot(Act_p[as.character(time_c) ],ylim =c(0,1),names.arg=time_c,xlab='Time',ylab='Probability',main="Probability of Activation")
  ylim = c((min(results$pos_NAct)-1),(max(results$pos_NAct)+1))
  plot(results$pos_NAct,main="Number of Activation",xlab="Iteration", ylab="",ylim = ylim ,yaxt = "n")
  axis(2, at = ylim[1]:ylim[2])
  
  
  return(Act_p)
}

