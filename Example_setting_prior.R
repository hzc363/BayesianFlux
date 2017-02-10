source("Rfunctions.r")

# simulate 3 data, all resting, all activating, and mixed
S_active = list("Ca2"=rnorm(60,mean=3,sd=0.5),"Time"=1:60) 
S_rest = list("Ca2"=rnorm(60,mean=1,sd=0.5),"Time"=1:60) 
S_mix = fluxSimu(jump=1,sdY = 0.5) 


# fit the data with the bayesian model
# Returns the posterior distribution of parameters described in the paper.
# Here, we beleve that the baseline should be somewhere 1, so we specify 
# Ea to be 1. And use a small standard deviation Sa=0.2 to indicate that we are 
# pretty sure the baseline should be around 1. 

results_active = Bflux(S_active$Ca2,S_active$Time,Ea=1,Sa=0.2) 
results_rest = Bflux(S_rest$Ca2,S_rest$Time,Ea=1,Sa=0.2) 
results_mix = Bflux(S_mix$Ca2,S_mix$Time,Ea=1,Sa=0.2) 


# plot the result from Bflux, returns the probability of activation at each time point
# With the right "prior knowledge", the program can tell "all-activation" from
# "all-resting". The Probability of Activation for the S_active are all above 0.5
# while most of the probability of activatio for the S_rest are lower than 0.5.
POA_active = fluxPlot(results_active,S_active$Ca2,S_active$Time)  
POA_rest = fluxPlot(results_rest,S_rest$Ca2,S_rest$Time)  
POA_mix = fluxPlot(results_mix,S_mix$Ca2,S_mix$Time)  

