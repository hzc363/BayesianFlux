source("Rfunctions.r")

# simulate data with SY/b = 0.5
S = fluxSimu(jump=1,sdY = 0.5) 

# fit the data with the bayesian model
# Returns the posterior distribution of parameters described in the paper.
results = Bflux(S$Ca2,S$Time) 

# plot the result from Bflux, returns the probability of activation at each time point
POA = fluxPlot(results,S$Ca2,S$Time)  
