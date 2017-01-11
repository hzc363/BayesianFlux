source("Rfunctions.r")

S = fluxSimu(jump=1,sdY = 0.5) # simulate data with SY/b = 0.5

results = Bflux(S$Ca2,S$Time) # fit the data with the bayesian model

p = fluxPlot(results,S$Ca2,S$Time) # plot the result from Bflu, returns the probability of 