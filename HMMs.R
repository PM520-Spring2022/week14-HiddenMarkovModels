# Based on example taken from http://web.stanford.edu/class/stats366/hmmR2.html
#install.packages('HMM')
library('HMM') 
set.seed(985)
require(HMM)
nSim = 2000
States = c("Fair", "Unfair")
Symbols = 1:6
# Define the transition matrix
transProbs = matrix(c(0.99, 0.01, 0.02, 0.98), c(length(States),length(States)), byrow = TRUE)
# Define the emission probabilities
emissionProbs = matrix(c(rep(1/6, 6), c(rep(0.1, 5), 0.5)),c(length(States), length(Symbols)), byrow = TRUE)
#Set up the HMM
hmm = initHMM(States, Symbols, transProbs = transProbs, emissionProbs = emissionProbs)
# Simulate some test data to play with
sim = simHMM(hmm, nSim)
# the resulting die-rolls are stored in sim$observation
plot(sim$observation,pch='.')
# the sequence recording which die was rolled is stored in sim$states. Let's look at the first few values
sim$states[1:20]

# find the most probable sequence of underlying states
vit = viterbi(hmm, sim$observation)   # the first argument is the name of your HMM; the second is the sequence of observed states you want to analyze

# Compute the (log) forward probabilities. The forward probability for state X up to observation at time k is defined as the probability of observing the 
# sequence of observations e_1, ... ,e_k and that the state at time k is X. That is:
# f[X,k] := Prob(E_1 = e_1, ... , E_k = e_k , X_k = X).
# Where E_1...E_n = e_1...e_n is the sequence of observed emissions and X_k is a random variable that represents the state at time k.
f = forward(hmm, sim$observation)

# Compute the (log) backward probabilities. The backward probability for state X and observation at time k is defined as the probability of observing the 
# sequence of observations e_k+1, ... ,e_n under the condition that the state at time k is X. That is:
# b[X,k] := Prob(E_k+1 = e_k+1, ... , E_n = e_n | X_k = X).
# Where E_1...E_n = e_1...e_n is the sequence of observed emissions and X_k is a random variable that represents the state at time k.
b = backward(hmm, sim$observation)

# calculate posterior probabilities of being in underlying states
i <- f[1, nSim]
j <- f[2, nSim]
probObservations = (i + log(1 + exp(j - i)))
posterior = exp((f + b) - probObservations)

# that showed you how to do it 'long-hand'.
# The same results can be obtained directly using the "posterior()" function:
direct_posterior <- posterior(hmm,sim$observation)

# combine all the results into a single list
x = list(hmm = hmm, sim = sim, vit = vit, posterior = posterior)

# create some labels for the plots
mn = "Fair and unfair die"
xlb = "Throw nr."
ylb = ""

# plot the sequence of observed states
plot(x$sim$observation, ylim = c(-7.5, 6), pch = 3, main = mn,
     xlab = xlb, ylab = ylb, bty = "n", yaxt = "n")
axis(2, at = 1:6)

# plot the dice that was used at each throw - we use green to represent the fair dice
text(0, -1.2, adj = 0, cex = 0.8, col = "black", "True: green = fair die")
for (i in 1:nSim) {
  if (x$sim$states[i] == "Fair")
    rect(i, -1, i + 1, 0, col = "green", border = NA)
  else rect(i, -1, i + 1, 0, col = "red", border = NA)
}

########Plot the most probable sequence of underlying true states (i.e. which die was used at each roll) calculated by the viterbi call#######################
text(0, -3.2, adj = 0, cex = 0.8, col = "black", "Most probable path")
for (i in 1:nSim) {
  if (x$vit[i] == "Fair")
    rect(i, -3, i + 1, -2, col = "green", border = NA)
  else rect(i, -3, i + 1, -2, col = "red", border = NA)
}


################## Show where the last two plots differ
text(0, -5.2, adj = 0, cex = 0.8, col = "black", "Difference")
differing = !(x$sim$states == x$vit)
for (i in 1:nSim) {
  if (differing[i]) # the states differ
    rect(i, -5, i + 1, -4, col = rgb(0.3, 0.3, 0.3),border = NA)
  else rect(i, -5, i + 1, -4, col = rgb(0.9, 0.9, 0.9),border = NA)
}

#################Plot the marginal posterior probability that we are using the unfair dice at each step#########################
points(x$posterior[2, ] - 3, type = "l")

###############Show when the posterior probs infer the wrong underlying state:############
text(0, -7.2, adj = 0, cex = 0.8, col = "black", "Difference by posterior-probability")
differing = !(x$sim$states == x$vit)
for (i in 1:nSim) {
  if (posterior[1, i] > 0.5) {
    if (x$sim$states[i] == "Fair")
      rect(i, -7, i + 1, -6, col = rgb(0.9, 0.9, 0.9),
           border = NA)
    else rect(i, -7, i + 1, -6, col = rgb(0.3, 0.3, 0.3),
              border = NA)
  }
  else {
    if (x$sim$states[i] == "Unfair")
      rect(i, -7, i + 1, -6, col = rgb(0.9, 0.9, 0.9),
           border = NA)
    else rect(i, -7, i + 1, -6, col = rgb(0.3, 0.3, 0.3),
              border = NA)
  }
}

# Baum-Welch
bw3 = baumWelch(hmm,sim$observation,maxIterations=100,pseudoCount=0)
# let's compare actual and estimated transition probabilities/ Start by vectorizing both matrices
Act<-c(hmm$transProbs)
Est<-c(bw3$hmm$transProbs)
plot(Act,Est)
abline(0,1,lty=2)
