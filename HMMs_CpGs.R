
set.seed(98655)
require(HMM)

nSim = 5000
# we define a hidden layer that has 8 states:
# A+,C+,G+,T+     -these correspond to being in a CpG island
# A-,C-,G-,T-     -these correspond to being in a CpG island
States = c("A+","C+","G+","T+","A-","C-","G-","T-")
# but we only see the actual nucleotide, so our observed stats  (the "symbols") are as follows:
Symbols = c("A","C","G","T")

# we now need to build the transition matrix. We assume the following:
# nucleotide transition matrices for CpG islands and then non-CpG islands are as follows:
# (nucleotides are listed in order A,C,G,T.)
CpGProbs<-matrix(c(.180, .274, .426, .120, .171, .368, .274, .187, .161, .339, .375, .125, .079, .355, .384, .182),
                 nrow=4,ncol=4,byrow=TRUE)
NonCpGProbs<-matrix(c(.300,.205,.285,.210,.322,.298,.078,.302,.248,.246,.298,.208,.177,.239,.292,.292),
                    nrow=4,ncol=4,byrow=TRUE)
# We also need the transition probs for moving between CpG and non-CpG states
# we will assume that these are independent of the current nucleotide
ProbCpGToNonCpG<-0.01
ProbNonCpGToCpG<-0.01
# we now build the transition matrix for the hidden layer - start with the probs if we stay in (or out) of a CpG island
transProbs=matrix(nrow=8,ncol=8)
for (i in 1:4){
  for (j in 1:4){
    transProbs[i,j]<-(1-ProbCpGToNonCpG)*CpGProbs[i,j]
    transProbs[i+4,j+4]<-(1-ProbNonCpGToCpG)*NonCpGProbs[i,j]
  }
}
# Now define the probs when we transition between CpG and non-CpG (or vice versa)
# We will assume we choose the next nucleotide uniformly in that case
for (i in 1:4){
  for (j in 5:8){
    transProbs[i,j]<-ProbCpGToNonCpG/4
    transProbs[j,i]<-ProbNonCpGToCpG/4
  }
}
transProbs

# check we haven't messed up
rowSums(transProbs)

# Define the emission probabilities the prob that we emit an A,C,G, or T.
# We assume we emit the same value as the hidden state, but without the plus or minus
# the rows correspond to the hidden states; the columns correspond to the observed states
emissionProbs = matrix(nrow=8,ncol=4)
for (i in 1:8){
  for (j in 1:4){
    if ((i==j)||((i-4==j))){
      emissionProbs[i,j]<-0.997   # see below - I had to allow for errors here. Otherwise the model fit blew up
    }else
    {
      emissionProbs[i,j]<-0.001   # I got errors here if I set it to 0
    }
  }
}
emissionProbs

# we have now defined everything. So,...
#Set up the HMM
hmm = initHMM(States, Symbols, transProbs = transProbs, emissionProbs = emissionProbs)
# Simulate some test data to play with
sim = simHMM(hmm, nSim)
# the resulting nucleotides are stored in sim$observation

# the sequence recording the true underlying state is stored in sim$states. Look at the first few values
sim$states[1:20]
# Let's create a variable to record whether or not we are in a CpG Island
Island<-rep("FALSE",nSim)
for (i in 1:nSim){
  Island[i]<-((sim$states[i]=="A+")||(sim$states[i]=="C+")||(sim$states[i]=="G+")||(sim$states[i]=="T+"))  
}

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

# combine all the results into a single list
x = list(hmm = hmm, sim = sim, vit = vit, posterior = posterior)



# create some labels for the plots
mn = "CpGIsland or NoCpGIsland "
xlb = "Nucleotide nr."
ylb = ""

# plot the sequence of observed states
# create a variable that recodes the nucleotides as 1 through 4
obs<-as.numeric(factor(x$sim$observation))
#plot it
plot(obs, ylim = c(-7.5, 6), pch = 3, main = mn,
     xlab = xlb, ylab = ylb, bty = "n", yaxt = "n")
axis(2, at = 1:4)

# plot whether or not we were in a CpG island - we use green to represent being in an island
text(0, -1.2, adj = 0, cex = 0.8, col = "black", "True: green = CpG island")
for (i in 1:nSim) {
  if (Island[i] == TRUE)
    rect(i, -1, i + 1, 0, col = "green", border = NA)
  else rect(i, -1, i + 1, 0, col = "red", border = NA)
}

########Plot the most probable sequence of underlying true states (i.e.whether we were in a CpG island) calculated by the viterbi call#######################
text(0, -3.2, adj = 0, cex = 0.8, col = "black", "Most probable path")
for (i in 1:nSim) {
  if ((x$vit[i] == "A+")||(x$vit[i] == "C+")||(x$vit[i] == "G+")||(x$vit[i] == "T+"))
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

#################Plot the marginal posterior probability that we are in a CpG island at each step#########################
PosteriorPCpG<-colSums(x$posterior[1:4,])/colSums(x$posterior[1:8,])
points(((PosteriorPCpG-4.6)/colSums(x$posterior)) , type = "l")

###############Show when the posterior probs infer the wrong underlying state:############
text(0, -7.2, adj = 0, cex = 0.8, col = "black", "Difference by posterior-probability")
InferredCpG<-colSums(x$posterior[1:4,])>colSums(x$posterior[5:8, ])
#differing = !(x$sim$states == x$vit)
for (i in 1:nSim) {
  if (InferredCpG[i]) {
      if ((x$sim$states[i] == "A+")||(x$sim$states[i] == "C+")||(x$sim$states[i] == "G+")||(x$sim$states[i] == "T+"))
      rect(i, -7, i + 1, -6, col = rgb(0.9, 0.9, 0.9),
           border = NA)
    else rect(i, -7, i + 1, -6, col = rgb(0.3, 0.3, 0.3),
              border = NA)
  }
  else {
    if ((x$sim$states[i] == "A-")||(x$sim$states[i] == "C-")||(x$sim$states[i] == "G-")||(x$sim$states[i] == "T-"))
      rect(i, -7, i + 1, -6, col = rgb(0.9, 0.9, 0.9),
           border = NA)
    else rect(i, -7, i + 1, -6, col = rgb(0.3, 0.3, 0.3),
              border = NA)
  }
}





# Now let's pretend we didn't know the underlying transition matrix
# we use the Baum-Welch algorithm
BW<-baumWelch(hmm, observation=sim$observation, maxIterations=100, delta=1E-9, pseudoCount=0)

# let's compare actual and estimated transition probabilities/ Start by vectorizing both matrices
Act<-c(hmm$transProbs)
Est<-c(BW$hmm$transProbs)
plot(Act,Est)
abline(0,1,lty=2)

# Here's an example application of bw to a simpler HMM that works
# Initial HMM
hmm2 = initHMM(c("A","B"),c("L","R"),
              transProbs=matrix(c(.9,.1,.1,.9),2),
              emissionProbs=matrix(c(.5,.51,.5,.49),2))
print(hmm2)
# Sequence of observation
a = sample(c(rep("L",1000),rep("R",3000)))
b = sample(c(rep("L",3000),rep("R",1000)))
observation = c(a,b)
# Baum-Welch
bw2 = baumWelch(hmm2,observation,10)
print(bw2$hmm)
# let's compare actual and estimated transition probabilities/ Start by vectorizing both matrices
Act<-c(hmm2$transProbs)
Est<-c(bw2$hmm$transProbs)
plot(Act,Est)
abline(0,1,lty=2)
