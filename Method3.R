# Method 3:
# The function estimate_M3 evaluates the distance between two sequences. 
# The method uses the modulo of the number of transitions between the two sequences
# and the number of transitions per site evaluated from an uncensored phylogenetic tree 
# (for example, Phylotree for mtDNA). The distance calculated by the function is calibrated 
# so that the total edge length of the relevant phylogenetic tree is set to 1.
#
# This method uses a Bayesian approach, assuming a Gamma distribution prior to the transition rates.
#
# Input variables:
# xs = Observed number of transitions per site (evaluated from a full uncensored phylogenetic tree
# like Phylotree in the case of mtDNA).
# zs = modulo of the number of transitions between the two sequences, i.e., 0 for sites identical
# in both sequences and 1 for sites that differ by transition (A<->G, C<->T).
#
# Important note: sites with observed transversions (in the phylogenetic tree or between the 
# two sequences) should be removed from both xs and zs!

library (MASS)
estimate_M3 = function(xs,zs){
  mod = glm.nb (xs~1)
  alpha = mod$theta
  beta = alpha/exp(mod$coefficients[1])
  p = optimize(df_loglik_M3, xs = xs, zs = zs, alpha = alpha, beta = beta, interval = c(0.001,2), maximum = TRUE)
  return(p$maximum)
}

df_loglik_M3 = function(p,xs,zs, alpha, beta){
  tmp = sum(log(1 + ((-1)^zs) * ((beta + 1)/(beta + 1  + 2*p))^(xs+alpha) ))
 return(tmp) 
}


