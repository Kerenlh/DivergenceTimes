# "Method 2:"
# The function estimate_M2 evaluates the distance between two sequences. 
# The method uses the modulo of the number of transitions between the two sequences
# and the number of transitions per site evaluated from an uncensored phylogenetic tree 
# (for example, Phylotree for mtDNA). The distance calculated by the function is calibrated 
# so that the total edge length of the relevant phylogenetic tree is set to 1.
#
# The transition rates are estimated solely based on xs and then 
# the distance is estimated as if the transition rates are known, so only zs is considered.
# An additional assumption is that the number of transitions between the two sequences denoted as
# ys has a binomial distribution so that y_i|x_i ~ Bin(x_i,p). 
# The index i represents a specific site, and p is the distance between the sequences.
#
# Input variables:
# xs = Observed number of transitions per site (evaluated from a full uncensored phylogenetic tree
# like Phylotree in the case of mtDNA).
# zs = modulo of the number of transitions between the two sequences, i.e., 0 for sites identical
# in both sequences and 1 for sites that differ by transition (A<->G, C<->T).
#
# Important note: sites with observed transversions (in the phylogenetic tree or between the 
# two sequences) should be removed from both xs and zs!

estimate_M2 = function(xs,zs){
  epsilon = 0.0913
  xs[which(xs==0)] = epsilon
  p = optimize(df_loglik_M2, xs = xs, zs = zs, interval = c(0.00001,1))
  return(p$minimum)
}

df_loglik_M2 = function(p,xs,zs){
  tmp = abs(sum((1-2*p)^xs) - length(xs) + 2*sum(zs))
  return(tmp)
}
