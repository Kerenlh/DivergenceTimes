# "Method 1:"
# The function estimate_M1 uses a maximum likelihood approach to evaluate the distance between two 
# sequences. This method uses the modulo of the number of transitions between the two sequences
# and the number of transitions per site evaluated from an uncensored phylogenetic tree 
# (for example, Phylotree for mtDNA). The distance calculated by the function is calibrated 
# so that the total edge length of the relevant phylogenetic tree is set to 1.
# Input variables:
# xs = Observed number of transitions per site (evaluated from a full uncensored phylogenetic tree
# like Phylotree in the case of mtDNA).
# zs = modulo of the number of transitions between the two sequences, i.e., 0 for sites identical
# in both sequences and 1 for sites that differ by transition (A<->G, C<->T).
# lambda_min_init_val = minimal value for lambdas in the optimization part, default value is 0.0001.
#
# Important note: sites with observed transversions (in the phylogenetic tree or between the 
# two sequences) should be removed from both xs and zs!


# Functions:

estimate_M1 = function(xs,zs,lambda_min_init_val= 0.0001){
  all_x_z_pairs = apply(cbind(xs,zs),1,paste0,collapse = ".")
  x_z_pairs = unique(all_x_z_pairs)
  x_z_pairs_nums = 0*c(1:length(x_z_pairs))
  for (i in 1:length(x_z_pairs)){
    x_z_pairs_nums[i] = length(which(all_x_z_pairs==x_z_pairs[i]))
  }
  p = optimize(df_loglik_for_p, xs = xs, zs = zs, lambda_min_init_val = lambda_min_init_val,
               x_z_pairs = x_z_pairs, x_z_pairs_nums = x_z_pairs_nums, interval = c(0.0001,2),maximum = TRUE)
  return(p$maximum)
}
df_loglik_for_p = function(p,xs,zs,lambda_min_init_val,x_z_pairs,x_z_pairs_nums){
  loglik = solve_likelihood_for_lambdas(xs,zs,p,lambda_min_init_val,x_z_pairs,x_z_pairs_nums)
  return(loglik)
}

solve_likelihood_for_lambdas = function(xs, zs, p,lambda_min_init_val,x_z_pairs,x_z_pairs_nums){
  lambdas = 0*xs
  loglik = 0
  for (i in c(1:length(x_z_pairs))){
    tmp = strsplit(x_z_pairs[i],split = "\\.")[[1]]
    curr_x = as.numeric(tmp[1]); curr_z = as.numeric(tmp[2])
    lambda_i = solve_likelihood_for_lambda(curr_x, curr_z, p, lambda_min_init_val)
    loglik = loglik + x_z_pairs_nums[i]*(log(1 + (-1)^curr_z*exp(-2*p*lambda_i))) - x_z_pairs_nums[i]*lambda_i
    if (curr_x!=0){
      loglik = loglik + x_z_pairs_nums[i]*curr_x*log(lambda_i)
    }
  }
  return (loglik)
}

df_likelihood_lmd = function(lmd,x,p,z){
  logLik_lmd_i = -lmd
  if (x!=0){
    logLik_lmd_i = logLik_lmd_i + x*log(lmd)
  }
  if (z==0){
    logLik_lmd_i = logLik_lmd_i+log(1 + exp(-2*lmd*p))
  }else{
    logLik_lmd_i = logLik_lmd_i+log(1 - exp(-2*lmd*p))
  }
  if (lmd<0.1){
    print(paste(lmd, logLik_lmd_i, x, z))
  }
  return(logLik_lmd_i)
}

solve_likelihood_for_lambda = function(x, z, p,lambda_min_init_val){
  if (x==0 & z==0){
    lambda = 0
  }else if (x==0 & z==1){
    lambda = log(2*p+1)/(2*p)
  }else {
    lambda = optim(par = max(x,lambda_min_init_val),  x = x, z = z, p = p, fn = df_likelihood_lmd,
                   method = "BFGS",gr = df_likelihood_derivative_lmd, control=list(fnscale=-1))
    lambda = lambda$par
  }
  return (lambda)
}

df_likelihood_derivative_lmd = function(lmd,x,p,z){
  derivative = -1 
  if (x!=0){
    derivative = derivative + x/lmd
  }
  if (z==0){
    derivative = derivative - 2*p*exp(-2*lmd*p) /(1 + exp(-2*lmd*p)) 
  }else{
    derivative = derivative + 2*p*exp(-2*lmd*p) /(1 - exp(-2*lmd*p)) 
  }
  return (derivative)
}