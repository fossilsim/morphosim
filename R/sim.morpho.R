# https://github.com/KlausVigo/phangorn/blob/master/R/simSeq.R

sim.morpho = function(tree, k = 2, nchar = c(10), variable.coding = FALSE, partitions = 1, alpha = NULL){
  # Mk
  # + V
  # + P
  # + G
  # wish list
  # missing data
  # FBD range
  # note k = max state number

  #TODO check tree etc, e.g., tree must tree, tree must root, or we arbitrarily root the tree, check is fully bifurcating
  #TODO add check partitions and character make sense
  # if (p > 1 && nchar length == 1, nchar[[1]] / p)

  if(is.null(alpha)) rates = 1
  else rates = phangorn::discrete.gamma(alpha, 4) # TODO: NEED TO CHECK WITH MIKE

  for(p in partitions:1){ # for now, goes backwards from p max. Otherwise we have to deal with attr(seq, "levels"); or you can use cbind(x, seq)
    if(max(partitions) == 1) states = as.character(c(0:(k-1)))
    else states = as.character(c(0:p))

    for(i in 1:nchar[p]){
      if(nchar[p] == 0) next

      x = phangorn::simSeq(tree, l = 1, type = "USER", levels = states, rate = sample(rates, 1))

      if(variable.coding){
        while( length(unique(as.character(x))) == 1 ){
          x = phangorn::simSeq(tree, l = 1, type = "USER", levels = states, rate = sample(rates, 1))
        }
      }

      if(i == 1 && p == max(partitions)){
        seq = x
      }
      else seq = cbind(seq, x)
    }
  }
  seq
}


# Step 1. regular RB inference

# Step 2. simulate data in R using the var trees
  # 2a. draw values from Mike's prior

# Step 3. run inference on those replicates RB
  # Check we can estimate alpha

# ?Next steps... do the PP simulations in R?



v = c()
for(i in 1:100000){
  v = c(v, 1/rexp(1))
}
mean(v)
