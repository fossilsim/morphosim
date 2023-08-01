# https://github.com/KlausVigo/phangorn/blob/master/R/simSeq.R

sim.morpho = function(tree, k = 2, nchar = c(10), variable.coding = FALSE, partitions = 1, alpha = NULL){
  # Mk
  # + V
  # + P
  # + G

  if(is.null(alpha)) rates = 1
  else rates = phangorn::discrete.gamma(alpha, 4) # TODO: NEED TO CHECK WITH BEN

  for(p in partitions:1){ # for now, goes backwards from p max. Otherwise we have to deal with attr(seq, "levels"); or you can use cbind(x, seq)
    if(max(partitions) == 1) states = as.character(c(0:(k-1)))
    else states = as.character(c(0:p))

    for(i in 1:nchar[p]){

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
