#'#' @description
#' This function is used to colour the branches in the plotting function when plotting the reconstructed tree.
#'
#'
#' @import ape
#' @import FossilSim

get_colours<- function(data){

  recon <- FossilSim::reconstructed.tree.fossils.objects(data$fossil, data$tree)
  tps <- unname(recon$tree$tip.label)

  matches <- grepl("_1$", tps )

  # Step 2: Extract and clean elements that match
  reconTreeTips <- gsub("_1$", "", tps[matches])

  transition_history$tree$tip.label

  # which tips are not in the reconstructred tree
  setdiff(transition_history$tree$tip.label, reconTreeTips)

  b.colours <- numeric(length(t$edge.length))


  # Go through all the tips to determine are the in the reconstructed tree
  for ( Nt in 1:ape::Ntip(t)){


    if (any(t$tip.label[Nt] == reconTreeTips)) {
      branch.path <- ape::nodepath(t, from = Nt, to = ape::Ntip(t) + 1)

      for (b in 1:length(branch.path)){
        b.colours[which(t$edge[,2] == branch.path[b] )] = "black"
      }
    }
  }


  # now go through what is left to determine if there are fossils along the branch
   # this needs to be updated to change the colour along the branch
  remaining.branches <-which(b.colours == 0)

  for (rb in 1:length(remaining.branches)){

    if(any(transition_history$fossil$ape.branch == remaining.branches[rb])){
      Nt <- t$edge[remaining.branches[rb],][2]
      branch.path <- ape::nodepath(t, from = Nt, to = ape::Ntip(t) + 1)

      for (b in 1:length(branch.path)){
        b.colours[which(t$edge[,2] == branch.path[b] )] = "black"
      }

    }

  }

  b.colours[b.colours == "0"] <- "grey"
  return(b.colours)
}
