#' Color branches for plotting a reconstructed tree
#'
#' @description
#' This function generates colors for branches when plotting a reconstructed tree
#' from a morpho object containing fossil data. Branches that are part of the
#' reconstructed tree or have fossils along them are colored black; all others are grey.
#'
#' @param data A morpho object which contains fossil data and a time-calibrated tree.
#'
#' @return A list of length 2:
#'   \item{b.colours}{Vector of branch colors for plotting.}
#'   \item{rem}{Indices of branches with fossils.}
#'
#' @export
#' @examples
#' # Assuming `morpho_data` is a morpho object with fossil data:
#' # branch_colors <- reconstruct.tree(morpho_data)
reconstruct.tree <- function(data) {

  recon <- FossilSim::reconstructed.tree.fossils.objects(data$fossil, data$trees$TimeTree,
                                                         tip_order = "youngest_first")
  tps <- unname(recon$tree$tip.label)
  matches <- grepl("_1$", tps)

  # Extract elements that match
  reconTreeTips <- gsub("_1$", "", tps[matches])
  b.colours <- numeric(length(data$trees$TimeTree$edge.length))

  # Go through all the tips to determine if they are in the reconstructed tree
  for (Nt in 1:ape::Ntip(data$trees$TimeTree)) {
    if (any(data$trees$TimeTree$tip.label[Nt] == reconTreeTips)) {
      branch.path <- ape::nodepath(
        data$trees$TimeTree,
        from = Nt,
        to = ape::Ntip(data$trees$TimeTree) + 1
      )

      for (b in 1:length(branch.path)) {
        b.colours[which(data$trees$TimeTree$edge[,2] == branch.path[b])] = "black"
      }
    }
  }

  # Now go through remaining branches to check for fossils
  remaining.branches <- which(b.colours == 0)
  rem <- c()
  r.cols <- list(NA, NA)

  for (rb in 1:length(remaining.branches)) {
    if (any(data$fossil$ape.branch == remaining.branches[rb])) {
      Nt <- data$trees$TimeTree$edge[remaining.branches[rb],][2]
      branch.path <- ape::nodepath(
        data$trees$TimeTree,
        from = Nt,
        to = ape::Ntip(data$trees$TimeTree) + 1
      )

      for (b in 1:length(branch.path)) {
        b.colours[which(data$trees$TimeTree$edge[,2] == branch.path[b])] = "black"
        if (branch.path[b] <= ape::Ntip(data$trees$TimeTree)) rem <- rbind(rem, remaining.branches[rb])
      }
    }
  }

  remaining.branches <- which(b.colours == 0)
  for (rb in 1:length(remaining.branches)) {
    if (any(data$fossil$ape.branch == remaining.branches[rb])) {
      rem <- rbind(rem, remaining.branches[rb])
    }
  }

  b.colours[b.colours == "0"] <- "grey"
  r.cols[[1]] <- b.colours
  r.cols[[2]] <- rem

  return(r.cols)
}
