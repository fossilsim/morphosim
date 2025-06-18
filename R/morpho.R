#' Morpho object
#'
#' Create a morpho object.
#'
#' @param sequences A list containing all of the sequences simulated. This can contain sequences for taxa
#' at the tips or the tree, along the nodes, and if present, for sampled ancestors (SA)
#' @param trees A list containing the trees and branch lengths used for the simulation. EvolTree contains a
#' phylogenetic tree with branch lengths representing evolutionary distance. TimeTree (if present)
#' contains the same tree with branch lengths in unit of time. BrRates can either be a single value, when
#' simulating under a strict colck, or a vector of values representing the rate/branch
#' @param model A list containing all model attributes. Model specifies the components specified to simulate under.
#' RateVar containes the realtive values drawn from the specified distribution. RateVarTrait species the rate used
#' to simualte each trait
#' @param root.states A vector supplying the root state for each character
#' @param fossil Fossil object used to simulate data
#' @param transition_history The constant character transitions along the branches
#'
#' @export
#morpho <- function(data = NULL, tree = NULL, time.tree = NULL, model = NULL,
          #         node.seq = NULL, transition_history = NULL, root.states = NULL, fossil_sequence = NULL,
          #         fossil = NULL, ACRV_rate = NULL, gamma_rates = NULL){

  morpho <- function(data = NULL, trees = NULL, model = NULL,
                     transition_history = NULL, root.states = NULL,
                     fossil = NULL){

  # check the number of sequences match the number of tips in the tree
  #if (length(data) != length(tree$tip.label)) {
  #  stop("Number of sequences doesn't match the number of tips in the tree.")
 # }

  # check each sequence contains the same number of characters
  #if (length(unique(sapply(data, length))) != 1) {
  #  stop("Each sequence should contain the same number of characters.")
  #}

  # Create a list to store sequences, tree, and model
  morpho.list <- list(
    sequences = data,
    trees = trees,
    model = model,
    transition_history = transition_history,
    root.states = root.states,
    fossil = fossil
  )

  # assign class "morpho" to the object
  attr(morpho.list, "class") <- c("morpho", class(morpho.list))
  return(morpho.list)
}

#' @export
#' @aliases morpho
print.morpho <- function(x, max.length = 5, ...){

  # Convert sequences into a data frame
  seq.data <- t(as.data.frame(x$sequences$tips))

  if(ncol(seq.data) > 4){
  print(seq.data[,1:max.length])
  } else {
    print(seq.data)
    }

  # Print a summary of the morphological data
  cat("Morphological data for", length(x$sequences$tips), "taxa with",
      length(x$sequences$tips[[1]]), "traits per taxon and", sort(unique(as.vector(seq.data))),
      "as character states\n")
      cat("Showing maximum of", max.length, "traits here for now\n")
}

#' @export
#' @aliases morpho
summary.morpho <- function(object, max.length = 5, ...){

  # Use print.morpho to show the data
  print(object, max.length = max.length)

  # Additional information about the tree and model
  cat("Tree with", length(object$tree$tip.label), "tips\n")
  cat("Model:", object$model, "\n")
}

#' @export
#' @rdname morpho
as.morpho <- as.morpho <- function(data, trees, model = NULL,
                                   transition_history = NULL, root.states = NULL,
                                   fossil = NULL) UseMethod("as.morpho")


#' @export
as.morpho.default <- function(data, ...) {
  morpho(data, ...)
}

#' @export
#' @rdname morpho
is.morpho <- function(data) {
  inherits(data, "morpho")
}

