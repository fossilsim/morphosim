#' Mopho object
#'
#' Create a morpho.
#'
#' @param data ?
#' @param tree ?
#' @param tree Tree with branches that represent genetic distance associated with the character data.
#' @param time.tree Tree with branches that represent time associated with the character data.
#'
#' @export
morpho <- function(data = NULL, tree = NULL, time.tree = NULL, model = NULL){

  # check the number of sequences match the number of tips in the tree
  if (length(data) != length(tree$tip.label)) {
    stop("Number of sequences doesn't match the number of tips in the tree.")
  }

  # check each sequence contains the same number of characters
  if (length(unique(sapply(data, length))) != 1) {
    stop("Each sequence should contain the same number of characters.")
  }

  # Create a list to store sequences, tree, and model
  morpho.list <- list(
    sequences = data,
    tree = tree,
    model = model
  )

  # assign class "morpho" to the object
  attr(morpho.list, "class") <- c("morpho", class(morpho.list))
  return(morpho.list)
}

#' @export
#' @aliases morpho
print.morpho <- function(x, max.length = 5, ...){

  # Convert sequences into a data frame
  seq.data <- t(as.data.frame(x$sequences))

  print(seq.data[,1:max.length])

  # Print a summary of the morphological data
  cat("Morphological data for", length(x$sequences), "taxa with",
      length(x$sequences[[1]]), "traits per taxon and", sort(unique(c(seq.data[,1], seq.data[,2]))),
      "as character states\n")
      cat("Showing", max.length, "traits here for now\n")
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
as.morpho <- function(data, tree, model) UseMethod("as.morpho")


#' @export
as.morpho.default <- function(data, ...) {
  morpho(data, ...)
}

#' @export
#' @rdname morpho
is.morpho <- function(data) {
  inherits(data, "morpho")
}
