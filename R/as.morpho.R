# Define the morpho constructor
morpho <- function(data, tree, model) {
  
  # Create a list to store sequences, tree, and model
  morpho_list <- list(
    sequences = data,
    tree = tree,
    model = model
  )
  
  ## Checks 
  # 1. Does the number of sequences match the number of tips in the tree
  if (length(data) != length(tree$tip.label)) {
    stop("Number of sequences doesn't match the number of tips in the tree.")
  }
  
  # 2. Does each sequence contain the same number of characters?
  if (length(unique(sapply(data, length))) != 1) {
    stop("Each sequence should contain the same number of characters.")
  }
  
  # Assign class "morpho" to the object
  #attr(morpho, "class") <- c("morphology", class(morpho))
  class(morpho_list) <- "morpho"
  return(morpho_list)
}

# Define print method for the morpho object
### this was copied from fossil sim so might want to change it up at some stage!!
print.morpho <- function(x, max.length = 5, ...) {
  # Convert sequences into a data frame
  seq_data <- as.data.frame(x$sequences)
  
  # Print first `max.length` rows (or fewer if max.length exceeds the number of rows)
  print(head(seq_data, max.length))
  
  # Print a summary of the morphological data
  cat("Morphological data for", length(x$sequences), "taxa with", 
      length(x$sequences[[1]]), "traits per taxon and", unique(c(seq_data[,1], seq_data[,2])),
      "as character states\n")

# Define summary method for the morpho object
summary.morpho <- function(object, max.length = 5, ...) {
  # Use print.morpho to show the data
  print(object, max.length = max.length)
  
  # Additional information about the tree and model
  cat("Tree with", length(object$tree$tip.label), "tips.\n")
  cat("Model:", object$model, ".\n")
}

# Define the as.morpho generic function
as.morpho <- function(data, tree, model) {
  UseMethod("as.morpho")
}

# Default method for as.morpho
as.morpho.default <- function(data, tree, model) {
  morpho(data, tree, model)
}

# Define is.morpho function to check if the object is of class 'morpho'
is.morpho <- function(data) {
  inherits(data, "morpho")
}
