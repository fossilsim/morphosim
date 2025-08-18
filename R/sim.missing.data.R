#' Remove morphological character data
#'
#' @description
#' This function removes characters from a morphological matrix simulated using morphosim
#' @param data A `morpho` object with sequence data.
#' @param seq Character. Which sequence data to use: "tips", "nodes", or "SA".
#' @param method Character. Method for removing data. Options:
#'   - "random": removes characters randomly across the matrix.
#'   - "partition": removes characters by partition (probabilities per partition).
#'   - "rate": removes characters by rate category (probabilities per rate category).
#'   - "trait": removes characters from specific traits.
#'   - "taxa": removes characters from specific taxa.
#' @param probability Numeric. Probability of missing data (single value or vector depending on method).
#' @param traits When method = "trait", indices of traits to remove.
#' @param taxa When method = "taxa", indices of taxa to remove.
#'
#' @return An object of class morpho.
#'
#' @export

#' @examples
#'
#' # randomly remove data
#' missing.data <- sim.missing.data(data = transition_history,
#'                                   method = "random",
#'                                   seq = "tips",
#'                                   probability = 0.5)
#'
#'
#' # remove data based on the partition
#' missing.data <- sim.missing.data(data = transition_history,
#'                                  method = "partition",
#'                                  seq = "tips",
#'                                  probability = c(0.7, 0, 1, 0.5))
#'
#' # remove data based on the rate it was simulated under
#' missing.data <- sim.missing.data(data = transition_history,
#'                                         method = "rate",
#'                                         seq = "tips",
#'                                         probability = c(0,0,0.2,1))
#'
#' # remove  characters from specific traits
#' missing.data <- sim.missing.data(data = transition_history,
#'                                  method = "trait",
#'                                  seq = "tips",
#'                                  probability = 1,
#'                                  traits = c(1,2,5))
#'
#' # remove  characters from specific taxa
#' missing.data <- sim.missing.data(data = transition_history,
#'                                  method = "taxa",
#'                                  seq = "tips",
#'                                  probability = 1,
#'                                  traits = c(1,3,6))
#'



sim.missing.data <- function(data = NULL, seq = NULL, method = NULL, probability = NULL,
                             traits = NULL, taxa = NULL){


  #### Checks
  if (is.null(method)) {
    stop("You must specify a `method`: 'random', 'partition', 'rate', 'trait', or 'taxa'.")
  }
  if (is.null(data) || !inherits(data, "morpho")) {
    stop("`data` must be a morpho object.")
  }
  if (is.null(seq) || !seq %in% c("tips", "nodes", "SA")) {
    stop("`seq` must be one of 'tips', 'nodes', or 'SA'.")
  }
  if (is.null(probability)) {
    stop("You must provide a `probability` value.")
  }
  if (!is.numeric(probability) || any(probability < 0 | probability > 1)) {
    stop("`probability` must be between 0 and 1.")
  }



  ## Create data frame
  x <- t(as.data.frame(data$sequences[[seq]]))
  trait.num  <- length(x[1,])
  taxa.num <-   length(x[,1])

  ## Method: Random
if (method == "random"){

if (length(probability) > 1)
  stop("For method = 'random', provide a single probability.")
remove <- round((trait.num* taxa.num)* probability, 0)
total_cells <- taxa.num*trait.num
all_combinations <- expand.grid(Row = 1:taxa.num, Column = 1:trait.num)
random_cells<- all_combinations[sample(1:total_cells, remove, replace = FALSE), ]
for ( i in 1:remove){
  x[random_cells$Row[i], random_cells$Column[i]] <- "?"
  }
}


  ## Method: Rate
  if( method == "rate"){
  rates <-  data$model$RateVarTrait
  if (length(probability) != length(unique(rates[1, ]))) {
    stop("Length of `probability` must match the number of rate categories.")
  }

   for ( j in 1:max(rates[1,])){
     traits_per_rate <- which(rates == j)
     remove <- round((length(traits_per_rate)* taxa.num)* probability[j], 0)
     total_cells <- length(traits_per_rate)*taxa.num

     all_combinations <- expand.grid(Row = 1:taxa.num,
                                     Column =  traits_per_rate)
     random_cells <- all_combinations[sample(1:total_cells,
                                             remove, replace = FALSE), ]
       for ( i in 1:remove){
       x[random_cells$Row[i], random_cells$Column[i]] <- "?"

        }
      }
  }


  ## Method: Partition
  if(method == "partition"){

    if (length(probability) != length(data$model))
      stop("Vector of probabilities does not match the number of partitions")

    start_col <- 1
    for ( j in 1:length(data$model)){
    traits_per_partition <-  as.numeric(sub(".*Part:(\\d+).*", "\\1",
                                            data[["model"]][j]))
    remove <- round((traits_per_partition* taxa.num)* probability[j], 0)
    total_cells <- traits_per_partition*taxa.num

    all_combinations <- expand.grid(Row = 1:taxa.num,
                                    Column =  start_col:(start_col + traits_per_partition -1))
    random_cells <- all_combinations[sample(1:total_cells, remove, replace = FALSE), ]


    for ( i in 1:remove){
      x[random_cells$Row[i], random_cells$Column[i]] <- "?"
     }
    start_col <- start_col + traits_per_partition
   }
  }


  ## Method: Trait
  if ( method == "trait"){
    if (length(probability) > 1) {
      stop("For method = 'trait', provide a single probability.")
    }
    if (is.null(traits)) stop("For method = 'trait', you must specify `traits`.")

    remove <- round((length(traits)* taxa.num)* probability, 0)
    total_cells <- length(traits)*taxa.num
    all_combinations <- expand.grid(Row = 1:taxa.num, Column = traits)
    random_cells <- all_combinations[sample(1:total_cells, remove, replace = FALSE), ]
    for ( i in 1:remove){
      x[random_cells$Row[i], random_cells$Column[i]] <- "?"
    }
  }

  ## Method: Taxa
  if ( method == "taxa"){

    if (length(probability) > 1) {
      stop("For method = 'taxa', provide a single probability.")
    }
    if (is.null(taxa)) stop("For method = 'taxa', you must specify `taxa`.")
    remove <- round((length(taxa)* trait.num)* probability, 0)
    total_cells <- length(taxa)* trait.num
    all_combinations <- expand.grid(Column = 1:trait.num,  Row = taxa)
    random_cells <- all_combinations[sample(1:total_cells, remove, replace = FALSE), ]
    for ( i in 1:remove){
      x[random_cells$Row[i], random_cells$Column[i]] <- "?"
    }
 }


  ## update morpho object
 tax <- rownames(x)
  for ( i in 1:length(tax)){
    data$sequences[[seq]][[tax[i]]] <- x[tax[i],]

  }
  return(data)
}





