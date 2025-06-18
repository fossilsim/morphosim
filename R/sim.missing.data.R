#' Remove morphological character data
#'
#' @description
#' This function removes characters from a morphological matrix simulated using morphosim
#' @param data A morpho object with sequence data
#' @param seq Specify which sequence data you want to use ("tips", "nodes", "SA")
#' @param method Specify the method/variable controlling the removal of data. There are 5 options available. "Random",
#' "Partition", "Rate", "Trait", and "Taxa". Random removes characters at random across the entire matrix according
#' to a single probability. Partition allows you to specify different probabilities of being removed for each partition.
#' The number of probabilities provided here must match the number of partitions the data was simulated under. Rate
#' allows you to specify different probabilities of being removed for each rate category. The number of probabilities
#' provided here must match the number of rates the data was simulated under. Traits allows you to specify a probability
#'for specific traits. Taxa llows you to specify a probability for specific taxa.
#' @param probability The probability of missing data to simulate.
#' @param traits When method = trait, used to specify which trait/s you want to remove data from
#' @param taxa When method = taxa, used to specify which taxon/taxa you want to remove data from
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

  if(is.null(method)) stop("You must specify which method you would like to use")

  if(is.null(data) ||!inherits(data, "morpho")){
    stop("Provide a morphosim object to simulate the missing data")
  }

  x <- t(as.data.frame(data$sequences[[seq]]))
  trait.num  <- length(x[1,])
  taxa.num <-   length(x[,1])

if (method == "random"){

if (length(probability) > 1) stop ("Provide 1 probability to apply across the entire matrix")

remove <- round((trait.num* taxa.num)* probability, 0)
total_cells <- taxa.num*trait.num
all_combinations <- expand.grid(Row = 1:taxa.num, Column = 1:trait.num)
random_cells<- all_combinations[sample(1:total_cells, remove, replace = FALSE), ]


for ( i in 1:remove){
  x[random_cells$Row[i], random_cells$Column[i]] <- "?"
}


}


  if( method == "rate"){

   rates <-  data$model$RateVarTrait

   if (length(probability) != length(unique(rates[1,]))) stop("Vector of probabilities does not match the number of rate categories")

   for ( j in 1:max(rates[1,])){
     traits_per_rate <- which(rates == j)
     remove <- round((length(traits_per_rate)* taxa.num)* probability[j], 0)
     total_cells <- length(traits_per_rate)*taxa.num

     all_combinations <- expand.grid(Row = 1:taxa.num, Column =  traits_per_rate)
     random_cells <- all_combinations[sample(1:total_cells, remove, replace = FALSE), ]


     for ( i in 1:remove){
       x[random_cells$Row[i], random_cells$Column[i]] <- "?"

     }
      }
  }



  if(method == "partition"){

    if (length(probability) != length(data$model)) stop("Vector of probabilities does not match the number of partitions")

    start_col <- 1
    for ( j in 1:length(data$model)){


    traits_per_partition <-  as.numeric(sub(".*Part:(\\d+).*", "\\1", data[["model"]][j]))
    remove <- round((traits_per_partition* taxa.num)* probability[j], 0)
    total_cells <- traits_per_partition*taxa.num

    all_combinations <- expand.grid(Row = 1:taxa.num, Column =  start_col:(start_col + traits_per_partition -1))
    random_cells <- all_combinations[sample(1:total_cells, remove, replace = FALSE), ]


    for ( i in 1:remove){
      x[random_cells$Row[i], random_cells$Column[i]] <- "?"


  }
    start_col <- start_col + traits_per_partition

  }


  }


  if ( method == "trait"){

    if (length(probability) > 1) stop ("Provide 1 probability to apply across all traits")

    remove <- round((length(traits)* taxa.num)* probability, 0)
    total_cells <- length(traits)*taxa.num

    all_combinations <- expand.grid(Row = 1:taxa.num, Column = traits)
    random_cells <- all_combinations[sample(1:total_cells, remove, replace = FALSE), ]
    for ( i in 1:remove){
      x[random_cells$Row[i], random_cells$Column[i]] <- "?"


    }


  }


  if ( method == "taxa"){

    if (length(probability) > 1) stop ("Provide 1 probability to apply across all traits")

    remove <- round((length(taxa)* trait.num)* probability, 0)
    total_cells <- length(taxa)* trait.num

    all_combinations <- expand.grid(Column = 1:trait.num,  Row = taxa)

    random_cells <- all_combinations[sample(1:total_cells, remove, replace = FALSE), ]
    for ( i in 1:remove){
      x[random_cells$Row[i], random_cells$Column[i]] <- "?"


    }


  }




  tax <- rownames(x)
  for ( i in 1:length(tax)){
    data$sequences[[seq]][[tax[i]]] <- x[tax[i],]

  }
  return(data)



}





