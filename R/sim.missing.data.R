#' Remove morphological character data
#'
#' @description
#' This function removes characters from a morphological matrix simulated using morphosim
#' @param data A morpho object with sequence data
#' @param method Under what criteria data should be removed. Methods include "random" which randomly removes
#' characters across the entire matrix, "partition" which allows different partitions (that were used to simulate the data)
#' to loose varying amounts of characters, "rates" which removes characters based on the rate the were simulated under,
#' and "trait" where the user can specify particular traits to remove data.
#' @param probability The probability of missing data to simulate.
#' @param traits When method = trait, used to specify which trait you want to remove data from
#'
#' @return An object of class morpho.
#'
#' @export

#' @examples
#'# simulated tree
#' phy <- ape::rtree(10)
#'
#'# simulate characters along the branches of the tree
#'transition_history <-  sim.morpho(tree = phy,
#'                                          k = 3,
#'                                          trait.num = 30,
#'                                          ancestral = TRUE,
#'                                          ACRV = "gamma",
#'                                          variable = FALSE)
#'
#' # randomly remove data
#' missing.data <- sim.missing.data(data = transition_history,
#'                                   method = "random",
#'                                   probability = 1)
#'
#' # remove data based on the rate it was simulated under
#' missing.data <- sim.missing.data(data = transition_history,
#'                                         method = "rate",
#'                                         probability = c(0,0,0.2,1))
#'
#' # remove specific characters from specific traits
#' missing.data <- sim.missing.data(data = transition_history,
#'                                  method = "trait",
#'                                  probability = 1,
#'                                  traits = c(1,2,5))
#'
#' # remove data based on the partition
#' missing.data <- sim.missing.data(data = transition_history,
#'                                  method = "partition",
#'                                  probability = c(1, 0, 1))


sim.missing.data <- function(data = NULL, method = NULL, probability = NULL, traits = NULL){

  if(is.null(method)) stop("Please specify which method you would like to use")

  if(is.null(data) ||!inherits(data, "morpho")){
    stop("Please provide a morphosim object to simulate the missing data")
  }

  x <- t(as.data.frame(data$sequences))
  trait.num  <- length(x[1,])
  taxa.num <-   length(x[,1])

if (method == "random"){



remove <- round((trait.num* taxa.num)* probability, 0)
total_cells <- taxa.num*trait.num
all_combinations <- expand.grid(Row = 1:taxa.num, Column = 1:trait.num)
random_cells<- all_combinations[sample(1:total_cells, remove, replace = FALSE), ]


for ( i in 1:remove){
  x[random_cells$Row[i], random_cells$Column[i]] <- "?"
}


}


  if( method == "rate"){
   rates <-  data$ACRV_rate

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

    remove <- round((length(traits)* taxa.num)* probability, 0)
    total_cells <- length(traits)*taxa.num

    all_combinations <- expand.grid(Row = 1:taxa.num, Column = traits)
    random_cells <- all_combinations[sample(1:total_cells, remove, replace = FALSE), ]
    for ( i in 1:remove){
      x[random_cells$Row[i], random_cells$Column[i]] <- "?"


    }


  }

  taxa <- rownames(x)
  for ( i in 1:length(taxa)){
    data$sequences[[taxa[i]]] <- x[taxa[i],]

  }
  return(data)



}





