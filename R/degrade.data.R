#' Remove morphological character data
#'
#' @description
#' This function removes characters from a morphological matrix simulated using morphosim
#' @param data A morpho object with sequence data
#' @param method This specifies under what criteria data should be removed. Methods include "random" which randomly removes
#' characters across the entire matrix, "partition" which allows different partitions (that were used to simualte the data)
#' to loose varying amounts of characters, "rates" which removes characters depening on the rate the were simulated under,
#' and "character" where the user can specify particular characters to degrade the data.
#' @param percentage This specifies how much data to remove. The number of percentages supplied varies depending on the
#' method chosen. For random the user must supply 1 percentage, and for all other methods must supply 1 percentage
#' per category.
#'
#' @return An object of class morpho.
#'
#' @export

#' @examples
#'# simulated tree
#' phy <- ape::rtree(10)
#'
#'# simulate characters along the branches of the tree
#'transition_history <-  sim.morpho.history(tree = phy,
#'                                          k = 3,
#'                                          trait.num = 30,
#'                                          ancestral = TRUE,
#'                                          ACRV = "gamma",
#'                                          variable = FALSE)
#'
#' # randomly remove data
#' degraded <- degrade.data(data = transition_history, method = "random", percentage = 10)
#'
#'
#' # remove data based on the rate it was simulated under
#' degraded <- degrade.data(data = transition_history, method = "rate", percentage = c(0,0,20,100))


degrade.data <- function(data = NULL, method = NULL, percentage = NULL){

  x <- t(as.data.frame(data$sequences))
  trait.num  <- length(x[1,])
  taxa.num <-   length(x[,1])

if (method == "random"){


remove <- round((trait.num* taxa.num)* (percentage/100), 0)
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
     remove <- round((length(traits_per_rate)* taxa.num)* (percentage[j]/100), 0)
     total_cells <- length(traits_per_rate)*taxa.num

     all_combinations <- expand.grid(Row = 1:taxa.num, Column =  traits_per_rate)
     random_cells <- all_combinations[sample(1:total_cells, remove, replace = FALSE), ]


     for ( i in 1:remove){
       x[random_cells$Row[i], random_cells$Column[i]] <- "?"

     }
      }
     }

  taxa <- rownames(x)
  for ( i in 1:length(taxa)){
    data$sequences[[taxa[i]]] <- x[taxa[i],]

  }

  return(data)

}







