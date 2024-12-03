#' This function plots your morpho object to a character matrix
#' @description
#' This function creates a matrix from your morpho object highliting whic tips a which character state
#' @import ggplot2
#' @import reshape2
#' @param x An object of class \code{"morpho"}.
#' @param xlab A string for the x-axis label (default: "Characters").
#' @param ylab A string for the y-axis label (default: "Taxa").
#' @param name A string for the legend title (default: "Character State").
#' @param col A character vector of colors. The length must match the number 
#'   of unique character states in the data.
#' @export 
#' @examples
#' tree <- sim.bd.taxa(8,1,1,0.5)[[1]]
#' x <- sim.morpho.completeprocess(time.tree = tree, br.rates = 0.2, k = 2, num.traits = 5)
#' plot.morpho.grid(x = x)

plotMorphoGrid <- function(x, xlab = "", ylab = "",name = "" ,  col = c("white", "gray", "lightblue", "pink", "yellow", "green", "orange")){
  #tip/taxon labels
  tips<-x[[2]][[2]]
  
  # make empty container
  char_matrix<-matrix(nrow= length(tips),ncol = length(x[[1]][[1]]))
  rownames(char_matrix)<-tips
  
  # add character information to the data frame
  
  for (i in 1:length(x[[1]])){
    
    #traits for a given taxon
    
    traits<-x[[1]][[i]]
    
    #make sure the values are numeric
    
    char_matrix[i,]<-as.numeric(traits)
  }
  
  #character states: Question, do we want to use 0 and 1 always? Having issues coding
  #it to be whatever the character states are
  #charas<-as.numeric(unique(x[[1]][[1]]))
  
  #Plotting
  
  dat <- reshape2::melt(char_matrix)
  
  unique_states <- sort(unique(dat$value))
  color_mapping <- setNames(col[1:length(unique_states)], unique_states)
  
  #Question: What should we name the axes? And legend?
  ggplot2::ggplot(dat, aes(Var2, Var1, fill=factor(value))) +
    ggplot2::geom_tile(color="black") +
    ggplot2::geom_text(aes(label = value), color = "black", size = 4) +
    ggplot2::scale_fill_manual(values = color_mapping, name = name) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = xlab, y = ylab) +  # add it to function parameters :TB
    ggplot2::theme(
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = 10)
    )
}
