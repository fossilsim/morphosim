#'@param function to create a grid showing taxon character traits from simulated data
#'@param x morpho object
#'@param col1 Color for character state: "1"
#'@param col2 Color for character state: "2"

plot.morpho.grid <- function(x, col1 = "white", col2 ="gray"){

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

    char_matrix[i,]<-as.numeric(traits)}

  #character states: Question, do we want to use 0 and 1 always? Having issues coding
  #it to be whatever the character states are
  #charas<-as.numeric(unique(x[[1]][[1]]))

  #Plotting

  dat <- reshape2::melt(char_matrix)

  #Question: What should we name the axes? And legend?
  ggplot2::ggplot(dat, aes(Var2, Var1, fill=factor(value))) +
    ggplot2::geom_tile(color="black") +
    ggplot2::scale_fill_manual(values = c("1"=col1, "0"=col2), name="Character") +
    ggplot2::theme_minimal() +
    ggplot2::labs(x="", y="Species") +
    ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 1))

}
