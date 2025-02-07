#' Plot full evolutionary history
#'
#'#' @description
#' This function creates a plot showing continuous evolution of discrete traits
#' @param x A morpho object
#' @param timetree TRUE or FALSE Indicate whether you want to plot a time tree or not. default FALSE, uses distance tree if FALSE
#' @param trait What trait number you want to visualize
#' @param fossil Do you want to plot the fossil along the tree. Default set to FALSE
#' @param br.rates Required if you provide a time tree
#' @param root.edge If TRUE plot the root edge. Default = FALSE
#' @param col A vector of colors that should be the same length or longer than the number of different character states (k). if not specified, the traits from 0 to 6 can be differentiated
#' @param col.timescale a single color for the timescale, "darkgrey" by standard
#' @param ... other arguments to be passed to methods, such as graphical parameters
#' @import ape
#' @import FossilSim
#' @import stats
#' @import graphics
#' @importFrom methods is
#' @import TreeSim
#' @export

#' @examples

#' # Non-time tree
#' phy <- ape::rtree(10)
#'
#' # simulate characters along the branches of the tree
#' transition_history <-  sim.morpho(tree = phy,
#'                                   k = c(2,3,4),
#'                                    trait.num = 20,
#'                                   ancestral = TRUE,
#'                                   partition = c(10,5,5),
#'                                  ACRV = "gamma",
#'                                  variable = TRUE,
#'                                   ncats.gamma = 4)
#'
#' plot(x = transition_history, trait = 1, timetree = FALSE)
#'

plot.morpho <- function(x = NULL, trait = NULL, timetree = FALSE, br.rates = NULL,
                        fossil = FALSE, root.edge = FALSE, col = c("#fdfdfd", "lightgray", "lightblue", "pink", "yellow", "green", "orange"), col.timescale = "darkgrey", ...){

  data <- x
  ## Are we using a time tree?
  if (timetree) {
    if(root.edge){
    plot(data$time.tree, root.edge = T)
    } else { plot(data$time.tree, root.edge = F)}
  } else {
    plot(data$tree)
  }


   # get the plot information
  tree_plot_info <- get("last_plot.phylo", envir = .PlotPhyloEnv)

  edge_start_x <- tree_plot_info$xx[data$tree$edge[, 1]]
  edge_end_x <- tree_plot_info$xx[data$tree$edge[, 2]]

  yy <- tree_plot_info$yy

  edge <- data$tree$edge
  parent <- as.integer(edge[, 1])
  child <- as.integer(edge[, 2])

  # Identify the root
  root <- as.integer(parent[!match(parent, child, 0)][1])


  if(!is.null(trait)){
  df <- data$transition_history[trait][[1]]

  if (nrow(df) > 0) {
    # Loop over the dataframe to add points
    for (i in 1:nrow(df)) {
      branch <- as.numeric(df$edge[i])

      if (timetree) {
        actual_position <- as.numeric(df$hmin[i]) * br.rates
        position <- (actual_position / (data$tree$edge.length[branch] * br.rates))
      } else {
        actual_position <- as.numeric(df$hmin[i])
        position <- (actual_position / data$tree$edge.length[branch])
      }

      # Calculate the point's coordinates along the branch
      point_x <- edge_start_x[branch] + position * (edge_end_x[branch] - edge_start_x[branch])

      if (timetree) {
        point_y <- yy[data$time.tree[["edge"]][branch, 2]]
      } else {
        point_y <- yy[data$tree[["edge"]][branch, 2]]
      }
      paint <- as.numeric(df$state[i]) + 1
      points(point_x, point_y, pch = 22, col = "black", bg = col[paint], cex = 4)
      text(point_x, point_y, labels = as.numeric(df$state[i]))
    }
    points(0, yy[root], pch = 22, col = "black", bg = col[as.numeric(data$root.states[trait])+1], cex = 4)
    text(0, yy[root], label = as.numeric(data$root.states[trait]))

  } else {
    points(0, yy[root], pch = 22, col = "black", bg = col[as.numeric(data$root.states[trait])+1], cex = 4)
    text(0, yy[root], label = as.numeric(data$root.states[trait]))
   # message("No transitions in this state across taxa")
  }
  }
  ## add fossils
  if(fossil){

    tree.age <- max(ape::node.depth.edgelength(data$time.tree))
    #root.age <- tree.age - data$time.tree$root.edge
    for (fsl in 1:length(data$fossil$sp)){
      branch <- data$fossil$ape.branch[fsl]

      #how far along the tree is the fossil
      time <- as.numeric(data$fossil$hmax[fsl])
      if (root.edge){
      actual_position <- tree.age - time + data$time.tree$root.edge
      } else{
        actual_position <- tree.age - time
      }
      position <- actual_position / data$time.tree$edge.length[branch]
      # Calculate the point's coordinates along the branch
      point_x <-  position * (edge_end_x[branch] - edge_start_x[branch])
      point_y <- yy[data$time.tree[["edge"]][branch, 2]]

      points(point_x, point_y, pch = 18, col = "black", cex = 1)

    }
  }






  # Add a timescale below the plot

if (timetree) {

  if(root.edge){
    axis_labels <- c((tree.age + data$time.tree$root.edge),0)
    axis(1, at = c(0, (tree.age + data$time.tree$root.edge)), labels = round(axis_labels, 2),
                     line = 1,col = col.timescale, lwd = 3, cex.axis = 1.3, col.axis = col.timescale)
  } else {
    axis_labels <- c(tree.age,0)
    axis(1, at = c(0, tree.age), labels = round(axis_labels, 2),
         line = 1,col = col.timescale, lwd = 3, cex.axis = 1.3, col.axis = col.timescale)
    }


  mtext("Time before present", side = 1, line = 2.5, cex = 1.3, col = col.timescale)
}
}
