#' Plot full evolutionary history
#'
#'#' @description
#' This function creates a plot showing continuous evolution of discrete traits
#' @param data A morpho object
#' @param timetree TRUE or FALSE Indicate whether you want to plot a time tree or not. default FALSE, uses distance tree if FALSE
#' @param trait What trait number you want to visualize
#' @param fossil Plot the fossil along the tree. Default set to FALSE
#' @param br.rates Required if you provide a time tree
#' @param root.edge If TRUE plot the root edge. Default = FALSE
#' @param reconstructed Plot the reconstructed tree. Default = FALSE
#' @param col A vector of colors that should be the same length or longer than the number of different character states (k). if not specified, the traits from 0 to 6 can be differentiated
#' @param col.timescale a single color for the timescale, "darkgrey" by standard
#' @param ... other arguments to be passed to methods, such as graphical parameters
#' @import ape
#' @import FossilSim
#' @import stats
#' @import graphics
#' @import TreeSim
#' @export

#' @examples
#' # simulate a phylogenetic tree
#' phy <- ape::rtree(10)
#'
#' # simulate characters along the branches of the tree
#' transition_history <-  sim.morpho(tree = phy,
#'                                   k = c(2,3,4),
#'                                    trait.num = 20,
#'                                   ancestral = TRUE,
#'                                   partition = c(10,5,5),
#'                                   ACRV = "gamma",
#'                                   variable = TRUE,
#'                                   ACRV.ncats = 4,
#'                                   define.Q = NULL)
#'
#' plot(transition_history_1, trait= 4, timetree = T, fossil = T,
#' root.edge = T, reconstructed = T)
#'

plot <- function(data = NULL, trait = NULL, timetree = FALSE,
                        fossil = FALSE, reconstructed = FALSE, root.edge = FALSE, edges = 1,
                        label.offset = 0.05, e.cex = 0.5, f.cex = 1,col = c("#fdfdfd", "lightgray", "lightblue", "pink", "yellow", "green", "orange"), col.timescale = "darkgrey", ...){


  ## Are we using a time tree?
  if(reconstructed == FALSE){
  if (timetree) {
    if(root.edge){
    plot(data$trees$TimeTree, edge.width = edges, label.offset = label.offset, root.edge = T)
    } else { plot(data$trees$TimeTree, label.offset = label.offset, edge.width = edges, root.edge = F)}
  } else {
    plot(data$trees$EvolTree, label.offset = label.offset, edge.width = edges)
  }
  } else {
    b.cols <- reconstruct.tree(data)
    if (timetree) {
      if(root.edge){
      plot(data$trees$TimeTree, root.edge = T, label.offset = label.offset,edge.width = edges, edge.color = b.cols[[1]])
    } else { plot(data$trees$TimeTree, root.edge = F,label.offset = label.offset, edge.width = edges, edge.color = b.cols[[1]])}
  } else {
    plot(data$trees$EvolTree, edge.width = edges,label.offset = label.offset, edge.color = b.cols[[1]])
  }
  }


   # get the plot information
  tree_plot_info <- get("last_plot.phylo", envir = .PlotPhyloEnv)

  edge_start_x <- tree_plot_info$xx[data$trees$EvolTree$edge[, 1]]
  edge_end_x <- tree_plot_info$xx[data$trees$EvolTree$edge[, 2]]

  yy <- tree_plot_info$yy

  edge <- data$trees$EvolTree$edge
  parent <- as.integer(edge[, 1])
  child <- as.integer(edge[, 2])

  # Identify the root
  root <- as.integer(parent[!match(parent, child, 0)][1])

  if(!is.null(trait)){
  df <- data$transition_history[trait][[1]]


  if(root.edge){
    points(data$trees$TimeTree$root.edge, yy[root], pch = 22, col = "black", bg = col[as.numeric(data$root.states[trait])+1], cex = 4)
    text(data$trees$TimeTree$root.edge, yy[root], label = as.numeric(data$root.states[trait]))
  }else{
  points(0, yy[root], pch = 22, col = "black", bg = col[as.numeric(data$root.states[trait])+1], cex = 4)
  text(0, yy[root], label = as.numeric(data$root.states[trait]))}


  if (nrow(df) > 0) {
    # Loop over the dataframe to add points
    for (i in 1:nrow(df)) {
      branch <- as.numeric(df$edge[i])

      if (timetree) {
      position <- as.numeric(df$hmin[i]) / data$trees$BrRates[branch]
        } else {
      position <- as.numeric(df$hmin[i])
      }
      point_x <- edge_start_x[branch] + position


      if (timetree) {
        point_y <- yy[data$trees$TimeTree[["edge"]][branch, 2]]
      } else {
        point_y <- yy[data$trees$EvolTree[["edge"]][branch, 2]]
      }
      paint <- as.numeric(df$state[i]) + 1
      points(point_x, point_y, pch = 22, col = "black", bg = col[paint], cex = 4)
      text(point_x, point_y, labels = as.numeric(df$state[i]))
    }
 }
}

  tree.age <- max(ape::node.depth.edgelength(data$trees$TimeTree))
  ################
  ## colour branches
  ###############
  if(reconstructed){

    if (length(b.cols) == 2){

    for (p in 1:length(b.cols[[2]]) ){
      q <-  which(data$fossil$ape.branch == b.cols[[2]][p])
      if (length(q) > 1)  {
        fossil_pos <- data$fossil$hmin[which(data$fossil$hmin == min(data$fossil$hmin[q]))]
      } else {
        fossil_pos <- data$fossil$hmin[q]
      }

      ## add lines
     tree.age <- max(ape::node.depth.edgelength(data$trees$TimeTree))
     branch <- data$fossil$ape.branch[q]

      #how far along the tree is the fossil
      time <- fossil_pos
      if (root.edge){
        actual_position <- tree.age - time + data$trees$TimeTree$root.edge
      } else{
        actual_position <- tree.age - time
      }
      position <- actual_position / data$tree$TimeTree$edge.length[branch]
      # Calculate the point's coordinates along the branch
      point_x <-  position * (edge_end_x[branch] - edge_start_x[branch])
      y_vals <- yy[data$trees$TimeTree[["edge"]][b.cols[[2]][p], 2]]
      segments(point_x, y_vals, edge_end_x[b.cols[[2]][p]], y_vals, col = "grey")
    }
  }
}

  ################
  ## add fossils
  ###############
  if(fossil){

    tree.age <- max(ape::node.depth.edgelength(data$trees$TimeTree))
    #root.age <- tree.age - data$time.tree$root.edge
    for (fsl in 1:length(data$fossil$sp)){
      branch <- data$fossil$ape.branch[fsl]
      #how far along the tree is the fossil
      time <- as.numeric(data$fossil$hmax[fsl])
      if (root.edge){
        actual_position <- tree.age - time + data$trees$TimeTree$root.edge
      } else{
        actual_position <- tree.age - time
      }
      position <- actual_position / data$trees$TimeTree$edge.length[branch]
      # Calculate the point's coordinates along the branch
      point_x <-  position * (edge_end_x[branch] - edge_start_x[branch])
      point_y <- yy[data$trees$TimeTree[["edge"]][branch, 2]]

      if (data$fossil$hmax[fsl] == 0){

      points(point_x, point_y, pch = 16, col = "forestgreen", cex = e.cex)
      } else {
        points(point_x, point_y, pch = 18, col = "black", cex = f.cex)
      }
    }
  }

  # Add a timescale below the plot

if (timetree) {

  if(root.edge) {
    axis_labels <- c((tree.age + data$trees$TimeTree$root.edge),0)
    axis(1, at = c(0, (tree.age + data$trees$TimeTree$root.edge)), labels = round(axis_labels, 2),
                     line = 1,col = col.timescale, lwd = 3, cex.axis = 1.3, col.axis = col.timescale)
  } else {
    axis_labels <- c(tree.age,0)
    axis(1, at = c(0, tree.age), labels = round(axis_labels, 2),
         line = 1,col = col.timescale, lwd = 3, cex.axis = 1.3, col.axis = col.timescale)
    }
  mtext("Time before present", side = 1, line = 2.5, cex = 1.3, col = col.timescale)
  }
}
