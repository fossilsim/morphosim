#' Plots morphological matrix
#'
#' @description
#' This function plots the full morphological matrix assocaited with the character data
#' at the tips of a tree. Requires a morpho object as input.
#' @param data A morpho object
#' @param num.trait default is set to "all" which plots all traits in black font. If you
#' want to focus on a specific trait set it here, e.g. num.trait = 1 and this trait will
#' be highlighted
#' @param col A vector of colors that should be the same length or longer than the number of different character states (k). if not specified, the traits from 0 to 6 can be differentiated
#'
#' @export


#' @examples
#'  # simualte a tree
#' phy <- ape::rtree(10)
#'
#' # simulate character data for the tree
#' simulated_morpho <- sim.morpho(phy, k = 2, trait.num = 10)
#'
#' # plot the character matrix
#' plotMorphoGrid(data = simulated_morpho, num.trait = 1)
#'
plotMorphoGrid <- function(data = NULL, timetree = FALSE, num.trait = "all", col =  c("#fdfdfd", "lightgray", "lightblue", "pink", "yellow", "green", "orange")){

  x <- data

n.taxa <- length(x$tree$tip.label)
n.traits <- length(x$sequences[[1]])
ordered_tree <- reorder(x$tree)

## Are we using a time tree?
if (timetree) {
  tree_string <-  write.tree(x$time.tree)
  tip_labs <- regmatches(tree_string, gregexpr("t\\d+", tree_string))[[1]]
} else {
  tree_string <-  write.tree(x$tree)
  tip_labs <- regmatches(tree_string, gregexpr("t\\d+", tree_string))[[1]]
}


#ordered_tree <- reorder$tree.tip.labels

#tip_labs <- ordered_tree$tree$tip.labels
#ordered_tree <- reorder(x$tree)


par(xaxs = "i", yaxs = "i")
plot(x =c(0,1), y =c(0,1),xaxt = 'n', yaxt = 'n',bty = 'n', pch = '',
     ylab = '', xlab = '', main = "Morphological matrix", cex.main= 1, font.main = 6)



## use ablines instead of grid so there is a bit more control
# need to divide the plot into the number of traits and taxa
xx <- 1/n.traits
yy <- 1/n.taxa

## to get the center of the yy boxes
center_a <- yy * 0.5
center_final <- 1 - center_a
y_labs <- seq(center_a, center_final, yy )

axis(2, at=y_labs,labels=tip_labs,
     col.axis="black", las=2, cex.axis=0.8, lwd ="0")

## to get the center of the xx boxes
center_b <- xx * 0.5
center_final <- 1 - center_b
x_labs <- seq(center_b, center_final, xx )

axis(3, at=x_labs,labels=1:n.traits,
     col.axis="black", las=1, cex.axis=0.8, lwd ="0", pos = 0.95)



for (i in 1:n.traits) {
  for (j in 1:n.taxa) {
    state <- as.numeric(x$sequences[[tip_labs[j]]][i])
    bg_col <- col[state + 1]
# Draw a rectangle for each box
rect(
  xleft = x_labs[i]-xx/2, xright = x_labs[i]+xx/2,
  ybottom = y_labs[j]-yy/2, ytop = y_labs[j]+yy/2,
  col = bg_col, border = "black"
)
  }
}

# add border lines
abline(v=0, col= "grey")
abline(h=0, col= "grey")


## grid up the plot area
for ( i in 1:n.traits){
  abline(v = xx*i, col= "grey")
}


for ( j in 1:n.taxa){
  abline(h = yy*j, col= "grey")
}


## fill in the boxes with the state
for (i in 1:n.traits) {
  for (j in 1:n.taxa) {
    state <- as.numeric(x$sequences[[tip_labs[j]]][i])


    # Add the state text in the box
    if (i == num.trait || num.trait == "all") {
      text(x_labs[i], y_labs[j], state, cex = 1, col = "black")
    } else {
      text(x_labs[i], y_labs[j], state, cex = 1, col = "darkgrey")
    }
  }
}
par(xaxs = "r", yaxs = "r")
}







