plot.morpho = function(x, tree, show.tree = TRUE,
                     # age info/options
                     max.age = NULL, show.axis = TRUE,
                     # taxonomy
                     rho = 1,
                     # tree appearance
                     root.edge = TRUE, hide.edge = FALSE, edge.width = 1, show.tip.label = FALSE, align.tip.label = FALSE,
                     # taxa appearance
                     extant.col = 1, taxa.palette = "viridis",
                     col.axis = "gray35", cex = 1.2, pch = evol$bp_change, ...) {
  
  # hard coded options for tree appearance
  evol = x
  edge.color = "black"
  edge.lty = 1
  font = 3 # italic
  tip.color = "black"
  label.offset = 0.02
  align.tip.label.lty = 3
  underscore = FALSE
  plot = TRUE
  node.depth = 1
  no.margin = FALSE
  adj = NULL
  srt = 0
  strata = 1
  
  if(!show.tree) align.tip.label = TRUE
  
  if(class(evol)[1]!="evol")
    stop("x must be an object of class \"evol\"")
  
  if(!"phylo" %in% class(tree))
    stop("tree must be an object of class \"phylo\"")
  
  if(!all( as.vector(na.omit(evol$edge)) %in% tree$edge))
    stop("Mismatch between evol and tree objects")
  
  if(!(rho >= 0 && rho <= 1))
    stop("rho must be a probability between 0 and 1")
  
  # tolerance for extant tips and interval/ fossil age comparisons
  tol = min((min(tree$edge.length)/100), 1e-8)
  
  # If there are no extant samples, simulate extant samples
  #if(!any( abs(evol$hmax) < tol )){
  #do we need this?
  # evol = sim.extant.samples(evol, tree = tree, rho = rho)
  #}
  
  if(is.null(tree$root.edge))
    root.edge = FALSE
  
  if(is.null(max.age))
    ba = tree.max(tree, root.edge = root.edge)
  else ba = max.age
  
  offset = 0 # distance from youngest tip to present
  if(!is.null(tree$origin.time)) offset = min(n.ages(tree))
  
  # check the tree
  Ntip <- length(tree$tip.label)
  if (Ntip < 2) {
    warning("found less than 2 tips in the tree")
    return(NULL)
  }
  
  if (any(tabulate(tree$edge[, 1]) == 1))
    stop("there are single (non-splitting) nodes in your tree; you may need to use collapse.singles()")
  
  if(!ape::is.rooted(tree))
    stop("tree must be rooted")
  
  if(is.null(tree$edge.length))
    stop("tree must have edge lengths")
  
  Nnode <- tree$Nnode
  if (any(tree$edge < 1) || any(tree$edge > Ntip + Nnode))
    stop("tree badly conformed; cannot plot. Check the edge matrix.")
  ROOT <- Ntip + 1
  
  # collect data for plotting the tree
  type = "phylogram"
  direction = "rightwards"
  horizontal = TRUE # = "rightwards"
  
  xe = tree$edge # used in the last part of the fxn
  yy = numeric(Ntip + Nnode)
  TIPS = tree$edge[tree$edge[, 2] <= Ntip, 2]
  yy[TIPS] = 1:Ntip
  
  yy = ape::node.height(tree)
  xx = ape::node.depth.edgelength(tree)
  
  if(root.edge)
    xx = xx + tree$root.edge
  
  # x.lim is defined by max(ba, tree age)
  if(ba - offset > max(xx))
    xx = xx + (ba - offset - max(xx))
  
  x.lim = c(0, NA)
  pin1 = par("pin")[1]
  strWi = graphics::strwidth(tree$tip.label, "inches", cex = par("cex"))
  xx.tips = xx[1:Ntip] * 1.04
  alp = try(stats::uniroot(function(a) max(a * xx.tips + strWi) - pin1, c(0, 1e+06))$root, silent = TRUE)
  if (is.character(alp)) {
    tmp = max(xx.tips)
    if (show.tip.label)
      tmp = tmp * 1.5
  } else {
    tmp = if (show.tip.label)
      max(xx.tips + strWi/alp)
    else max(xx.tips)
  }
  if (show.tip.label)
    tmp = tmp + label.offset
  x.lim[2] = tmp
  
  y.lim = c(1, Ntip)
  
  use.species.ages = FALSE
  
  s1 = (ba - offset) / strata # horizon length (= max age of youngest horizon)
  horizons.max = seq(s1 + offset, ba, length = strata)
  horizons.min = horizons.max - s1
  s1 = horizons.max - horizons.min
  
    old.par = par("xpd", "mar")
    par(xpd=NA)
    par(mar=c(5,4,2.5,2) + 0.1) ## default is c(5,4,4,2) + 0.1
    # open a new plot window
    graphics::plot.default(0, type = "n", xlim = x.lim, ylim = y.lim, xlab = "", pch=pch, ylab = "", axes = FALSE, asp = NA, ...)
  
  if (plot) {
    
    # add colored strata
    # rect(xleft, ybottom, xright, ytop)
    if(show.axis){
      
      # y-axis:
      #y.bottom = 0
      #y.top = max(yy)+1
      y.bottom = 0.5
      y.top = max(yy) + 0.5
      
      # x-axis:
      #s1 = ba / strata
      x.left = 0 - (ba - offset - max(xx))
      x.right = x.left + s1[1]
      cc = 1 # color switch
      
      axis.strata = x.left
      
        rect(xleft = x.left, xright = x.right, ybottom = y.bottom, ytop = y.top, col="#ffffff", border=NA)
        x.left = x.right
        x.right = x.left + s1[1]
        cc = cc + 1
        axis.strata = c(axis.strata, x.left)
      
      x.labs = rev(round(c(offset, horizons.max),2))
      labs = x.labs
      
      
      if(show.axis){
        axis(1, col = col.axis, at = axis.strata, labels = labs, lwd = 2, line = 0.5, col.axis = col.axis, cex.axis = .9)
          mtext(1, col = 'grey35', text="Time before present", line = 3, cex = 1.2)
      }
    }
    
    # plot the tree
    if(show.tree)
      ape::phylogram.plot(tree$edge, Ntip, Nnode, xx, yy, horizontal, edge.color, edge.width, edge.lty)
    
    # format the root edge
    if (root.edge && show.tree && !hide.edge) {
      rootcol = if(length(edge.color) == 1)
        edge.color
      else "black"
      rootw = if(length(edge.width) == 1)
        edge.width
      else 1
      rootlty = if(length(edge.lty) == 1)
        edge.lty
      else 1
      
      # plot the root edge
      segments(xx[ROOT] - tree$root.edge, yy[ROOT], xx[ROOT], yy[ROOT], col = rootcol, lwd = rootw, lty = rootlty)
    }
    
    # format tip labels
    if (show.tip.label) {
      
      # x and y co-ordinates
      adj = 0
      MAXSTRING <- max(graphics::strwidth(tree$tip.label, cex = cex))
      loy <- 0
      lox <- label.offset + MAXSTRING * 1.05 * adj
      
      if (is.expression(tree$tip.label))
        underscore <- TRUE
      if (!underscore)
        tree$tip.label <- gsub("_", " ", tree$tip.label)
      
      if (align.tip.label) {
        xx.tmp <- max(xx[1:Ntip])
        yy.tmp <- yy[1:Ntip]
        segments(xx[1:Ntip], yy[1:Ntip], xx.tmp, yy.tmp, lty = align.tip.label.lty)
      }
      else {
        xx.tmp <- xx[1:Ntip]
        yy.tmp <- yy[1:Ntip]
      }
      text(xx.tmp + lox, yy.tmp + loy, tree$tip.label,adj = adj, font = font, srt = srt, cex = cex, col = tip.color)
    }
    evol$r = max(xx) - evol$hmax + offset 
    evol = na.omit(evol)
    points(evol$r,yy[evol$edge], cex = cex, pch = pch)
  }
  
  if(!is.na(old.par[1])) par(old.par)
  L <- list(type = type, use.edge.length = TRUE,
            node.pos = NULL, node.depth = node.depth, show.tip.label = show.tip.label,
            show.node.label = FALSE, font = font, cex = cex,
            adj = adj, srt = srt, no.margin = no.margin, label.offset = label.offset,
            x.lim = x.lim, y.lim = y.lim, direction = direction,
            tip.color = tip.color, Ntip = Ntip, Nnode = Nnode, root.time = tree$root.time,
            align.tip.label = align.tip.label)
  assign("last_plot.phylo", c(L, list(edge = xe, xx = xx, yy = yy)),
         envir = ape::.PlotPhyloEnv)
  invisible(L)
  
  
}
