#' Write reconstructed tree to file
#'
#' @description
#' Write the reconstructed tree to Newick string
#'
#' @param data Morpho object
#' @param file File name
#'
#' @export
#'
write.recon.tree <- function (data, file) {

r_tree <- FossilSim::reconstructed.tree.fossils.objects(fossils  = data$fossil,
                                                        tree = data$trees$TimeTree)
ape::write.tree(r_tree$tree, file = file)
}

#' Write reconstructed Matrix to file
#'
#' @description
#' Write the reconstructed matrix to a nexus file
#'
#' @param data Morpho object
#' @param file File name
#' @param keep_matrix Logical. If TRUE, returns a matrix showing the naming transformations
#' between `morphosim` and `fossilsim` of Sampled ancestors.
#'
#' @export
#'
write.recon.matrix <- function (data, file = NULL, keep_matrix = F) {
  r_tree <- FossilSim::reconstructed.tree.fossils.objects(fossils  = data$fossil,
                                                          tree = data$trees$TimeTree)
  SA_tips <- c()
  tps <- unname(r_tree$tree$tip.label)
  matches <- grepl("_1$", tps )
  # Extract elements that match
  reconTreeTips <- gsub("_1$", "", tps[matches])

  ## add these tip labels to the file + plus all sampled ancestor
  seq_tips <- which(names(data$sequences$tips) %in% reconTreeTips)

  matches <- grepl("_2$", tps )
  reconSA <-  gsub("_2$", "", tps[matches])

  transformation <- matrix(ncol = 2, nrow = length(reconSA))
  colnames(transformation) <- c("Morphosim", "Fossilsim")
  if (length(reconSA) > 0){

  for (l in 1:length(reconSA)){
    t_label <- which(data$trees$TimeTree$tip.label == reconSA[l])
    b_num <- which(data$trees$TimeTree$edge[,2] == t_label)
    spec_min <- min(data$fossil$hmin[data$fossil$ape.branch == b_num])
    spec_num <- data$fossil$specimen[ data$fossil$hmin == spec_min ]
    SA_tips <- rbind(SA_tips, c(paste0(spec_num, "_", b_num)))

    transformation[l,"Morphosim"] <- SA_tips[l]
    transformation[l,"Fossilsim"] <- reconSA[l]

  }

  total_tips <- c(data$sequences$tips[c(seq_tips)], data$sequences$SA[c(SA_tips)])
  } else {
    total_tips <- data$sequences$tips[c(seq_tips)]
  }



  ## need to change the sequence names to match the reconstructed tree
  for ( l in 1:length(seq_tips)){
    current <- names(data$sequences$tips[c(seq_tips)[l]])
    names(total_tips)[names(total_tips) == current ] <- paste0(current, "_1")
  }

  if (length(reconSA) > 0){
  for( l in 1:length(SA_tips)){
    rematch <- unname(transformation[l,"Fossilsim"])
    names(total_tips)[names(total_tips) == SA_tips[l]] <- paste0(rematch, "_2")
  }
  }

  if(!is.null(file)) {
  ape::write.nexus.data(total_tips , file = file)
  }

  if(keep_matrix == T){
    return(transformation)
  }
}


#' Write the taxa ages
#'
#' @description
#' Writes the ages of the specimen in the true tree to a file. The tsv format used
#' here is directly compatible with RevBayes
#'
#' @param data Morpho object
#' @param file File name
#' @param uncertainty Numeric. Adds uncertainty to fossil ages in the morpho object.
#'  The ages in the object are point estimates by default; setting `uncertainty`
#'  will create an age range of ± this value (in millions of years).
#'
#'
#' @export

write.tsv <- function (data, file, uncertainty = 0) {

  ## ages of full tree
  tip_depths <- ape::node.depth.edgelength(data$trees$TimeTree)[1:length(data$trees$TimeTree$tip.label)]
  tree_height <- max(node.depth.edgelength(data$trees$TimeTree))
  tip_ages <-   round(abs(tip_depths - tree_height),3)
 # extant_tips <- data$trees$TimeTree$tip.label[abs(tip_depths - tree_height) < 1e-8]

  cat("taxon", "min_age", "max_age", sep = "\t", "\n", file = file)
  for ( i in 1:length(tip_ages)){
    if (tip_ages[i] == 0){
    cat(data$trees$TimeTree$tip.label[i], tip_ages[i] ,
        tip_ages[i], sep = "\t", file = file, append = T )
      cat("\n", file = file, append = TRUE)
    } else {
      cat(data$trees$TimeTree$tip.label[i], (tip_ages[i] -  uncertainty) ,
          (tip_ages[i] + uncertainty), sep = "\t", file = file, append = T )
      cat("\n", file = file, append = TRUE)
    }
  }

}

#' Write the taxa ages of reconstructed tree
#'
#' @description
#' Writes the ages of the specimen in the reconstructed tree to a file. The tsv format used
#' here is directly compatible with RevBayes
#'
#' @param data Morpho object
#' @param file File name
#' @param uncertainty Numeric. Adds uncertainty to fossil ages in the morpho object.
#'  The ages in the object are point estimates by default; setting `uncertainty`
#'  will create an age range of ± this value (in millions of years).

write.recon.tsv <- function (data, file, uncertainty = 0){

  r_tree <- FossilSim::reconstructed.tree.fossils.objects(fossils  = data$fossil,
                                                          tree = data$trees$TimeTree)
  transformations <- write.recon.matrix(data, keep_matrix = T)

  cat("taxon", "min_age", "max_age", sep = "\t", "\n", file = file)


  ## true tree tips
  tps <- unname(r_tree$tree$tip.label)
  matches <- grepl("_1$", tps )
  # Extract elements that match
  reconTreeTips <- gsub("_1$", "", tps[matches])

  ## add these tip labels to the file + plus all sampled ancestor
  seq_tips <- which(names(data$sequences$tips) %in% reconTreeTips)


  for ( i in 1:length(reconTreeTips)){
  ord <-  which(data$trees$TimeTree$tip.label ==reconTreeTips[i])
  node_pos <- ape::node.depth.edgelength(data$trees$TimeTree)[ord]
  tree_height <- max(node.depth.edgelength(data$trees$TimeTree))
  tip_ages <-   round(abs(node_pos - tree_height),3)

  if(tip_ages == 0){
    nm <- paste0(reconTreeTips[i], "_1")
    cat(nm, tip_ages, tip_ages,
        sep = "\t", file = file, append = T )
    cat("\n", file = file, append = TRUE)
  } else {
    nm <- paste0(reconTreeTips[i], "_1")
    if (tip_ages - uncertainty < 0){
      min_age <- 0
    } else {
      min_age <- tip_ages - uncertainty
    }
    cat(nm,min_age, (tip_ages + uncertainty),
        sep = "\t", file = file, append = T )
    cat("\n", file = file, append = TRUE)

  }

  }

  ## sampled ancestors

  if (length(transformations[,"Morphosim"]) > 0){
  for (i in 1:length(transformations[,1])){

    parts <- as.numeric(strsplit(transformations[i, "Morphosim"], "_")[[1]])
    specimen_num <- parts[1]
    branch_num   <- parts[2]

    # Subset the data frame to get hmin
    hmin <- data$fossil$hmin[data$fossil$ape.branch == branch_num &
        data$fossil$specimen  == specimen_num
    ]


    nm <- paste0(transformations[i, "Fossilsim"], "_2")
    if (hmin - uncertainty < 0){
      min_age <- 0
    } else {
      min_age <- hmin - uncertainty
    }
    cat(nm,min_age, (hmin + uncertainty),
        sep = "\t", file = file, append = T )
    cat("\n", file = file, append = TRUE)
  }
  }

}







