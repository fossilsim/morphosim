## MorphoSim 

R package for simulating discrete character data along rooted phylogenies. 

The latest version (when made public) can be installed in R using the following commands

    library(devtools)
    install_github("fossilsim/morphosim")
    
    
### Quick Start 
To simulate traits at internal nodes 

```
phy <- ape::rtree(10)
# simulate characters at internal nodes of the tree
traits <- sim.morpho(phy, k=4, trait.num =12)
```


To simulate the continuous evolution of characters along the branches

```
# simulate characters along the branches of the tree
traits <- sim.morpho.completeprocess(phy, k=4, trait.num =12)
```


