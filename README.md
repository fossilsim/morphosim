## MorphoSim 

R package for simulating discrete character data along rooted phylogenies. 

The latest version  can be installed in R using the following commands

    library(devtools)
    devtools::install_github("https://github.com/fossilsim/morphosim")
    
    
### Quick Start 
To simulate traits at internal nodes 

```
phy <- ape::rtree(10)
# simulate characters at internal nodes of the tree
traits <- sim.morpho(phy, k=4, trait.num =12)
```

For a general introduction see [here](https://github.com/fossilsim/morphosim/blob/main/intro.html)


