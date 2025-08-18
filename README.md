## MorphoSim 

R package for simulating discrete character data along rooted phylogenies. By integrating with established R packages, MorphoSim enables the simulation of biologically meaningful datasets across a variety of scenarios. 

The latest version  can be installed in R using the following commands

``` r   
library(devtools)
    devtools::install_github("https://github.com/fossilsim/morphosim")
```
    
### Example workflow for data simulation using morphoSim 

##### Step 1: Simulate a birth death tree using Tree sim.
```r
lambda = 0.1
mu = 0.05
tips = 5
t = TreeSim::sim.bd.taxa(n = tips, numbsim = 1, 
                         lambda = lambda, mu = mu)[[1]]
 plot(t)                        
```
                         

 ##### Step 2: Simulate fossil/extant sampling along the tree using fossil sim
 
```r
rate = 0.3
f = FossilSim::sim.fossils.poisson(rate = rate, tree = t, root.edge = F)
rho = 0.3
f2 = FossilSim::sim.extant.samples(fossils = f, tree = t, rho = rho)
```

##### Step 3: Simulate branch lengths using simclock. Use the following code to install simclock

```r
devtools::install_github("dosreislab/simclock")

```
Simulate branch lengths using a independent log normal relaxed clock model

```r
reltt <- simclock::relaxed.tree(t, model="iln", r=.04e-1, s2=.26e-2)
```

##### Step 4: Simulate Morphological data using morphosim. 

```r
morpho_data <-  sim.morpho(k = (3),
                time.tree = t,
                trait.num = 10,
                ancestral = TRUE,
                br.rates = reltt$edge.length,
                #partition = c(10,10),
                ACRV = "gamma",
                fossil = f2,
                variable = TRUE,
                ncats.gamma = 4)

```
### Further steps
##### Simulate missing data

```r
missing.morpho.data <- sim.missing.data(data = morpho_data, seq = "tips", 
                 method = "random", probability = 0.1)
```

##### Plotting data

```r
plot(morpho_data, trait= 10, timetree = T,
     fossil = T, root.edge = F, reconstructed = T)

```
##### Summaries data

```r
sum <- stats_morpho(data = morpho_data)
sum$Tree
unique(sum$Convergent_Traits$trait)
```

### Writing out simulated data for phylogenetic inference.  

```r
ape::write.nexus.data(morpho_data$sequences$tips, file = "../tip.nex")

write.recon.tree(data = morpho_data, file = "../example.txt")
write.recon.matrix(data = morpho_data, file = "../example_mat.nex")
write.recon.tsv(data = morpho_data, file = "../example.tsv")
```
