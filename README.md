## MorphoSim 

R package for simulating discrete character data along rooted phylogenies. Through integration with established R packages, MorphoSim enables the simulation of biologically meaningful datasets across a variety of scenarios. 

The latest version can be installed from github using `devtools` 

``` r   
library(devtools)
devtools::install_github("https://github.com/fossilsim/morphosim")
```
    


## Ex 1: Simulating Data Along a Phylogenetic Tree - Basic set up
To simualte morphological data you will first need a phylogenetic tree. You can simulate a tree using the existing R package `Treesim`

```r
install.packages("Treesim")
library(Treesim)
```

##### Simulate a birth death tree using Treesim
```r
lambda = 0.1         # speciation rate
mu = 0.05            # exintction rate
tips = 5             # number of tips

tree = TreeSim::sim.bd.taxa(n = tips, 
                         numbsim = 1, 
                         lambda = lambda, 
                         mu = mu)[[1]]                       
```
##### Simulate Morphological data using morphosim. 
MorphoSim lets you specify a variety of parameters for your simulations. In this example, we will focus on a few key parameters to get started.
- **k**: the maximum number of character states. This can be a single integer greater than 1 or a vector of integers if you want to simulate more than 1 partitions. So in the following example, we are simulating two partitions, the first with 2 states and the second with 3.
- **timetree**: if using a time tree specify here.
- **tree**: if using a genetic distance tree specify here.
- **partitions**: specify the number of traits per partition. This must match the number of states specificed in k.
- **triat.num**: specify the total number of traits 
- **br.rates**: This can be a single integer (strict clock) or a vector of integers (relaxed clock). This must be provided when using a time tree to convert the branch lengths into genetic distance.
- **ancestral**: return the states at the internal nodes in the tree.

> Setting up a simulation as follows will simulate characters under an Mk (equal rates) model with a strict clock.


```r
morpho_data <-  sim.morpho(k = c(2,3), 
                           time.tree = tree,
                           tree = NULL,
                           partition = c(10,10),
                           trait.num = 20,
                           br.rates = 0.1,
                           ancestral = TRUE)  

```
##### Explore the simluated data
Morphosim has a plotting function which allows you to plot the simualted data along the branches of the tree. 
```r
plot(morpho_data, 
     trait= 1) # specify which trait you want to plot along the tree
```
This function returns a number of characterstics about your data. It calculates the consistency index, the rentention index, identifies convergent traits, and the number of extinct and extant taxa.
```r
stats <- stats_morpho(morpho_data)
```




## Example 2: Simulating Biologically Realistic Data Sets
The previous example, simulated characters according to an Mk model under a strict clock. However, there are more complex models we may want to use for simulations. Firstly, with MorphoSim we can use a relaxed clock model where each branch has a different rate. We can use an existing R package, `simclock` to simulate rates along our tree which we can input into Morphosim 


#####  Simulate branch lengths using a independent log normal relaxed clock model

```r

devtools::install_github("dosreislab/simclock")
reltt <- simclock::relaxed.tree(t, 
                                model="iln", 
                                r=.04e-1,     # mean mutation rate
                                s2=.26e-2).   # diffusion rate
```
For the character data we can relax the assumptions of the Mk model in a number of ways as described below

-**ACRV**: Here you can specify if you would like to model among character rate variation. This is commonly used for morphological data sets as traits may evolve at different rates. There are three input options here, `gamma`, `lgn`, and `user`. Gamma uses a discrete gamma distribution, lgn will use a discrete lognorma, and user allows the user to specify rates.

-**ACRV.ncats**: the number of rate categories you want to simualte under

-**alpha.gamma**: the shape of the alpha distribution

-**variable**: When set to true this will only return traits which vary acrosss taxa, i.e., the MkV model

```r
morpho_data <-  sim.morpho(k = c(2,3), 
                           time.tree = tree,
                           tree = NULL,
                           partition = c(10,10),
                           trait.num = 20,
                           br.rates = reltt$edge.length,
                           ACRV = "gamma",
                           alpha.gamma = 1,
                           ACRV.ncats = 4
                           variable = TRUE)

```

> Note: If you choose to simulate data under an ACRV using a lognormal distribution, you will need to provide values for the `meanlog` and the `sdlog`. Similarly for the user option the rates should be passed through `define.ACRV.rates`

### Ordered characters
In morphological datasets, transitions between character states are often restricted, meaning that a change can only occur between specific states. For example, a character might be able to transition from 0 → 1, but not directly from 1 → 2; the species must pass through the intermediate state. These are known as correlated characters. In morphosim, we can simulate such characters by defining a Q matrix that enforces these allowed transitions.

Define a Q matrix the allows 0 &harr; 1 &harr; 2

```r
ord_Q <- matrix(c(
  -0.5, 0.5, 0.0,
  0.3333333, -0.6666667, 0.3333333,
  0.0, 0.5, -0.5
), nrow = 3, byrow = TRUE)r
```
This Q matrix can then be used to simulate characters

```r
morpho_data <-  sim.morpho(k = 3, 
                           define.Q = ord_Q
                           time.tree = tree,
                           tree = NULL,
                           trait.num = 10,
                           br.rates = reltt$edge.length,
                           ACRV = "gamma",
                           alpha.gamma = 1,
                           ACRV.ncats = 4
                           variable = TRUE)
```
