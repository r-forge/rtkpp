pkgname <- "rtkpp"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('rtkpp')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("ClusterAlgo-class")
### * ClusterAlgo-class

flush(stderr()); flush(stdout())

### Name: ClusterAlgo
### Title: ['ClusterAlgo'] class for clustering algorithms.
### Aliases: ClusterAlgo ClusterAlgo-class

### ** Examples

new("ClusterAlgo")
new("ClusterAlgo", algo="SEM", nbIteration=1000)
getSlots("ClusterAlgo")



cleanEx()
nameEx("ClusterCategoricalModel-class")
### * ClusterCategoricalModel-class

flush(stderr()); flush(stdout())

### Name: ClusterCategoricalModel-class
### Title: Definition of the ['ClusterCategoricalModel'] class
### Aliases: ClusterCategoricalModel-class

### ** Examples

getSlots("ClusterCategoricalModel")
  new("ClusterCategoricalModel", data=iris[1:4])



cleanEx()
nameEx("ClusterDiagGaussianModel-class")
### * ClusterDiagGaussianModel-class

flush(stderr()); flush(stdout())

### Name: ClusterDiagGaussianModel
### Title: Definition of the ['ClusterDiagGaussianModel'] class
### Aliases: ClusterDiagGaussianModel ClusterDiagGaussianModel-class

### ** Examples

getSlots("ClusterDiagGaussianModel")
  new("ClusterDiagGaussianModel", data=iris[1:4])



cleanEx()
nameEx("ClusterGammaModel-class")
### * ClusterGammaModel-class

flush(stderr()); flush(stdout())

### Name: ClusterGammaModel
### Title: Definition of the ['ClusterGammaModel'] class
### Aliases: ClusterGammaModel ClusterGammaModel-class

### ** Examples

getSlots("ClusterGammaModel")
  new("ClusterGammaModel", data=iris[1:4])



cleanEx()
nameEx("ClusterInit-class")
### * ClusterInit-class

flush(stderr()); flush(stdout())

### Name: ClusterInit
### Title: ['ClusterInit'] class
### Aliases: ClusterInit ClusterInit-class

### ** Examples

getSlots("ClusterInit")
  new("ClusterInit")



cleanEx()
nameEx("ClusterStrategy-class")
### * ClusterStrategy-class

flush(stderr()); flush(stdout())

### Name: ClusterStrategy
### Title: Constructor of ['ClusterStrategy'] class
### Aliases: ClusterStrategy ClusterStrategy-class

### ** Examples

new("ClusterStrategy")
  new("ClusterStrategy", shortAlgo=clusteringAlgo("SEM",1000))
  getSlots("ClusterStrategy")



cleanEx()
nameEx("IClusterModel-class")
### * IClusterModel-class

flush(stderr()); flush(stdout())

### Name: IClusterModel
### Title: Interface Class ['IClusterModel'] for clustering models.
### Aliases: IClusterModel IClusterModel-class

### ** Examples

getSlots("IClusterModel")



cleanEx()
nameEx("IClusterModelNames-class")
### * IClusterModelNames-class

flush(stderr()); flush(stdout())

### Name: IClusterModelNames
### Title: Interface base class ['IClusterModelNames']
### Aliases: IClusterModelNames IClusterModelNames-class

### ** Examples

getSlots("IClusterModelNames")



cleanEx()
nameEx("clusteringAlgo")
### * clusteringAlgo

flush(stderr()); flush(stdout())

### Name: clusteringAlgo
### Title: Create an instance of ['ClusterAlgo'] class
### Aliases: clusteringAlgo

### ** Examples

clusteringAlgo()
clusteringAlgo(algo="SEM", nbIteration=50)
clusteringAlgo(algo="CEM", epsilon = 1e-06)



cleanEx()
nameEx("clusteringInit")
### * clusteringInit

flush(stderr()); flush(stdout())

### Name: clusteringInit
### Title: Create an instance of ['ClusterInit'] class
### Aliases: clusteringInit

### ** Examples

clusteringInit(method = "class", nbInit=1, algo="CEM",nbIteration=50, epsilon=0.00001)
clusteringInit(nbIteration=0) # no algorithm



cleanEx()
nameEx("clusteringStrategy")
### * clusteringStrategy

flush(stderr()); flush(stdout())

### Name: clusteringStrategy
### Title: Create an instance of ['ClusterStrategy'] class
### Aliases: clusteringStrategy

### ** Examples

clusteringStrategy()
   clusteringStrategy(longRunAlgo= "CEM", nbLongIteration=100)
   clusteringStrategy(nbTry = 3, nbInit= 1, shortRunAlgo= "SEM", nbShortIteration=100)



cleanEx()
nameEx("print-methods")
### * print-methods

flush(stderr()); flush(stdout())

### Name: print,ClusterAlgo-method
### Title: Print a rtkpp class to standard output.
### Aliases: print print,ClusterAlgo-method print,ClusterInit-method
###   print,ClusterStrategy-method print,IClusterModelNames-method
###   print-algo,ClusterAlgo,ClusterAlgo-method
###   print-init,ClusterInit,ClusterInit-method
###   print-strategy,ClusterStrategy,ClusterStrategy-method
###   print-strategy,IClusterModelNames-method

### ** Examples

## for strategy
  strategy <- clusteringStrategy()
  print(strategy)



cleanEx()
nameEx("show-methods")
### * show-methods

flush(stderr()); flush(stdout())

### Name: show,ClusterAlgo-method
### Title: Show description of a rtkpp class to standard output.
### Aliases: show show,ClusterAlgo-method show,ClusterInit-method
###   show,ClusterStrategy-method show,IClusterModelNames-method
###   show-algo,ClusterAlgo,ClusterAlgo-method
###   show-init,ClusterInit,ClusterInit-method
###   show-strategy,ClusterStrategy,ClusterStrategy-method
###   show-strategy,IClusterModelNames-method

### ** Examples

## for strategy
  strategy <- clusteringStrategy()
  show(strategy)



### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
