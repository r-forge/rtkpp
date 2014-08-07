pkgname <- "rtkpp"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('rtkpp')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("ClusteringAlgo-class")
### * ClusteringAlgo-class

flush(stderr()); flush(stdout())

### Name: ClusteringAlgo
### Title: ['ClusteringAlgo'] class for clustering algorithms.
### Aliases: ClusteringAlgo ClusteringAlgo-class

### ** Examples

new("ClusteringAlgo")
new("ClusteringAlgo", algo="SEM", nbIteration=1000)
getSlots("ClusteringAlgo")



cleanEx()
nameEx("ClusteringCategoricalModel-class")
### * ClusteringCategoricalModel-class

flush(stderr()); flush(stdout())

### Name: ClusteringCategoricalModel-class
### Title: Definition of the ['ClusteringCategoricalModel'] class
### Aliases: ClusteringCategoricalModel-class

### ** Examples

getSlots("ClusteringCategoricalModel")
  new("ClusteringCategoricalModel", data=iris[1:4])



cleanEx()
nameEx("ClusteringDiagGaussianModel-class")
### * ClusteringDiagGaussianModel-class

flush(stderr()); flush(stdout())

### Name: ClusteringDiagGaussianModel
### Title: Definition of the ['ClusteringDiagGaussianModel'] class
### Aliases: ClusteringDiagGaussianModel ClusteringDiagGaussianModel-class

### ** Examples

getSlots("ClusteringDiagGaussianModel")
  new("ClusteringDiagGaussianModel", data=iris[1:4])



cleanEx()
nameEx("ClusteringGammaModel-class")
### * ClusteringGammaModel-class

flush(stderr()); flush(stdout())

### Name: ClusteringGammaModel
### Title: Definition of the ['ClusteringGammaModel'] class
### Aliases: ClusteringGammaModel ClusteringGammaModel-class

### ** Examples

getSlots("ClusteringGammaModel")
  new("ClusteringGammaModel", data=iris[1:4])



cleanEx()
nameEx("ClusteringInit-class")
### * ClusteringInit-class

flush(stderr()); flush(stdout())

### Name: ClusteringInit
### Title: ['ClusteringInit'] class
### Aliases: ClusteringInit ClusteringInit-class

### ** Examples

getSlots("ClusteringInit")
  new("ClusteringInit")



cleanEx()
nameEx("ClusteringStrategy-class")
### * ClusteringStrategy-class

flush(stderr()); flush(stdout())

### Name: ClusteringStrategy
### Title: Constructor of ['ClusteringStrategy'] class
### Aliases: ClusteringStrategy ClusteringStrategy-class

### ** Examples

new("ClusteringStrategy")
  new("ClusteringStrategy", shortAlgo=clusteringAlgo("SEM",1000))
  getSlots("ClusteringStrategy")



cleanEx()
nameEx("IClusteringModel-class")
### * IClusteringModel-class

flush(stderr()); flush(stdout())

### Name: IClusteringModel
### Title: Interface Class ['IClusteringModel'] for clustering models.
### Aliases: IClusteringModel IClusteringModel-class

### ** Examples

getSlots("IClusteringModel")



cleanEx()
nameEx("IClusteringModelNames-class")
### * IClusteringModelNames-class

flush(stderr()); flush(stdout())

### Name: IClusteringModelNames
### Title: Interface base class ['IClusteringModelNames']
### Aliases: IClusteringModelNames IClusteringModelNames-class

### ** Examples

getSlots("IClusteringModelNames")



cleanEx()
nameEx("clusteringAlgo")
### * clusteringAlgo

flush(stderr()); flush(stdout())

### Name: clusteringAlgo
### Title: Create an instance of ['ClusteringAlgo'] class
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
### Title: Create an instance of ['ClusteringInit'] class
### Aliases: clusteringInit

### ** Examples

clusteringInit(method = "class", nbInit=1, algo="CEM",nbIteration=50, epsilon=0.00001)
clusteringInit(nbIteration=0) # no algorithm



cleanEx()
nameEx("clusteringStrategy")
### * clusteringStrategy

flush(stderr()); flush(stdout())

### Name: clusteringStrategy
### Title: Create an instance of ['ClusteringStrategy'] class
### Aliases: clusteringStrategy

### ** Examples

clusteringStrategy()
   clusteringStrategy(longRunAlgo= "CEM", nbLongIteration=100)
   clusteringStrategy(nbTry = 3, nbInit= 1, shortRunAlgo= "SEM", nbShortIteration=100)



cleanEx()
nameEx("print-methods")
### * print-methods

flush(stderr()); flush(stdout())

### Name: print,ClusteringAlgo-method
### Title: Print a rtkpp class to standard output.
### Aliases: print print,ClusteringAlgo-method print,ClusteringInit-method
###   print,ClusteringStrategy-method print,IClusteringModelNames-method
###   print-algo,ClusteringAlgo,ClusteringAlgo-method
###   print-init,ClusteringInit,ClusteringInit-method
###   print-strategy,ClusteringStrategy,ClusteringStrategy-method
###   print-strategy,IClusteringModelNames-method

### ** Examples

## for strategy
  strategy <- clusteringStrategy()
  print(strategy)



cleanEx()
nameEx("show-methods")
### * show-methods

flush(stderr()); flush(stdout())

### Name: show,ClusteringAlgo-method
### Title: Show description of a rtkpp class to standard output.
### Aliases: show show,ClusteringAlgo-method show,ClusteringInit-method
###   show,ClusteringStrategy-method show,IClusteringModelNames-method
###   show-algo,ClusteringAlgo,ClusteringAlgo-method
###   show-init,ClusteringInit,ClusteringInit-method
###   show-strategy,ClusteringStrategy,ClusteringStrategy-method
###   show-strategy,IClusteringModelNames-method

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
