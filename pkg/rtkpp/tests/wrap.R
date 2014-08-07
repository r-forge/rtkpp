# TODO: Add comment
# 
# Author: iovleff
###############################################################################
library(rtkpp)
data(iris)
mat <- as.matrix(iris[1:4])[1:10,]
mat[10,4] = NA;
xem <- .Call("wrapper", mat, PACKAGE="rtkpp")
xem
mat
