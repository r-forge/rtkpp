# TODO: Add comment
#
# Author: iovleff
###############################################################################
library(rtkpp)
data(iris)
mat1 <- as.matrix(iris[1:4])[1:10,]
mat2 <- as.matrix(iris[1:4])[11:20,]
mat3 <- as.matrix(iris[1:4])[21:30,]
mat1[10,4] = NA;
mat2[10,4] = NA;
mat3[10,4] = NA;

xem <- .Call("wrapper", mat1, mat2, mat3, PACKAGE="rtkpp")

mat1
mat2
mat3
xem

