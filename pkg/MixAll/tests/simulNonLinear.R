simulNonLinear <- function(n, var = 0.25)
{
  # simul class
  z <- sapply(runif(n), FUN= function(x){ return(floor(x+0.5))})
  tau <- runif(n, -4,4)
  eta <- rnorm(n, 0, sd=sqrt(var))
  x<- matrix(ncol=2, nrow = n);
  for (i in 1:n)
  {
    if (z[i] == 0)
    { x [i,] <- t(c(-1 + tau[i]+eta[i],  8 - tau[i]*tau[i]/2 + eta[i])); }
    else
    { x [i,] <- t(c( 1 + tau[i]+eta[i], -8 + tau[i]*tau[i]/2 + eta[i]));}
  }
  list(x=x, z=z)
}
