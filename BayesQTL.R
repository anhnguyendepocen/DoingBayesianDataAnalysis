## Multiple Interval Mapping
## http://www.genetics.org/content/152/3/1203.long



for ( s in 1: 25000 )
{
  ## update x
  lpy.x <- NULL; 
  for( x in 1:100 ){ 
    lpy.x <- c(lpy.x, lpy.theta(y, G, x ,mu, s2 ) )
  }
  
  x <- sample( 1:100, 1 , prob = exp(lpy.x - max(lpy.x) ) )

    ## update gx
  pg1.x <- prhet.sG( x, G, mpos )
  py.g1 <- dnorm ( y, mu[2], sqrt( s2 ) )
  py.g0 <- dnorm ( y, mu[1], sqrt( s2 ) )
  pg1.yx <- py.g1 * pg1.x / ( py.g1 * pg1.x + py.g0 * (1 - pg1.x ) )
  gx <- rbinom (n, 1, pg1.yx )
  
  ## update s2
  s2 <- 1/rgamma ( 1, (nu0 + n ) / 2, (nu0 * s20 + sum ( (y - mu[gx + 1] )^2) ) / 2 )
  
  ## update mu
  mu <- rnorm ( 2, ( mu0 * k0 + tapply(y, gx, sum) ) / ( k0 + table(gx) ), sqrt(s2 / (k0 + table(gx))) )
}