## Monte Carlo Simulation Methods Using R Book ##

library(mcsm)
library(datasets)
rm(list = ls()); gc(reset = TRUE)
x <- rnorm(1)
for( t in 2:10^3 ) {
  x <- c(x,.09*x[t-1]+rnorm(1))
}
plot(x,type="l",xlab="time",ylab="x",lwd=2,lty=2, col="steelblue",ylim=range(cumsum(x)))
lines(cumsum(x),lwd=2,col="orange3")

runif(100, min=2, max=5)


Nsim <- 10^4 #number of random numbers
x <- runif(Nsim)
x1 <- x[-Nsim] #vectors to plot
x2 <- x[-1] #adjacent pairs
par(mfrow = c(1, 3))
hist(x)
plot(x1,x2)
acf(x)



Nsim <- 10^4 #number of random variables
U <- runif(Nsim)
X <- -log(U) #transforms of uniforms
Y <- rexp(Nsim) #exponentials from R
par(mfrow=c(1,2)) #plots
hist(X, freq = FALSE, main = "Exp from Uniform")
hist(Y, freq = FALSE, main = "Exp from R")



  # Example 2.2: Generate Chi-Squre with 6 df

U <- runif(3*10^4)
U <- matrix(data = U, nrow = 3) #matrix for sums
X <- -log(U) #uniform to exponential
X <- 2* apply(X, 2, sum) #sum up to get chi square


  # Example 2.5 Generate a Poisson from Lambda

Nsim <- 10^4; lambda <- 100
spread <- 3*sqrt(lambda)
t <- round(seq(max(0, lambda - spread), lambda + spread, 1))
prob <- ppois(t, lambda)
X <- rep(0, Nsim)
for( i in 1 : Nsim ) {
  u <- runif(1)
  X[i] <- t[1] + sum(prob < u) 
}



  # Section 2.2.3 Mixture Representations

  # Example 2.6 Negative Binomial has Mixture representation between gamma and poisson

Nsim <- 10^4
n <- 6; p <- 0.3
y <- rgamma(Nsim, n, rate = p/(1-p))
x <- rpois(Nsim, y)
hist(x, main = "", freq = F, col = "grey", breaks = 40)
lines(1 : 50, dnbinom(1 : 50, n, p), lwd = 2,col = "sienna")



  # Accept Reject Method
    # randg(1) is a 'made-up' function as a place holder of an R function represnetating a distribution

u <- runif(1) * M
y <- randg(1)
while ( u > f(y) / g(y) ) {
  u <- runif(1) * M
  y <- randg(1)
}



  # Example 3.1

ch <- function( la ) {
  integrate(function( x ) {x ^ (la - 1) * exp(-x) }, 0, Inf)$val
}
plot(lgamma(seq(0.01, 10, le = 100)),log(apply(as.matrix(
  seq(0.01, 10, le = 100)), 1, ch)), xlab = "log(integrate(f))",
  ylab = expression(log(Gamma(lambda))), pch = 19, cex = 0.6)



  # Example 3.3: Monte Carlo Integrating h(x) = [cos(50x) + sin(20x)]^2

h <- function(x) { (cos(50 * x) + sin(20 * x))^2 }
par(mar = c(2, 2, 2, 1), mfrow = c(2, 1))
curve(h, xlab = "Function", ylab = "", lwd = 2)
integrate(h, 0, 1)
# 0.965201 with absolute error < 1.9e-10
x <- h(runif(10^4))
estint <- cumsum(x)/(1 : 10^4)
esterr <- sqrt(cumsum((x - estint)^2))/(1 : 10^4)
plot(estint, xlab = "Mean and error range", type = "l", lwd = 2,
     ylim = mean(x) + 20 * c(-esterr[10^4], esterr[10^4]), ylab = "")
lines(estint + 2 * esterr, col = "gold", lwd = 2)
lines(estint - 2 * esterr, col = "gold", lwd = 2)



  # Example 3.5 Importance Sampling

Nsim=10^3
y=rexp(Nsim)+4.5
weit=dnorm(y)/dexp(y-4.5)
plot(cumsum(weit)/1:Nsim,type="l")
abline(a=pnorm(-4.5),b=0,col="red")


# Example 3.6 Importance Sampling
par(mfrow = c(1,2))
f=function(a,b){
  exp(2*(lgamma(a+b)-lgamma(a)-lgamma(b))+
        a*log(.3)+b*log(.2))}
aa=1:150 #alpha grid for image
bb=1:100 #beta grid for image
post=outer(aa,bb,f)
image(aa,bb,post,xlab=expression(alpha),ylab=" ")
contour(aa,bb,post,add=T)

x=matrix(rt(2*10^4,3),ncol=2) #T sample
E=matrix(c(220,190,190,180),ncol=2) #Scale matrix
image(aa,bb,post,xlab=expression(alpha),ylab=" ")
y=t(t(chol(E))%*%t(x)+c(50,45))
points(y,cex=.6,pch=19)


ine=apply(y,1,min)
y=y[ine>0,]
x=x[ine>0,]
normx=sqrt(x[,1]^2+x[,2]^2)
f=function(a) exp(2*(lgamma(a[1]+a[2])-lgamma(a[1])-lgamma(a[2]))+a[1]*log(.3)+a[2]*log(.2))
h=function(a) exp(1*(lgamma(a[1]+a[2])-lgamma(a[1])-lgamma(a[2]))+a[1]*log(.5)+a[2]*log(.5))
den=dt(normx,3)
mean(f(y)/den)/mean(h(y)/den)

mean(y[,1]*apply(y,1,f)/den)/mean(apply(y,1,h)/den)
mean(y[,2]*apply(y,1,f)/den)/mean(apply(y,1,h)/den)



par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
weit = (apply(y, 1, f)/den)/mean(apply(y, 1, h)/den)
image(aa, bb, post, xlab = expression(alpha), ylab = expression(beta))
points(y[sample(1:length(weit), 10^3, rep = TRUE, pro = weit),], cex = 0.6, pch = 19)
boxplot(weit, ylab = "importance weight")
plot(cumsum(weit)/(1:length(weit)), type = "l", xlab = "simulations", ylab = "marginal likelihood")
boot = matrix(0, ncol = length(weit), nrow = 100)
for( t in 1:100 ) boot[t, ] = cumsum(sample(weit))/(1:length(weit))
uppa = apply(boot, 2, quantile, 0.95)
lowa = apply(boot, 2, quantile, 0.05)
polygon(c(1 : length(weit), length(weit) : 1), c(uppa, rev(lowa)), col = "gold")
lines(cumsum(weit)/(1 : length(weit)), lwd = 2)
plot(cumsum(weit)^2/cumsum(weit^2), type = "l", xlab = "simulations", ylab = "Effective sample size", lwd = 2)

  # Example 3.8 Importance Sampling, Infinite Variance

x=rnorm(10^6)
wein=dcauchy(x)/dnorm(x)
boxplot(wein/sum(wein))
plot(cumsum(wein*(x>2)*(x<6))/cumsum(wein),type="l")
abline(a=pcauchy(6)-pcauchy(2),b=0,col="sienna")



sam1=rt(.95*10^4,df=2)
sam2=1+.5*rt(.05*10^4,df=2)^2
sam=sample(c(sam1,sam2),.95*10^4)
weit=dt(sam,df=2)/(0.95*dt(sam,df=2)+.05*(sam>0)*dt(sqrt(2*abs(sam-1)),df=2)*sqrt(2)/sqrt(abs(sam-1)))
plot(cumsum(h(sam1))/(1:length(sam1)),ty="l")
lines(cumsum(weit*h(sam))/1:length(sam1),col="blue")



glm(type ~ bmi,data=Pima.tr,family=binomial(link="probit"))

post=function(beda){
  mia=mean(Pima.tr$bmi)
  prod(pnorm(beda[1]+(Pima.tr$bm[Pima.tr$t=="Yes"]-
                        mia)*beda[2]))*
    prod(pnorm(-beda[1]-(Pima.tr$bm[Pima.tr$t=="No"]
                         -mia)*beda[2]))/exp(sum(beda^2)/200)
}

sim <- cbind(rnorm(10^3, mean = -0.4, sd = 0.04), rnorm(10^3, mean = 0.065, sd = 0.005))
weit <- apply(sim, 1, post)/(dnorm(sim[, 1], mean = -0.4, sd = 0.04) * dnorm(sim[, 2], mean = 0.065, sd = 0.005))


sim=rbind(sim[1:(.95*10^3),],cbind(rnorm(.05*10^3,sd=10), rnorm(.05*10^3,sd=10)))
weit=apply(sim,1,post)/(.95*dnorm(sim[,1],m=-.4,sd=.081)*dnorm(sim[,2],m=0.065,sd=.01) + 
                            0.05*dnorm(sim[,1],sd=10)*dnorm(sim[,2],sd=10))


  # Example 4.1 Monitoring Variation

  # h function from 3.3
h <- function(x) { (cos(50 * x) + sin(20 * x))^2 }

x <- matrix(h(runif(200 * 10^4)), ncol = 200)
estint <- apply(x, 2, cumsum)/(1 : 10^4)

plot(estint[,1],ty = "l", col = 0, ylim = c(0.8, 1.2))
y <- apply(estint,1, quantile, c(.025,.975))
polygon(c(1:10^4,10^4:1), c(y[1,], rev(y[2,])), col = "wheat")

boot <- matrix(sample(x[,1],200*10^4, rep = TRUE), nrow = 10^4, ncol = 200)
bootit <- apply(boot, 2, cumsum)/(1 : 10^4)
bootup <- apply(bootit, 1, quantile, 0.975)
bootdo <- apply(bootit, 1, quantile, 0.025)


  # Example 4.2 / 4.3: Effective Sample Size and Perplexity

norma <- matrix(rnorm(500*10^4), ncol = 500) + 2.5
weit <- 1/(1+norma^2)
esti <- apply(norma * weit, 2, cumsum)/apply(weit, 2, cumsum)
plot(esti[,1], type="l", col = "white", ylim = c(1.7, 1.9))
band <- apply(esti, 1, quantile, c(0.025, 0.975))
polygon(c(1 : 10^4, 10^4 : 1), c(band[1, ], rev(band[2, ])))

vare=cumsum(weit[,1]*norma[,1]^2)/cumsum(weit[,1])-esti[,1]^2
lines(esti[,1]+2*sqrt(vare/(1:10^4)))
lines(esti[,1]-2*sqrt(vare/(1:10^4)))

cocha=matrix(rcauchy(500*10^4),ncol=500)
range(cocha)
wach=dnorm(cocha,mean=2.5)
range(wach)
wachd=wach
wachd[apply(wachd,2,cumsum)<10^(-10)]=10^(-10)

ess=apply(weit,2,cumsum)^2/apply(weit^2,2,cumsum)
essbo=apply(ess,1,quantile,c(.025,.975))
ech=apply(wachd,2,cumsum)^2/apply(wachd^2,2,cumsum)
echbo=apply(ech,1,quantile,c(.025,.975))

sumweit=apply(weit,2,cumsum)
plex=(apply(weit*log(weit),2,cumsum)/sumweit)-log(sumweit)
chumweit=apply(wachd,2,cumsum)
plech=(apply(wachd*log(wachd),2,cumsum)/chumweit) - log(chumweit)
plob=apply(exp(plex),1,quantile,c(.025,.975))
ploch=apply(exp(plech),1,quantile,c(.025,.975))



## Example 4.5

Nsim <- 10^4
norma <- rnorm(Nsim)+2.5
hnorm <- norma*dcauchy(norma)
munorm <- mean(hnorm)
sdnorm <- sd(hnorm)
f <- function(x) (cumsum(hnorm))[round(Nsim * x)]/round(x * Nsim)
curve(munorm + (0.1 + 3.15 * sqrt(x)) * sdnorm * 10^2/round(x * Nsim), lwd = 2, from = 0,to = 1)
curve(munorm - ((0.1 + 3.15 * sqrt(x)) * sdnorm * 10^2/round(x * Nsim)), lwd = 2,from = 0,to = 1,add = T)
curve(f,lwd=2,from=0.001,to=1,col="steelblue", add = T)

norma <- rcauchy(Nsim)
hnorm <- norma * dnorm(norma-2.5)


## Example 4.6: Rao-Blackwellization and deconditioning (Dickey's decomposition)

nu <- 5
mu <- 3
sigma <- 0.5

y <- sqrt(rchisq(Nsim,df = nu)/nu)
x <- rnorm(Nsim, mu, sigma/y)

d1 <- cumsum(exp(-x^2))/(1 : Nsim)
d2 <- cumsum(exp(-mu^2/(1 + 2 * (sigma/y)^2))/sqrt(1 + 2 * (sigma/y)^2))/(1 : Nsim)

plot(d1, type = "l")
lines(d2, col = "goldenrod")