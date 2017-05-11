
library(vars)
data(Canada)
VAR(Canada, p = 2)
summary(VAR(Canada, p = 2))
VAR(Canada, p = 2, type = "none")
VAR(Canada, p = 2, type = "trend")
VAR(Canada, p = 2, type = "both")
VAR(Canada, lag.max = 4, ic = "AIC") 
VAR(Canada, lag.max = 4, ic = "SC") 
var.2c = VAR(Canada, p = 2, type = "const")
stab = stability(var.2c, type = "OLS-CUSUM")
plot(stab)
normality.test(var.2c)
normality.test(var.2c, multivariate.only = FALSE)
serial.test(var.2c, lags.pt = 16, type = "PT.adjusted")
arch.test(var.2c)
arch.test(var.2c, multivariate.only = FALSE)

#-----------------------------------------------------------------------------#

library(lars)
library(GeneNet)
data(arth800)
subset = c(60, 141, 260, 333, 365, 424, 441, 512, 521, 578, 789, 799)
arth12 = arth800.expr[, subset]
x = arth12[1:(nrow(arth12) - 2), ]     
y = arth12[-(1:2), "265768_at"]       
lasso.fit = lars(y = y, x = x, type = "lasso") 
fit.all = lapply(colnames(arth12),
   function(gene) { 
     y = arth12[-(1:2), gene]
     lars(y = y, x = x, type = "lasso") 
   })
plot(lasso.fit)
coef(lasso.fit)
lasso.cv = cv.lars(y = y, x = x, mode = "fraction") 
frac = lasso.cv$index[which.min(lasso.cv$cv)]
predict(lasso.fit, s = frac, type = "coef", mode = "fraction")
predict(lasso.fit, s = 3, type = "coef", mode = "step")$coefficients
predict(lasso.fit, s = 0.2, type = "coef", mode = "lambda")$coefficients
lar.fit = lars(y = y, x = x, type = "lar") 
plot(lar.fit)
lar.cv = cv.lars(y = y, x = x, type = "lar")
step.fit = lars(y = y, x = x, type = "stepwise") 
plot(step.fit)
step.cv = cv.lars(y = y, x = x, type = "stepwise")     

#-----------------------------------------------------------------------------#

library(simone)
simone(arth12, type = "time-course") 
ctrl = setOptions(clusters.crit = "BIC")
simone(arth12, type = "time-course", clustering = TRUE, control = ctrl) 
plot(simone(arth12, type = "time-course",
       clustering = TRUE, control = ctrl), output = "BIC")
plot(simone(arth12, type = "time-course", 
       clustering = TRUE, control = ctrl))

#-----------------------------------------------------------------------------#

library(GeneNet)
dyn = ggm.estimate.pcor(arth800.expr, method = "dynamic")
arth.arcs = network.test.edges(dyn)
arth.net = extract.network(arth.arcs, method.ggm = "number", cutoff.ggm = 10)
arth.net = extract.network(arth.arcs, method.ggm = "prob", cutoff.ggm = 0.05)

#-----------------------------------------------------------------------------#

library(G1DBN)
data(arth800line)
subset = c(60, 141, 260, 333, 365, 424, 441, 512, 521, 578, 789, 799)
arth12 = as.matrix(arth800line[, subset])
step1 = DBNScoreStep1(arth12, method = "ls")
round(step1$S1ls, 2)[1:6, 1:6]
edgesG1 = BuildEdges(score = step1$S1ls, threshold = 0.50,  prec = 6)
nrow(edgesG1)
step2 = DBNScoreStep2(step1$S1ls, data = arth12,
       method = "ls", alpha1 = 0.50)
edgesG = BuildEdges(score = step2, threshold = 0.05, prec = 6)

#-----------------------------------------------------------------------------#

library(ARTIVA)
data(simulatedProfiles)
targets = c("1", "10", "20", "TF3", "45", "50")
parents = c("TF1", "TF2", "TF3", "TF4", "TF5")
DBN = ARTIVAnet(
        targetData = simulatedProfiles[targets, ],
        parentData = simulatedProfiles[parents, ], 
        targetNames = targets,
        parentNames = parents,
        niter = 50000,
        savePictures = FALSE)
head(ARTIVAtest1[, -7])

