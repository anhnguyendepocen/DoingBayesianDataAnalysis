
library(snow)
library(rsprng)
cl = makeCluster(2, type = "MPI")
clusterSetupSPRNG(cl)
parApply(cl, X = marks, MARGIN = 2, FUN = mean)
mean(marks)

clusterExport(cl, list("marks"))
clusterEvalQ(cl, ls())
parLapply(cl, c("pearson", "spearman"), function(m) { cor(marks, method = m) } )
slave.id = function(id) { 
  assign("id", value = id, envir = .GlobalEnv)
}
parSapply(cl, 1:2, slave.id)
calls = list(mean, sd)
clusterExport(cl, list("calls"))
clusterEvalQ(cl, calls[[id]](marks))
stopCluster(cl)

#-----------------------------------------------------------------------------#

data(hailfinder)
cl = makeCluster(2, type = "MPI")
unlist(clusterEvalQ(cl, .test.counter))
.test.counter
cl = makeCluster(3, type = "MPI")
res = gs(hailfinder, cluster = cl)
unlist(clusterEvalQ(cl, .test.counter))
stopCluster(cl)
cl = makeCluster(4, type = "MPI")
res = gs(hailfinder, cluster = cl)
unlist(clusterEvalQ(cl, .test.counter))
stopCluster(cl)
system.time(gs(hailfinder))
ntests(gs(hailfinder, optimized = TRUE))
ntests(gs(hailfinder, optimized = FALSE))

#-----------------------------------------------------------------------------#

cl = makeCluster(4, type = "MPI")
clusterEvalQ(cl, library(bnlearn))
start = random.graph(names(hailfinder), num = 4, method = "melancon")
parallel.multistart = function(net) {
  hc(hailfinder, start = net)
}
netlist = parLapply(cl, start, parallel.multistart)
unlist(lapply(netlist, score, data = hailfinder))
parallel.multistart = function(net) {
  tabu(hailfinder, start = net)
}
netlist = parLapply(cl, start, parallel.multistart)
unlist(lapply(netlist, score, data = hailfinder))
s0 = random.graph(names(hailfinder), method = "melancon")
system.time(tabu(hailfinder, start = s0))
system.time(parLapply(cl, start, parallel.multistart))

#-----------------------------------------------------------------------------#

sparse = bn.boot(hailfinder, algorithm = "hc", R = 200, statistic = narcs)
summary(unlist(sparse))
quantile(unlist(sparse), c(0.05, 0.95))
system.time(bn.boot(hailfinder, algorithm = "hc", R = 200, statistic = narcs))
cl = makeCluster(2, type = "MPI")
system.time(bn.boot(hailfinder, algorithm = "hc", R = 200, statistic = narcs, cluster = cl))
stopCluster(cl)

#-----------------------------------------------------------------------------#

bn.cv(hailfinder, 'mmhc', loss = "pred", loss.args = list(target = "CompPlFcst"))
system.time(bn.cv(hailfinder, 'mmhc', loss = "pred", loss.args = list(target = "CompPlFcst")))
cl = makeCluster(2, type = "MPI")
system.time(bn.cv(hailfinder, 'mmhc', loss = "pred", loss.args = list(target = "CompPlFcst"), cluster = cl))
stopCluster(cl)
cl = makeCluster(3, type = "MPI")
system.time(bn.cv(hailfinder, 'mmhc', loss = "pred", loss.args = list(target = "CompPlFcst"), cluster = cl))
stopCluster(cl)
cl = makeCluster(4, type = "MPI")
system.time(bn.cv(hailfinder, 'mmhc', loss = "pred", loss.args = list(target = "CompPlFcst"), cluster = cl))
stopCluster(cl)
naive = naive.bayes(training = "CompPlFcst", data = hailfinder)
bn.cv(hailfinder, naive, loss = "pred")

#-----------------------------------------------------------------------------#

fitted = bn.fit(mmhc(hailfinder), hailfinder)
n = nrow(hailfinder)
summary(hailfinder[, "CompPlFcst"]) / n
cp = cpdist(fitted, nodes = "CompPlFcst", (InsInMt == "Strong") & (CldShadeConv == "Marked"), n = 10^7)
n = nrow(cp)
summary(cp[, CompPlFcst]) / n
system.time(cpdist(fitted, nodes = "CompPlFcst", (InsInMt == "Strong") & (CldShadeConv == "Marked"), n = 10^7, batch = 10^6))
cl = makeCluster(2, "SOCK")
system.time(cpdist(fitted, nodes = "CompPlFcst", (InsInMt == "Strong") & (CldShadeConv == "Marked"), n = 10^7, cluster = cl))
stopCluster(cl)
cpquery(fitted, (WindFieldMt == "Westerly"), (WindFieldPln == "E_NE"), n = 10^7)
n = nrow(hailfinder)
summary(hailfinder[, "WindFieldMt"]) / n
system.time(cpquery(fitted, (WindFieldMt == "Westerly"), (WindFieldPln == "E_NE"), n = 10^7, batch = 10^6))
cl = makeCluster(2, type = "MPI")
system.time(cpquery(fitted, (WindFieldMt == "Westerly"), (WindFieldPln == "E_NE"), n = 10^7, cluster = cl))
stopCluster(cl)


