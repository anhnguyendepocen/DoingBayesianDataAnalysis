
isachs = read.table("sachs.interventional.txt",
           header = TRUE, colClasses = "factor")
library(gRain)
library(bnlearn)
val.str = paste("[PKC][PKA|PKC][praf|PKC:PKA]",
           "[pmek|PKC:PKA:praf][p44.42|pmek:PKA]",
           "[pakts473|p44.42:PKA][P38|PKC:PKA]",
           "[pjnk|PKC:PKA][plcg][PIP3|plcg]",
           "[PIP2|plcg:PIP3]")
val = model2network(val.str)
isachs = isachs[, 1:11]
for (i in names(isachs))
  levels(isachs[, i]) = c("LOW", "AVG", "HIGH")
fitted = bn.fit(val, isachs, method = "bayes")
jtree = compile(as.grain(fitted))
jprop = setFinding(jtree, nodes = "p44.42",
          states = "LOW")
querygrain(jtree, nodes = "pakts473")$pakts473
querygrain(jprop, nodes = "pakts473")$pakts473
querygrain(jtree, nodes = "PKA")$PKA
querygrain(jprop, nodes = "PKA")$PKA
names(which.max(querygrain(jprop,
                    nodes = c("PKA"))$PKA))
particles = cpdist(fitted, nodes = "pakts473", 
              evidence = (p44.42 == "LOW"))
prop.table(table(particles))
particles = cpdist(fitted, nodes = "PKA",
              evidence = (p44.42 == "LOW"))
prop.table(table(particles))
cpquery(fitted,
  event = (pakts473 == "LOW") & (PKA != "HIGH"),
  evidence = (p44.42 == "LOW") | (praf == "LOW"))

#-----------------------------------------------------------------------------#

x = arth12[1:(nrow(arth12) - 2), ]
y = arth12[-(1:2), "265768_at"]
lasso.fit = lars(y = y, x = x, type = "lasso")
lasso.cv = cv.lars(y = y, x = x, mode = "fraction")
frac = lasso.cv$index[which.min(lasso.cv$cv)]
lasso.est = predict(lasso.fit, type = "fit",
              newx = x, s = frac, 
              mode = "fraction")$fit
lasso.est
lasso.pred = predict(lasso.fit, type = "fit",
               newx = arth12[c("24-1", "24-2"), ],
               s = frac, mode = "fraction")$fit
lasso.pred
library(penalized)
lambda = optL1(response = y, penalized = x)$lambda
lasso.t = penalized(response = y, penalized = x,
            lambda1 = lambda)
coef(lasso.t)
dbn1 = 
  model2network("[245094_at][265768_at|245094_at]")
xp.mean = mean(x[, "245094_at"])
xp.sd = sd(x[, "245094_at"])
dbn1.fit = 
  custom.fit(dbn1, 
    dist = list("245094_at" = list(coef = xp.mean, 
             sd = xp.sd), "265768_at" = lasso.t))
cpquery(dbn1.fit, event = (`265768_at` > 8),
                  evidence = (`245094_at` > 8))
cpquery(dbn1.fit, event = (`265768_at` > 8),
                  evidence = (`245094_at` < 8))
dist.low = cpdist(dbn1.fit, node = "265768_at", 
             evidence = (`245094_at` < 8))
dist.high = cpdist(dbn1.fit, node = "265768_at",
              evidence = (`245094_at` > 8))
y = arth12[-(1:2), "245094_at"]
colnames(x)[12] = "245094_at1"
lambda = optL1(response = y, penalized = x)$lambda
lasso.s = penalized(response = y, penalized = x,
             lambda1 = lambda)
coef(lasso.s)
dbn2 = empty.graph(c("265768_at", "245094_at", 
        "258736_at", "257710_at", "255070_at",
        "245319_at", "245094_at1"))
dbn2 = set.arc(dbn2, "245094_at", "265768_at")
for (node in names(coef(lasso.s))[-c(1, 6)])
  dbn2 = set.arc(dbn2, node, "245094_at")
dbn2 = set.arc(dbn2, "245094_at1", "245094_at")
dbn2.data = as.data.frame(x[, nodes(dbn2)[1:6]])
dbn2.data[, "245094_at"] = y
dbn2.data[, "245094_at1"] = x[, "245094_at"]
dbn2.fit = bn.fit(dbn2, dbn2.data)
dbn2.fit[["265768_at"]] = lasso.t 
dbn2.fit[["245094_at"]] = lasso.s
cpquery(dbn2.fit, event = (`265768_at` > 8), 
  evidence = (`245094_at` > 8) & (`245094_at1` > 8))
cpquery(dbn2.fit, event = (`265768_at` > 8),
 evidence = (`245094_at1` > 7) & (`245094_at1` < 8))
cpquery(dbn2.fit, event = (`265768_at` > 8),
  evidence = TRUE)
cpd = 
  cpdist(dbn2.fit, node = "245094_at", evidence = 
    (`245094_at1` > 6.5) & (`245094_at1` < 7.5) & 
    (`265768_at` > 7) & (`265768_at` < 8))
summary(cpd)

