library(diversitree)

tr <- read.tree('carex-all-0.001-dated-966-ultrametric-2017-08-02.tre')

states <- read.delim("carex-two-regions-broad.txt", row.names = 1, as.is = TRUE)
tip.states <- states[tr$tip.label, 'state'] #A=Nonhimalaya (1,if present only in nonhimalaya), B=himalaya (2,if present in himalaya), AB=both (0,if present in both)
names(tip.states) <- tr$tip.label
model.full.4c <- make.geosse(tr,  tip.states, sampling.f = c(128/239, (966-128-51)/(1912-239-88), 51/88))
							 # sampled 128 out of estimated 239 himalaya spp in which 88 are endemic in himalaya, of which here 51 species are included, and 841 non-himalayan specie out of 1912 
model.4c.no.sAB <- constrain(model.full.4c, sAB ~ 0)
model.4c.noDiff <- constrain(model.full.4c, sA ~ sB, xA ~ xB)

model.4c.p <- starting.point.geosse(tr, eps=0.5)
model.4c.fitted.full <- find.mle(model.full.4c, model.4c.p)
model.4c.fitted.no.sAB <- find.mle(model.4c.no.sAB, model.4c.p)
model.4c.fitted.noDiff <- find.mle(model.4c.noDiff, model.4c.p)

geosse.4c.results <- round(rbind(full = c(coef(model.4c.fitted.full), lnL = logLik(model.4c.fitted.full), aic = AIC(model.4c.fitted.full)),
            no.sAB = c(coef(model.4c.fitted.no.sAB, T), lnL = logLik(model.4c.fitted.no.sAB), aic = AIC(model.4c.fitted.no.sAB)),
            noDiff = c(coef(model.4c.fitted.noDiff, T), lnL = logLik(model.4c.fitted.noDiff), aic = AIC(model.4c.fitted.noDiff))), 6)
write.csv(geosse.4c.results, 'Himalaya_broad_geosse.4c.MLfit.966.csv')

geosse.4c.anova = anova(model.4c.fitted.full, no.sAB = model.4c.fitted.no.sAB, noDiff = model.4c.fitted.noDiff)

## MCMC, GEOSSE, NO sAB
p.4c.mcmc.no.sAB <- coef(model.4c.fitted.no.sAB)
prior <- make.prior.exponential(1/2)

set.seed(1)
tmp <- mcmc(model.4c.no.sAB, p.4c.mcmc.no.sAB, nsteps=100, prior=prior, w=1, print.every=0)
w <- diff(sapply(tmp[2:7], quantile, c(0.025, 0.975)))
rm(tmp)

mcmc.4c.no.sAB <- mcmc(model.4c.no.sAB, p.4c.mcmc.no.sAB, nsteps = 100000, prior = prior, w = w)

mcmc.4c.diff <- with(mcmc.4c.no.sAB, data.frame(s.diff=sA-sB,
  x.diff=xA-xB, d.diff=dA-dB, div.A=sA-xA, div.B=sB-xB, div.diff=(sA-xA)-(sB-xB)))

colMeans(mcmc.4c.diff > 0)

pdf('mcmc.4c.plot.no.sAB.longRun.broad.966.div.pdf')
col1 <- c("red", "green", "blue", "purple", "black", "brown", "orange", "grey", "pink")
col2 <- col1[c(1, 3, 5, 7, 8, 9 )]
mcmc.4c.diff <- with(mcmc.4c.no.sAB, data.frame(div.A=sA-xA, div.B=sB-xB, s.diff=sA-sB, x.diff=xA-xB, d.diff=dA-dB, div.diff=(sA-xA)-(sB-xB)))

par(mfrow=c(2,1), mar=c(3, 4, 0, 1))
profiles.plot(mcmc.4c.no.sAB[2:7], col.line=col1, n.br = 100, xlab="", ylab="")
legend("top", "right", argnames(model.4c.no.sAB), col=col1, lty=1, lwd=1)
profiles.plot(mcmc.4c.diff, col.line=col2, n.br = 100, xlab="", ylab="")
legend("top", "right", colnames(mcmc.4c.diff), col=col2, lty=1, lwd=1)
title(xlab="Rate", ylab="Posterior probability density", outer=T, line=-1)
dev.off()