library(diversitree)

tr <- read.tree('carex-all-0.001-dated-966-ultrametric-2017-08-02.tre')

states <- read.delim("carex-two-regions-narrow.txt", row.names = 1, as.is = TRUE)
tip.states <- states[tr$tip.label, 'state'] #A=Nonhimalaya (1,if present only in nonhimalaya), B=himalaya (2,if present in himalaya), AB=both (0,if present in both)
names(tip.states) <- tr$tip.label
model.full.4c <- make.geosse(tr,  tip.states, sampling.f = c(105/189, (966-105-10)/(1960-189-40), 10/40))
							 # sampled 105 out of estimated 189 himalaya spp in which 40 are endemic in himalaya, of which here 10 species are included 
model.4c.no.sAB <- constrain(model.full.4c, sAB ~ 0)
model.4c.noDiff <- constrain(model.full.4c, sA ~ sB, xA ~ xB)

model.4c.p <- starting.point.geosse(tr, eps=0.5)
model.4c.fitted.full <- find.mle(model.full.4c, model.4c.p)
model.4c.fitted.no.sAB <- find.mle(model.4c.no.sAB, model.4c.p)
model.4c.fitted.noDiff <- find.mle(model.4c.noDiff, model.4c.p)

geosse.4c.results <- round(rbind(full = c(coef(model.4c.fitted.full), lnL = logLik(model.4c.fitted.full), aic = AIC(model.4c.fitted.full)),
            no.sAB = c(coef(model.4c.fitted.no.sAB, T), lnL = logLik(model.4c.fitted.no.sAB), aic = AIC(model.4c.fitted.no.sAB)),
            noDiff = c(coef(model.4c.fitted.noDiff, T), lnL = logLik(model.4c.fitted.noDiff), aic = AIC(model.4c.fitted.noDiff))), 6)
write.csv(geosse.4c.results, 'Himalaya_narrow_966_geosse.4c.MLfit.csv')

geosse.4c.anova = anova(model.4c.fitted.full, no.sAB = model.4c.fitted.no.sAB, noDiff = model.4c.fitted.noDiff)

dev.off()