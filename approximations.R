# set wd, load and clean data ---------------------------------------

setwd("/Users/Bryan/Dropbox/frogs/data/Moen et al. 2013/")
library(ape); library(geiger); library(RRPP); library(ggplot2); library(geomorph); library(reshape2)

dat <- read.csv("MoenFrogs.csv")
dat <- dat[, c(1, 2, 7:8, 10, 13:17, 29:32, 33:34, 36:38, 53, 56)]
colnames(dat)[c(20:21)] <- c("v", "P.Mass")
dat <- dat[-c(93, 195, 207),] #remove individual with no performance
rem <- which(dat$Species == "Allobates_femoralis" | dat$Species == "Calluella_yunnanensis" | dat$Species == "Dendropsophus_rhodopeplus" | dat$Species == "Litoria_inermis" | dat$Species == "Nanorana_yunnanensis" | dat$Species == "Sphaenorhynchus_lacteus") #species with only 1 observation for either morphology, performance, or both
dat <- dat[-rem, ]
attach(dat)

# Mass regression ---------------------------------------------------------

anova(lm(log(MASS_specimen) ~ log(MASS_live)))
summary(lm(log(MASS_specimen) ~ log(MASS_live)))
confint(lm(log(MASS_specimen) ~ log(MASS_live)))

# Calc. proxies, performance, clean ---------------------------------------

#Mean species performance
Mass <- MASS_live / 1000 #convert to kilograms
Svl <- SVL / 1000 #convert to m

#calculate species means
V <- tapply(v, Species, mean) #m / s
E <- tapply(0.5 * Mass * v ^ 2, Species, mean) #Joules
P <- tapply(P.Mass * Mass, Species, mean) #Watts
svl <- tapply(Svl, Species, mean) # m
mass <- tapply(Mass, Species, mean) #kg

#remove NAs
V <- V[complete.cases(V)]
E <- E[complete.cases(E)]
P <- P[complete.cases(P)]
svl <- svl[complete.cases(svl)]
mass <- mass[complete.cases(mass)]

#Size-specific performance
V <- V / svl
E <- E / mass
P <- P / mass

#Mean species proxies
musc.dat <- read.csv("/Users/Bryan/Dropbox/frogs/Data/interspecific_proxies/Intraspecific_dataset_Mendoza et al.csv")
rem.spp <- setdiff(musc.dat$Species, dat$Species) #length 30
for(i in 1:30){
 musc.dat <- musc.dat[-which(musc.dat$Species == rem.spp[i]), ]   
} #keep only leg masses that are in Moen et al. 2013
musc.dat <- musc.dat[match(dat$Museum_tag, musc.dat$Specimen.tag..), ]
leg.mass <- musc.dat$Live_mass..g. * musc.dat$Muscle.mass.Body.mass.ratio / 1000 #calculate live muscle mass and convert to kilograms

FPCSA <- (leg.mass) ^ (2/3) #kg
L <- (FEMUR + TIBL + TARS + FOOT) / 1000 #m
v.proxy <- sqrt(2 * L * FPCSA / Mass) 

L.solo  <- tapply(L, Species, mean) #m
leg.mass.solo <- tapply(leg.mass, Species, mean) #kg
V.proxy <- tapply(v.proxy, Species, mean)
E.proxy <- tapply(L * FPCSA, Species, mean)
P.proxy <- tapply(sqrt(2 * L * (FPCSA ^ 3) / Mass), Species, mean)



L.solo <- L.solo[complete.cases(L.solo)]
leg.mass.solo <- leg.mass.solo[complete.cases(leg.mass.solo)]
V.proxy <- V.proxy[complete.cases(V.proxy)]
E.proxy <- E.proxy[complete.cases(E.proxy)]
P.proxy <- P.proxy[complete.cases(P.proxy)]

V.proxy <- V.proxy / svl
E.proxy <- E.proxy / mass
P.proxy <- P.proxy / mass

#Stick into data frame and remove rows with NAs
df.V <- data.frame(V.proxy = V.proxy, V = V, L = L.solo, leg.mass = leg.mass.solo, svl = svl, mass = mass)[complete.cases(cbind(V.proxy, V, L.solo, leg.mass.solo, svl, mass)),]
df.E <- data.frame(E.proxy = E.proxy, E = E, L = L.solo, leg.mass = leg.mass.solo, svl = svl, mass = mass)[complete.cases(cbind(E.proxy, E, L.solo, leg.mass.solo, svl, mass)),]
df.P <- data.frame(P.proxy = P.proxy, P = P, L = L.solo, leg.mass = leg.mass.solo, svl = svl, mass = mass)[complete.cases(cbind(P.proxy, P, L.solo, leg.mass.solo, svl, mass)),]

#make substitutions in Moen et al. 2013 data to match names from Jetz and Pyron 2018
rownames(df.V)[5] <- "Syncope_bassleri"
rownames(df.V)[18] <- "Cyclorana_australis"
rownames(df.V)[22] <- "Cyclorana_longipes"

rownames(df.E)[5] <- "Syncope_bassleri"
rownames(df.E)[18] <- "Cyclorana_australis"
rownames(df.E)[22] <- "Cyclorana_longipes"

rownames(df.P)[5] <- "Syncope_bassleri"
rownames(df.P)[18] <- "Cyclorana_australis"
rownames(df.P)[22] <- "Cyclorana_longipes"

#
dat2 <- read.csv("/Users/Bryan/Dropbox/frogs/Data/interspecific_proxies/Data_Astley_2016.csv")
dat2 <- dat2[complete.cases(dat2), ]
dat2 <- dat2[-c(15:17), ]#remove Heterixalus
spp <- dat2$Species

Mass2 <- dat2$mass_g / 1000 #convert to kilograms
musc.mass <- dat2$muscle_mass_percent * Mass2
V2 <- tapply(dat2$V_svl_s * dat2$SVL_mm / 1000, spp, mean)
E2 <- tapply(dat2$E_J_kg_muscle * musc.mass, spp, mean)
P2 <- tapply(dat2$P_W_kg_muscle * musc.mass, spp, mean)

svl2 <- tapply(dat2$SVL_mm / 1000, spp, mean)
mass2 <- tapply(Mass2, spp, mean)

#remove NAs
V2 <- V2[complete.cases(V2)]
E2 <- E2[complete.cases(E2)]
P2 <- P2[complete.cases(P2)]
svl2 <- svl2[complete.cases(svl2)]
mass2 <- mass2[complete.cases(mass2)]

#Size-specific performance
V2 <- V2 / svl2
E2 <- E2 / mass2
P2 <- P2 / mass2

#Mean species proxies
L2 <- dat2$leg_svl * dat2$SVL_mm / 1000
CSA2 <- musc.mass ^ (2/3)

L2.solo <- tapply(L2, spp, mean)
leg.mass2.solo <- tapply(musc.mass, spp, mean)
V2.proxy <- tapply(sqrt(2 * L2 * CSA2 / Mass2), spp, mean)
E2.proxy <- tapply(L2 * CSA2, spp, mean)
P2.proxy <- tapply(sqrt(2 * L2 * (CSA2 ^ 3) / Mass2), spp, mean)

L2.solo <- L2.solo[complete.cases(L2.solo)]
leg.mass2.solo <- leg.mass2.solo[complete.cases(leg.mass2.solo)]
V2.proxy <- V2.proxy[complete.cases(V2.proxy)]
E2.proxy <- E2.proxy[complete.cases(E2.proxy)]
P2.proxy <- P2.proxy[complete.cases(P2.proxy)]

V2.proxy <- V2.proxy / svl2
E2.proxy <- E2.proxy / mass2
P2.proxy <- P2.proxy / mass2

#Stick into data frame and remove rows with NAs
df.V2 <- data.frame(V.proxy = V2.proxy, V = V2, L = L2.solo, leg.mass = leg.mass2.solo, svl = svl2, mass = mass2)[complete.cases(cbind(V2.proxy, V2, L2.solo, leg.mass2.solo, svl2, mass2)),]
df.E2 <- data.frame(E.proxy = E2.proxy, E = E2, L = L2.solo, leg.mass = leg.mass2.solo, svl = svl2, mass = mass2)[complete.cases(cbind(E2.proxy, E2, L2.solo, leg.mass2.solo, svl2, mass2)),]
df.P2 <- data.frame(P.proxy = P2.proxy, P = P2, L = L2.solo, leg.mass = leg.mass2.solo, svl = svl2, mass = mass2)[complete.cases(cbind(P2.proxy, P2, L2.solo, leg.mass2.solo, svl2, mass2)),]

#rename Astley data to match tree
rownames(df.V2)[5] <- "Rana_pipiens"
rownames(df.E2)[5] <- "Rana_pipiens"
rownames(df.P2)[5] <- "Rana_pipiens"

# PGLS ---------------------------------------------------

tree <- read.nexus("/Users/Bryan/Dropbox/frogs/Data/interspecific_proxies/post.tree.branches.nex")

# Calculate independent contrasts, prune tree and data
treedat.V <- treedata(tree, as.matrix(rbind(df.V, df.V2)))
phy <- treedat.V$phy
df1.V <- rrpp.data.frame(V = treedat.V$data[, 2], V.proxy = treedat.V$data[, 1])
fit.V <- lm.rrpp(V.proxy ~ V, data = df1.V, Cov = vcv.phylo(phy), iter = 9999)
anova(fit.V)

treedat.E <- treedata(tree, as.matrix(rbind(df.E, df.E2)))
df1.E <- rrpp.data.frame(E = treedat.E$data[, 2], E.proxy = treedat.E$data[, 1])
fit.E <- lm.rrpp(E.proxy ~ E, data = df1.E, Cov = vcv.phylo(phy), iter = 9999)
anova(fit.E)

treedat.P <- treedata(tree, as.matrix(rbind(df.P, df.P2)))
df1.P <- rrpp.data.frame(P = treedat.P$data[, 2], P.proxy = treedat.P$data[, 1])
fit.P <- lm.rrpp(P.proxy ~ P, data = df1.P, Cov = vcv.phylo(phy), iter = 9999)
anova(fit.P)

    # obtain r^2 from regressions and round to 2 decimal places
rsqs.pgls <- round( c(anova(fit.V)$table$Rsq[1],
                      anova(fit.E)$table$Rsq[1],
                      anova(fit.P)$table$Rsq[1]), 2)
rsqs.pgls

#Analysis of Moen et al. data, N = 38

treedat.V2 <- treedata(tree, as.matrix(rbind(df.V)))
phy <- treedat.V2$phy
df2.V <- rrpp.data.frame(V = treedat.V2$data[, 2], V.proxy = treedat.V2$data[, 1])
fit2.V <- lm.rrpp(V.proxy ~ V, data = df2.V, Cov = vcv.phylo(phy), iter = 9999)
anova(fit2.V)

treedat.E2 <- treedata(tree, as.matrix(rbind(df.E)))
df2.E <- rrpp.data.frame(E = treedat.E2$data[, 2], E.proxy = treedat.E2$data[, 1])
fit2.E <- lm.rrpp(E.proxy ~ E, data = df2.E, Cov = vcv.phylo(phy), iter = 9999)
anova(fit2.E)

treedat.P2 <- treedata(tree, as.matrix(rbind(df.P)))
df2.P <- rrpp.data.frame(P = treedat.P2$data[, 2], P.proxy = treedat.P2$data[, 1])
fit2.P <- lm.rrpp(P.proxy ~ P, data = df2.P, Cov = vcv.phylo(phy), iter = 9999)
anova(fit2.P)

rsqs.pgls2 <- round( c(anova(fit2.V)$table$Rsq[1],
                       anova(fit2.E)$table$Rsq[1],
                       anova(fit2.P)$table$Rsq[1]), 2) 
rsqs.pgls2 # obtain r^2 from regressions and round to 2 decimal places

# Analysis of individual size-standardized variables
phy <- treedat.V$phy # Get tree
df3.V <- rrpp.data.frame(V =        treedat.V$data[, 2],
                         V.proxy =  treedat.V$data[, 1], 
                         L =        treedat.V$data[, 3] / treedat.V$data[, 5], 
                         leg.mass = treedat.V$data[, 4] / treedat.V$data[, 5], 
                         mass =     treedat.V$data[, 6] / treedat.V$data[, 5]) 

fit3.V <- lm.rrpp(V ~ V.proxy,  Cov = vcv.phylo(phy), data = df3.V, iter = 9999) 
fit4.V <- lm.rrpp(V ~ L,        Cov = vcv.phylo(phy), data = df3.V, iter = 9999) 
fit5.V <- lm.rrpp(V ~ leg.mass, Cov = vcv.phylo(phy), data = df3.V, iter = 9999) 
fit6.V <- lm.rrpp(V ~ mass,     Cov = vcv.phylo(phy), data = df3.V, iter = 9999) 

model.comparison(fit3.V, 
                 fit4.V, 
                 fit5.V, 
                 fit6.V, type = "logLik", tol = 0.01) 
    # V.proxy is prefered, deltaAIC is 56.

    # obtain r^2 from regressions and round to 2 decimal places
round(c(anova(fit3.V)$table[1, 4], 
        anova(fit4.V)$table[1, 4], 
        anova(fit5.V)$table[1, 4], 
        anova(fit6.V)$table[1, 4]), 2) 

df3.E <- rrpp.data.frame(E =        treedat.E$data[, 2], 
                         E.proxy =  treedat.E$data[, 1], 
                         L =        treedat.E$data[, 3] / treedat.E$data[, 6], 
                         leg.mass = treedat.E$data[, 4] / treedat.E$data[, 6])

fit3.E <- lm.rrpp(E ~ E.proxy,  Cov = vcv.phylo(phy), data = df3.E, iter = 9999) 
fit4.E <- lm.rrpp(E ~ L,        Cov = vcv.phylo(phy), data = df3.E, iter = 9999) 
fit5.E <- lm.rrpp(E ~ leg.mass, Cov = vcv.phylo(phy), data = df3.E, iter = 9999) 

model.comparison(fit3.E, 
                 fit4.E, 
                 fit5.E, type = "logLik", tol = 0.01)

    # obtain r^2 from regressions and round to 2 decimal places
round(c(anova(fit3.E)$table[1, 4], 
        anova(fit4.E)$table[1, 4], 
        anova(fit5.E)$table[1, 4]), 2) 


df3.P <- rrpp.data.frame(P =        treedat.P$data[, 2], 
                         P.proxy =  treedat.P$data[, 1], 
                         L =        treedat.P$data[, 3] / treedat.P$data[, 6], 
                         leg.mass = treedat.P$data[, 4] / treedat.P$data[, 6])

fit3.P <- lm.rrpp(P ~ P.proxy,  Cov = vcv.phylo(phy), data = df3.P, iter = 9999) 
fit4.P <- lm.rrpp(P ~ L,        Cov = vcv.phylo(phy), data = df3.P, iter = 9999) 
fit5.P <- lm.rrpp(P ~ leg.mass, Cov = vcv.phylo(phy), data = df3.P, iter = 9999) 

model.comparison(fit3.P, 
                 fit4.P, 
                 fit5.P, type = "logLik", tol = 0.01)

    # obtain r^2 from regressions and round to 2 decimal places
round(c(anova(fit3.P)$table[1, 4], 
        anova(fit4.P)$table[1, 4], 
        anova(fit5.P)$table[1, 4]), 2)

# PGLS Plots -------------------------------------------------------------------

setwd("/Users/Bryan/Dropbox/frogs/Figures/interspecific_proxies")
pdf("Fig_1_velocity.pdf", width = 4, height = 4)

upper <- predict(fit.V)$ucl
lower <- predict(fit.V)$lcl
plot.df.V <- cbind(data.frame(V = treedat.V$data[, 2], V.proxy = treedat.V$data[, 1]), upper, lower)

ggplot(data = plot.df.V, aes(x = V, y = V.proxy)) +
    geom_point(size = 2) +
    geom_abline(intercept = fit.V$LM$gls.coefficients[[1]], slope = fit.V$LM$gls.coefficients[[2]], color = "black", size = 1.1) +
    geom_line(aes(y = lower), color = "black", linetype = "dashed") +
    geom_line(aes(y = upper), color = "black", linetype = "dashed") +
    theme(axis.text.x = element_text(size = 8, color = "black")) +
    theme(axis.text.y = element_text(size = 8, color = "black")) +
    theme(axis.title.x = element_text(size = 13)) +
    theme(axis.title.y = element_text(size = 13)) +
    xlab("Observed Velocity (Body Lengths/s)") + ylab("Velocity Appx. (Body Lengths/s)") +

    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"))

dev.off()

pdf("Fig_2_energy.pdf", width = 4, height = 4)

upper <- predict(fit.E)$ucl
lower <- predict(fit.E)$lcl
plot.df.E <- cbind(data.frame(E = treedat.E$data[, 2], E.proxy = treedat.E$data[, 1]), upper, lower)

ggplot(data = plot.df.E, aes(x = E, y = E.proxy)) +
    geom_point(size = 2) +
    geom_abline(intercept = fit.E$LM$gls.coefficients[[1]], slope = fit.E$LM$gls.coefficients[[2]], color = "black", size = 1.1) +
    geom_line(aes(y = lower), color = "black", linetype = "dashed") +
    geom_line(aes(y = upper), color = "black", linetype = "dashed") +
    theme(axis.text.x = element_text(size = 8, color = "black")) +
    theme(axis.text.y = element_text(size = 8, color = "black")) +
    theme(axis.title.x = element_text(size = 13)) +
    theme(axis.title.y = element_text(size = 13)) +
    xlab("Observed Mass-Specific Energy (J/kg)") + ylab("Mass-Specific Energy Appx. (J/kg)") +

    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"))

dev.off()

pdf("Fig_3_power.pdf", width = 4, height = 4)

upper <- predict(fit.P)$ucl
lower <- predict(fit.P)$lcl
plot.df.P <- cbind(data.frame(P = treedat.P$data[, 2], P.proxy = treedat.P$data[, 1]), upper, lower)

ggplot(data = plot.df.P, aes(x = P, y = P.proxy)) +
    geom_point(size = 2) +
    geom_abline(intercept = fit.P$LM$gls.coefficients[[1]], slope = fit.P$LM$gls.coefficients[[2]], color = "black", size = 1.1) +
    geom_line(aes(y = lower), color = "black", linetype = "dashed") +
    geom_line(aes(y = upper), color = "black", linetype = "dashed") +
    theme(axis.text.x = element_text(size = 8, color = "black")) +
    theme(axis.text.y = element_text(size = 8, color = "black")) +
    theme(axis.title.x = element_text(size = 13)) +
    theme(axis.title.y = element_text(size = 13)) +
    xlab("Observed Mass-Specific Power (W/kg)") + ylab("Mass-Specific Power Appx. (W/kg)") +

    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"))

dev.off()

# OLS ----------------------------------------------------------------

#Regressions

fit.V <- lm.rrpp(V.proxy ~ V, data = df1.V, iter = 9999)
anova(fit.V)

fit.E <- lm.rrpp(E.proxy ~ E, data = df1.E, iter = 9999)
anova(fit.E)

fit.P <- lm.rrpp(P.proxy ~ P, data = df1.P, iter = 9999)
anova(fit.P)

anova(fit.V)$table$Rsq[1]
anova(fit.E)$table$Rsq[1]
anova(fit.P)$table$Rsq[1]

# OLS Plots -------------------------------------------------------------------

setwd("/Users/Bryan/Dropbox/frogs/Figures/interspecific_proxies")
pdf("Fig_1_velocity.pdf", width = 3, height = 3)

upper <- predict(fit.V)$ucl
lower <- predict(fit.V)$lcl
plot.df.V <- cbind(data.frame(V = treedat.V$data[, 2], V.proxy = treedat.V$data[, 1]), upper, lower)

ggplot(data = plot.df.V, aes(x = V, y = V.proxy)) +
    geom_point(size = 2) +
    geom_smooth(method = "lm", se = F, cex = 0.5, col = "black") +
    geom_line(aes(y = lower), color = "black", linetype = "dashed") +
    geom_line(aes(y = upper), color = "black", linetype = "dashed") +
    theme(axis.text.x = element_text(size = 6, color = "black")) +
    theme(axis.text.y = element_text(size = 6, color = "black")) +
    theme(axis.title.x = element_text(size = 9)) +
    theme(axis.title.y = element_text(size = 9)) +
    xlab("Velocity (m/s) / SVL (mm)") + ylab("Velocity Proxy / SVL (mm)") +

    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"))

dev.off()

pdf("Fig_2_energy.pdf", width = 3, height = 3)

upper <- predict(fit.E)$ucl
lower <- predict(fit.E)$lcl
plot.df.E <- cbind(data.frame(E = treedat.E$data[, 2], E.proxy = treedat.E$data[, 1]), upper, lower)

ggplot(data = plot.df.E, aes(x = E, y = E.proxy)) +
    geom_point(size = 2) +
    geom_smooth(method = "lm", se = F, cex = 0.5, col = "black") +
    geom_line(aes(y = lower), color = "black", linetype = "dashed") +
    geom_line(aes(y = upper), color = "black", linetype = "dashed") +
    theme(axis.text.x = element_text(size = 6, color = "black")) +
    theme(axis.text.y = element_text(size = 6, color = "black")) +
    theme(axis.title.x = element_text(size = 9)) +
    theme(axis.title.y = element_text(size = 9)) +
    xlab("Mass-Specific Energy (J/kg)") + ylab("Energy Proxy / Mass (kg)") +

    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"))

dev.off()

pdf("Fig_3_power.pdf", width = 3, height = 3)

upper <- predict(fit.P)$ucl
lower <- predict(fit.P)$lcl
plot.df.P <- cbind(data.frame(P = treedat.P$data[, 2], P.proxy = treedat.P$data[, 1]), upper, lower)

ggplot(data = plot.df.P, aes(x = P, y = P.proxy)) +
    geom_point(size = 2) +
    geom_smooth(method = "lm", se = F, cex = 0.5, col = "black") +
    geom_line(aes(y = lower), color = "black", linetype = "dashed") +
    geom_line(aes(y = upper), color = "black", linetype = "dashed") +
    theme(axis.text.x = element_text(size = 6, color = "black")) +
    theme(axis.text.y = element_text(size = 6, color = "black")) +
    theme(axis.title.x = element_text(size = 9)) +
    theme(axis.title.y = element_text(size = 9)) +
    xlab("Mass-Specific Power (W/kg)") + ylab("Power Proxy / Mass (kg)") +

    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"))

dev.off()



# Phylogenetic signal ----------------------------------------------------


appxs <- as.matrix(cbind(c(V.proxy, V2.proxy), c(E.proxy, E2.proxy), c(P.proxy, P2.proxy)))
rownames(appxs)[5] <- "Syncope_bassleri"
rownames(appxs)[18] <- "Cyclorana_australis"
rownames(appxs)[22] <- "Cyclorana_longipes"
rownames(appxs)[43] <- "Rana_pipiens"
physignal(scale(appxs), phy, iter = 9999) 

appxs <- as.matrix(cbind(c(V, V2), c(E, E2), c(P, P2)))
rownames(appxs)[5] <- "Syncope_bassleri"
rownames(appxs)[18] <- "Cyclorana_australis"
rownames(appxs)[22] <- "Cyclorana_longipes"
rownames(appxs)[43] <- "Rana_pipiens"
physignal(scale(appxs), phy, iter = 9999) 

appxs <- as.matrix(cbind(c(V.proxy, V2.proxy)))
rownames(appxs)[5] <- "Syncope_bassleri"
rownames(appxs)[18] <- "Cyclorana_australis"
rownames(appxs)[22] <- "Cyclorana_longipes"
rownames(appxs)[43] <- "Rana_pipiens"
physignal(appxs, phy, iter = 9999) 

appxs <- as.matrix(cbind(c(V, V2)))
rownames(appxs)[5] <- "Syncope_bassleri"
rownames(appxs)[18] <- "Cyclorana_australis"
rownames(appxs)[22] <- "Cyclorana_longipes"
rownames(appxs)[43] <- "Rana_pipiens"
physignal(appxs, phy, iter = 9999) 

appxs <- as.matrix(cbind(c(E.proxy, E2.proxy)))
rownames(appxs)[5] <- "Syncope_bassleri"
rownames(appxs)[18] <- "Cyclorana_australis"
rownames(appxs)[22] <- "Cyclorana_longipes"
rownames(appxs)[43] <- "Rana_pipiens"
physignal(appxs, phy, iter = 9999) 

appxs <- as.matrix(cbind(c(E, E2)))
rownames(appxs)[5] <- "Syncope_bassleri"
rownames(appxs)[18] <- "Cyclorana_australis"
rownames(appxs)[22] <- "Cyclorana_longipes"
rownames(appxs)[43] <- "Rana_pipiens"
physignal(appxs, phy, iter = 9999) 

appxs <- as.matrix(cbind(c(P.proxy, P2.proxy)))
rownames(appxs)[5] <- "Syncope_bassleri"
rownames(appxs)[18] <- "Cyclorana_australis"
rownames(appxs)[22] <- "Cyclorana_longipes"
rownames(appxs)[43] <- "Rana_pipiens"
physignal(appxs, phy, iter = 9999) 

appxs <- as.matrix(cbind(c(P, P2)))
rownames(appxs)[5] <- "Syncope_bassleri"
rownames(appxs)[18] <- "Cyclorana_australis"
rownames(appxs)[22] <- "Cyclorana_longipes"
rownames(appxs)[43] <- "Rana_pipiens"
physignal(appxs, phy, iter = 9999) 


# Phylomorphospace -------------------------------------------------------

appxs <- as.matrix(scale(cbind(c(V.proxy, V2.proxy), c(E.proxy, E2.proxy), c(P.proxy, P2.proxy))))
rownames(appxs)[5] <- "Syncope_bassleri"
rownames(appxs)[18] <- "Cyclorana_australis"
rownames(appxs)[22] <- "Cyclorana_longipes"
rownames(appxs)[43] <- "Rana_pipiens"

pca.phylo <- gm.prcomp(appxs, phy = phy)
summary(pca.phylo)
pca.phylo$rotation
pca.phylo$x


micro <- dat$microhab5[match(unique(dat$Species), dat$Species)]
names(micro) <- unique(dat$Species)
final.micro <- array(dim = c(length(labels(pca.phylo$x)[[1]]), 1), dimnames = list(labels(pca.phylo$x)[[1]]))

for(i in 1:length(names(micro))){
    final.micro[which(rownames(final.micro) == names(micro)[i])] <- as.character(micro[[i]])   
}
final.micro[45] <- "terrestrial" # Melanophryniscus stelzneri
final.micro[47] <- "semi-aquatic" #Phrynoidis aspera
final.micro[39] <- "terrestrial" #Anaxyrus fowleri
final.micro[46] <- "arboreal" #Osteopilus septentrionalis
final.micro[22] <- "burrowing" #Cyclorana longipes, Moen et al. 2013
final.micro[18] <- "burrowing" #Cyclorana australis, Moen et al. 2013
final.micro[44] <- "semi-aquatic" #Litorea aurea
final.micro[49] <- "arboreal" #Phyllomedusa hypochondrialis
final.micro[41] <- "terrestrial" #Kaloula pulchra
final.micro[5] <- "terrestrial" #Syncope bassleri, Moen et al. 2013
final.micro[48] <- "terrestrial" #Phrynomantis bifasciatus
final.micro[42] <- "terrestrial" #Kassina senegalensis
final.micro[50] <- "arboreal" #Polypedates leucomystax
final.micro[43] <- "semi-aquatic" #Rana pipiens
final.micro[51] <- "burrowing" #Scaphiopus holbrookii
final.micro[40] <- "semi-aquatic" #Bombina orientalis

final.micro <- data.frame(final.micro)
gps <- final.micro[,1]
final.micro
as.factor(final.micro)
palette(c("green", "black", "blue", "brown", "purple"))
pal <- palette(c("green", "black", "blue", "brown", "purple"))

plot.labs <- sort(labels(pca.phylo$x)[[1]])
num <- (1:51)[-c(12, 4, 25, 48, 38, 32, 51, 36, 43, 34, 27, 28, 31, 10, 11, 1, 24, 13, 50)]

for(i in 1:length(num)){
    plot.labs[num[i]] <- ""
}

gps2 <- gps
gps2[which(gps2 == "terrestrial")] <- "black"
gps2[which(gps2 == "torrential")] <- "purple"
gps2[which(gps2 == "semi-aquatic")] <- "blue"
gps2[which(gps2 == "arboreal")] <- "green"
gps2[which(gps2 == "burrowing")] <- "brown"
gps2

setwd("/Users/Bryan/Dropbox/frogs/Figures/interspecific_proxies")
pdf("Fig_4_phylomorphospace.pdf", width = 7, height = 5)

plot(pca.phylo, pch=22, phylo = TRUE, bg = gps2,
      phylo.par = list(edge.color = "darkgray", tip.labels = F, node.labels = F, anc.states = F))
 text(pca.phylo$x[match(sort(rownames(pca.phylo$x)), rownames(pca.phylo$x)),], labels = gsub("_", " ", plot.labs), pos = 3, font = 4, cex = .4, offset = .25)  #cex .4
 legend(2, 2.25, bty = "n", cex = .75, pt.cex = 1, x.intersp = 0.5, y.intersp = 1, fill = c("green",  "black", "blue", "brown", "purple"), legend = levels(as.factor(gps)))

dev.off()

# Individal Regressions ------------------------------------------------------------------

#obtain individual performance and proxy estimates
Mass <- MASS_live / 1000 #convert to kilograms
Svl <- SVL / 1000 #convert to m

#calculate species means
V <- v #m / s
E <- 0.5 * Mass * v ^ 2 #Joules
P <- P.Mass * Mass #Watts
svl <- Svl # m
mass <- Mass #kg

#remove NAs
V <- V[complete.cases(V)]
E <- E[complete.cases(E)]
P <- P[complete.cases(P)]
svl <- svl[complete.cases(svl)]
mass <- mass[complete.cases(mass)]

names(V) <- dat$Species
names(E) <- dat$Species
names(P) <- dat$Species

#Mean species proxies
musc.dat <- read.csv("/Users/Bryan/Dropbox/frogs/Data/interspecific_proxies/Intraspecific_dataset_Mendoza et al.csv")

#keep only leg masses that are in Moen et al. 2013
rem.spp <- setdiff(musc.dat$Species, dat$Species) #length 30
for(i in 1:30){
 musc.dat <- musc.dat[-which(musc.dat$Species == rem.spp[i]), ]   
}

musc.dat <- musc.dat[match(dat$Museum_tag, musc.dat$Specimen.tag..), ]

leg.mass <- musc.dat$Live_mass..g. * musc.dat$Muscle.mass.Body.mass.ratio / 1000 #calculate live muscle mass and convert to kilograms

FPCSA <- (leg.mass) ^ (2/3)
L <- (FEMUR + TIBL + TARS + FOOT) / 1000

v.proxy <- sqrt(2 * L * FPCSA / Mass)
V.proxy <- v.proxy
E.proxy <- L * FPCSA
P.proxy <- sqrt(2 * L * (FPCSA ^ 3) / Mass)

V.proxy <- V.proxy[complete.cases(V.proxy)]
E.proxy <- E.proxy[complete.cases(E.proxy)]
P.proxy <- P.proxy[complete.cases(P.proxy)]

names(V.proxy) <- dat$Species
names(E.proxy) <- dat$Species
names(P.proxy) <- dat$Species

dat2 <- read.csv("/Users/Bryan/Dropbox/frogs/Data/interspecific_proxies/Data_Astley_2016.csv")
dat2 <- dat2[complete.cases(dat2), ]
dat2 <- dat2[-c(15:17), ]#remove Heterixalus
spp <-dat2$Species

Mass2 <- dat2$mass_g / 1000 #convert to kilograms
musc.mass <- dat2$muscle_mass_percent * Mass2
V2 <- dat2$V_svl_s * dat2$SVL_mm / 1000
E2 <- dat2$E_J_kg_muscle * musc.mass
P2 <- dat2$P_W_kg_muscle * musc.mass

#remove NAs
V2 <- V2[complete.cases(V2)]
E2 <- E2[complete.cases(E2)]
P2 <- P2[complete.cases(P2)]

names(V2) <- dat2$Species
names(E2) <- dat2$Species
names(P2) <- dat2$Species

#Mean species proxies
L2 <- dat2$leg_svl * dat2$SVL_mm / 1000
CSA2 <- musc.mass ^ (2/3)

V2.proxy <- sqrt(2 * L2 * CSA2 / Mass2)
E2.proxy <- L2 * CSA2
P2.proxy <- sqrt(2 * L2 * (CSA2 ^ 3) / Mass2)

V2.proxy <- V2.proxy[complete.cases(V2.proxy)]
E2.proxy <- E2.proxy[complete.cases(E2.proxy)]
P2.proxy <- P2.proxy[complete.cases(P2.proxy)]

names(V2.proxy) <- dat2$Species
names(E2.proxy) <- dat2$Species
names(P2.proxy) <- dat2$Species

#Stick into data frame and remove rows with NAs
df.V <- rrpp.data.frame(V.proxy = c(V.proxy, V2.proxy), V = c(V, V2), Species = c(as.character(Species), as.character(dat2$Species)))
df.E <- rrpp.data.frame(E.proxy = c(E.proxy, E2.proxy), E = c(E, E2), Species = c(as.character(Species), as.character(dat2$Species)))
df.P <- rrpp.data.frame(P.proxy = c(P.proxy, P2.proxy), P = c(P, P2), Species = c(as.character(Species), as.character(dat2$Species)))


lm.V <- lm.rrpp(V ~ V.proxy/Species, data = df.V, iter = 9999)
anova(lm.V, effect.type = "F", error = c("V.proxy:Species", "Residuals"))

lm.E <- lm.rrpp(E ~ E.proxy/Species, data = df.E, iter = 9999)
anova(lm.E, effect.type = "F", error = c("E.proxy:Species", "Residuals"))

lm.P <- lm.rrpp(P ~ P.proxy/Species, data = df.P, iter = 9999)
anova(lm.P, effect.type = "F", error = c("P.proxy:Species", "Residuals"))

# Sex-specific predictions --------------------------------------------------------------------

svl <- tapply(SVL, list(Species, Sex), mean)[complete.cases(tapply(SVL, list(Species, Sex), mean)),]
mass <- tapply(Mass, list(Species, Sex), mean)[complete.cases(tapply(Mass, list(Species, Sex), mean)),]

V <- tapply(v, list(Species, Sex), mean)[complete.cases(tapply(v, list(Species, Sex), mean)),]
V <- melt(V / svl, varnames = c("Species", "Sex"))

E <- tapply(0.5 * Mass * v ^ 2, list(Species, Sex), mean)[complete.cases(tapply(0.5 * Mass * v ^ 2, list(Species, Sex), mean)),]
E <- melt(E / mass, varnames = c("Species", "Sex"))

P <- tapply(P.Mass * Mass, list(Species, Sex), mean)[complete.cases(tapply(P.Mass * Mass, list(Species, Sex), mean)),]
P <- melt(P / mass, varnames = c("Species", "Sex"))

V.proxy <- tapply(2 * L * FPCSA / Mass, list(Species, Sex), mean)[complete.cases(tapply(2 * L * FPCSA / Mass, list(Species, Sex), mean)),]
V.proxy <- melt(V.proxy / svl, varnames = c("Species", "Sex"))

E.proxy <- tapply(0.5 * Mass * v.proxy ^ 2, list(Species, Sex), mean)[complete.cases(tapply(0.5 * Mass * v.proxy ^ 2, list(Species, Sex), mean)),]
E.proxy <- melt(E.proxy / mass, varnames = c("Species", "Sex"))

P.proxy <- tapply(Mass * v.proxy ^ 3 / (2 * L), list(Species, Sex), mean)[complete.cases(tapply(Mass * v.proxy ^ 3 / (2 * L), list(Species, Sex), mean)),]
P.proxy <- melt(P.proxy / mass, varnames = c("Species", "Sex"))

df.V <- rrpp.data.frame(V.proxy = V.proxy$value, V = V$value, Species = V$Species, Sex = V$Sex)
df.E <- rrpp.data.frame(E.proxy = E.proxy$value, E = E$value, Species = E$Species, Sex = V$Sex)
df.P <- rrpp.data.frame(P.proxy = P.proxy$value, P = P$value, Species = P$Species, Sex = P$Sex)

lm.rrpp.V <- lm.rrpp(V.proxy ~ V * Sex, iter = 9999, SS.type = "II", data = df.V)
anova(lm.rrpp.V) #intx of proxy and sex insignificant

lm.rrpp.E <- lm.rrpp(E.proxy ~ E * Sex, iter = 9999, SS.type = "II", data = df.E)
anova(lm.rrpp.E) #intx of proxy and sex insignificant

lm.rrpp.P <- lm.rrpp(P.proxy ~ P * Sex, iter = 9999, SS.type = "II", data = df.P)
anova(lm.rrpp.P) #intx of proxy and sex insignificant



# VertLife and consensus tree -----------------------------------------------------

library(ape)
post.tree.all <- read.tree("/Users/Bryan/Downloads/amph_shl_new_Consensus_7238.tre") #load Jetz and Pyron 2018 tree

dat <- read.csv("/Users/Bryan/Dropbox/frogs/data/Moen et al. 2013/MoenFrogs.csv")
dat2 <- read.csv("/Users/Bryan/Dropbox/frogs/Data/interspecific_proxies/Data_Astley_2016.csv")
spp1 <- sort(unique(c(as.character(dat$Species), as.character(dat2$Species)))) # make list of all species in Moen et al. 2013 and Astley 2016 dataset

#make following substitutions for data to match Jetz and Pyron 2018
spp1[10] <- "Syncope_bassleri" 
spp1[27] <- "Rana_pipiens"
spp1[29] <- "Cyclorana_australis"
spp1[34] <- "Cyclorana_longipes"

match(spp1, post.tree.all$tip.label)
cat(gsub("_", " ", spp1), file = "spp_list_posterior.txt", sep = "\n") #write text file with names to submit to vertlife.org


#load trees from VertLife and obtain consensus tree

post.trees <- read.nexus("/Users/Bryan/Downloads/vert_life_tree-pruner-ecb9d133-5b08-4795-8584-f1c8ade20c4f/output.nex")
con.tree <- consensus(post.trees, p = 0.5, check.labels = TRUE)
setwd("/Users/Bryan/Dropbox/frogs/Data/interspecific_proxies/")
write.nexus(con.tree, file = "consensus.tree.nex")
