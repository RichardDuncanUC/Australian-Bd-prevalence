
# Bd Prevalence Analysis in Australian Frogs

# Load required libraries
library(tidyverse)
library(lubridate)
library(terra)
library(sf)
library(jagsUI)
library(data.table)
library(ape)
library(caper)

# read in helper functions
source("./RScripts/Helper functions.R")

#-------------------------------------------------------------------------------
# Read spatial and prevalence data
#-------------------------------------------------------------------------------
# Read Australian mainland shapefile
mainvec <- vect("./data/Australia shape vector.shp")

# Read in tidy Bd prevalence data
# these data were generated from the data in Murray et al. 2010
# downloaded from https://figshare.com/collections/The_distribution_and_host_range_of_the_pandemic_disease_chytridiomycosis_in_Australia_spanning_surveys_from_1956_2007/3302961
# and data in the Amphibian Disease portal
# downloaded from https://amphibiandisease.org/query/
# specifying disease = Bd, for a box bounding Australia
# the data were filtered to exclude uncertain records, nomenclature was standardised and duplicates were removed

all.dat <- read.csv("./data/Tidy Bd prevalence data.csv", strip.white = T) |>
  glimpse()

# column names:
# species = frog species tested for Bd
# year = year of Bd testing
# lat = latitude of sample location
# long = longitude of sample location
# family = frog species family
# n.ind = number of individuals tested for Bd in the sample
# n.pos = number of Bd positive individuals in the sample
# decline = declined due to Bd (yes) or not (no), taken from Scheele et al. 2017
# temp.overall = the overall (average) mean annual temperature at the sample location
# temp.year = the mean annual temperature at the sample location in the year of sampling
# derived from maximum and minimum temperature data over time for Australia downloaded from
# https://s3-ap-southeast-2.amazonaws.com/silo-open-data/Official/annual/min_temp/
# https://s3-ap-southeast-2.amazonaws.com/silo-open-data/Official/annual/max_temp/

# Check for sample duplicates
# by creating an identifier for each sample location
all.dat$loc <- paste(all.dat$long, all.dat$lat)
table(duplicated(cbind(all.dat$species, all.dat$loc, all.dat$year)))

# Calculate total individuals sampled and number of locations sampled per species
all.dat <- all.dat |>
  group_by(species) |>
  mutate(total.ind = sum(n.ind),
         total.loc = n()) |>
  glimpse()

# Create species list
sp.list <- levels(factor(all.dat$species))

#-------------------------------------------------------------------------------
# Preliminary Data Analysis
#-------------------------------------------------------------------------------
# First year of infection
first.year <- min(all.dat$year[all.dat$n.pos > 0])
first.year

# Year range
range(all.dat$year)

# Infection prevalence over time
pt <- all.dat |>
  group_by(year) |>
  summarise(nind = sum(n.ind),
            npos = sum(n.pos)) |>
  mutate(prop = npos/nind) |>
  glimpse()

# Prevalence per species
ps <- all.dat |>
  mutate(prev = n.pos / n.ind) |>
  group_by(species) |>
  summarise(mean.prev = mean(prev),
            n = n(),
            n.ind = sum(n.ind)) |>
  glimpse() 

#-------------------------------------------------------------------------------
# Visualise Infection Patterns
#-------------------------------------------------------------------------------
pdf("./Figures/Figure 2.pdf")

par(mfrow = c(2, 2))

# Spatial distribution of infections
plot(mainvec, ylab = "Latitude", xlab = "Longitude")
points(lat ~ long, data = filter(all.dat, n.pos == 0), pch = 16, col = rgb(0, 0, 1, 0.3), cex = 1)
points(lat ~ long, data = filter(all.dat, n.pos > 0), pch = 16, col = rgb(1, 0, 0, 0.3), cex = 1)
mtext("A", adj = 0)

# Mean prevalence per species
par(mar = c(4, 4, 1, 1))
hist(ps$mean.prev, breaks = seq(0, 1, 0.1), xlab = "Mean prevalence", main = "",
     cex.lab = 1.2)
text(0.8, 12, "Species, n = 42", cex = 1.2, xpd = NA)
mtext("C", adj = 0, xpd = NA)

# Infection prevalence over time
# scale the size of points
ss <- log(pt$nind) / 4
ss <- sc.range(pt$nind, old.min = min(pt$nind), old.max = max(pt$nind), new.min = 0.5, new.max = 2.5)

plot(prop ~ year, pch = 19, data = pt, cex = ss, bty = "l", 
     ylab = "Mean prevalence", xlab = "Year", cex.lab = 1.2)

# add legend for point sizes
p <- c(10, 100, 1000)
lp <- sc.range(p, old.min = min(pt$nind), old.max = max(pt$nind), new.min = 0.5, new.max = 2.5)
legend(2000, 1.1, legend = c("10", "100", "1000"), pch = 19, pt.cex = lp, cex = 0.8,
       title = "Number of individuals", xpd = NA)
mtext("B", adj = 0, xpd = NA)

# Prevalence per sample location
hist(all.dat$n.pos / all.dat$n.ind, breaks = seq(0, 1, 0.1), xlab = "Prevalence", main = "",
     cex.lab = 1.2)
text(0.8, 800, "Samples, n = 1455", cex = 1.2, xpd = NA)
mtext("D", adj = 0, xpd = NA)

dev.off()

#-------------------------------------------------------------------------------
# Data Preparation for Temperature Analysis
#-------------------------------------------------------------------------------
# Use annual temperature in each year for analysis
all.dat$temp <- all.dat$temp.year

#-------------------------------------------------------------------------------
# Descriptive Statistics
#-------------------------------------------------------------------------------
# Total individuals sampled
sum(all.dat$n.ind)

# Total positive samples
sum(all.dat$n.pos)
sum(all.dat$n.pos) / sum(all.dat$n.ind)

# Number of species
length(table(all.dat$species))

# Number of sampling events
nrow(all.dat)

# Number of locations
length(table(all.dat$loc))

# Summary by species
sum.dat <- all.dat |>
  group_by(species) |>
  summarise(n.loc = n(),
            n.ind = sum(n.ind)) |>
  arrange(-n.ind)

sum.dat

#-------------------------------------------------------------------------------
# Prepare Data for JAGS Model
#-------------------------------------------------------------------------------
y <- all.dat$n.pos                           # number of Bd positive individuals
n.ind <- all.dat$n.ind                       # total number of individuals tested
spp <- as.numeric(factor(all.dat$species))   # indicator for species
loc <- as.numeric(factor(all.dat$loc))       # indicator for sample location

temp <- all.dat$temp                         # mean annual temperature per location-year
temp[is.na(temp)] <- mean(temp, na.rm = T)   # one missing value - set to the mean

# centre and scale temperature
sc.temp <- sc(temp)

N <- length(y)
N.spp <- max(spp)
N.loc <- max(loc)

# continuous value for year scaled
year <- sc(all.dat$year)
table(year)

# indicator for year
fac.year <- as.numeric(factor(year))
N.year <- max(fac.year)

# observation level random effect to account for overdispersion
od <- 1:length(y)

#-------------------------------------------------------------------------------
# Read fitted temperature curves
t.curve <- read.csv("./data/Fitted frog temperature curves.csv")
b.curve <- read.csv("./data/Fitted Bd temperature curves.csv")

# Bd thermal affinity
med.temp.chy <- b.curve$med.temp

# calculate log-transformed thermal mismatch ratio (LTMR) and LTMR:MAT slope for each species
# limit temperature range to 90% quantiles
# tmr_func is a function in the Helper function.R file

tmr <- matrix(nrow = length(sp.list), ncol = 2)

for(k in 1:length(sp.list)) {
  tmr[k, ] <- tmr_func(q1 = t.curve$q1.temp[t.curve$species == sp.list[k]],
                       q9 = t.curve$q9.temp[t.curve$species == sp.list[k]],
                       a = t.curve$a[t.curve$species == sp.list[k]],
                       b = t.curve$b[t.curve$species == sp.list[k]],
                       c = t.curve$c[t.curve$species == sp.list[k]],
                       med.temp = t.curve$med.temp[t.curve$species == sp.list[k]],
                       plot = F)
}

# calculate thermal affinity difference = species thermal affinity - Bd thermal affinity
ta_dif <- t.curve$med.temp - b.curve$med.temp

# create species level dataset
frog <- t.curve
frog$tmr <- tmr[, 1]           # thermal mismatch ratio
frog$tmr_slope <- tmr[, 2]     # slope of TMR:MAT relationship
frog$ta_dif <- ta_dif          # thermal affinity difference

glimpse(frog)

#-------------------------------------------------------------------------------
# Add decline status and number of individuals sampled per species
#-------------------------------------------------------------------------------
dec <- all.dat |>
  dplyr::select(species, decline, n.ind) |>
  group_by(species, decline) |>
  summarise(n.ind = sum(n.ind)) |>
  glimpse()

frog <- left_join(frog, dec) |>
  glimpse()

#-------------------------------------------------------------------------------
# use thermal affinity difference in model and put on same scale as temp
sc.ta_dif <- (frog$ta_dif - mean(temp)) / sd(temp)

#-------------------------------------------------------------------------------
# JAGS Model Specification
# code is linked to equation numbers in the text
#-------------------------------------------------------------------------------
mod <- "model {
  for(i in 1:N) {
    # likelihood for prevalence
    y[i] ~ dbin(p[i], n.ind[i])                                                   # eq 2
    logit(p[i]) <- int + b.temp[spp[i]]*temp[i] + b.year.trend*year[i] +
                   b.year[fac.year[i]] + b.spp[spp[i]] + b.loc[loc[i]] + b.od[i]  # eq 3
    
    # overdispersion
    b.od[i] ~ dnorm(0, tau.od)                                                    # eq 9
  }

  for(i in 1:N.spp) {
    b.spp[i] ~ dnorm(0, tau.spp)                                                  # eq 4
    b.temp[i] ~ dnorm(mu.temp[i], tau.temp)                                       # eq 5
    mu.temp[i] <- int.temp + slope.temp * ta_dif[i]                               # eq 6
  }

  for(i in 1:N.loc) {
    b.loc[i] ~ dnorm(0, tau.loc)                                                  # eq 7
  }
  
  for(i in 1:N.year) {
    b.year[i] ~ dnorm(0, tau.year)                                                # eq 8
  }

  # priors
  int ~ dnorm(0, 0.7)         # close to uniform prior on overall prevalence on probability scale
  int.temp ~ dnorm(0, 0.01)
  slope.temp ~ dnorm(0, 0.01)
  b.year.trend ~ dnorm(0, 0.01)

  tau.spp <- 1/(sigma.spp * sigma.spp)
  sigma.spp ~ dunif(0, 10)

  tau.loc <- 1/(sigma.loc * sigma.loc)
  sigma.loc ~ dunif(0, 10)
  
  tau.year <- 1 / (sigma.year * sigma.year)
  sigma.year ~ dunif(0, 10)

  tau.od <- 1 / (sigma.od * sigma.od)
  sigma.od ~ dunif(0, 10)

  tau.temp <- 1 / (sigma.temp * sigma.temp)
  sigma.temp ~ dunif(0, 10)

}"  

write(mod, "model.txt")

#-------------------------------------------------------------------------------
# Model Execution
#-------------------------------------------------------------------------------
# Sample model in JAGS with three chains
set.seed(234)
mod.jags <- jags(model = "model.txt",
                  data = list(y = y, n.ind = n.ind, N = N, spp = spp, loc = loc, N.spp = N.spp, 
                              N.loc = N.loc, year = year, fac.year = fac.year, N.year = N.year,
                              temp = sc.temp, ta_dif = sc.ta_dif),
                  param = c("int", "int.temp", "slope.temp", "b.year.trend", "b.year", 
                            "sigma.od", "sigma.year", "sigma.spp", "sigma.loc", "sigma.temp", 
                            "b.temp", "b.spp","b.loc"),
                  n.chains = 3,
                  n.iter = 12000,
                  n.burnin = 2000,
                  parallel = T)

# Extract model summary
jags.sum <- mod.jags$summary
jags.sum[1:50, c(1, 2, 3, 7, 8, 9)]
tail(jags.sum)[, 1:9]

#-------------------------------------------------------------------------------
# Parameter Analysis
#-------------------------------------------------------------------------------
# Standard deviation parameters
sigma.spp <- mod.jags$sims.list$sigma.spp         # species-level intercepts
sigma.loc <- mod.jags$sims.list$sigma.loc         # sample locations
sigma.year <- mod.jags$sims.list$sigma.year       # sample years
sigma.od <- mod.jags$sims.list$sigma.od           # overdispersion
sigma.temp <- mod.jags$sims.list$sigma.temp       # species-level slopes

quantile(sigma.spp, c(0.025, 0.5, 0.975))
quantile(sigma.loc, c(0.025, 0.5, 0.975))
quantile(sigma.year, c(0.025, 0.5, 0.975))
quantile(sigma.od, c(0.025, 0.5, 0.975))
quantile(sigma.temp, c(0.025, 0.5, 0.975))

#-------------------------------------------------------------------------------
# Prevalence Estimation
#-------------------------------------------------------------------------------
# Extract model parameters
overall.int <- mod.jags$sims.list$int
spp.int <- mod.jags$sims.list$b.spp
slope <- mod.jags$sims.list$b.temp
b.year.trend <- mod.jags$sims.list$b.year.trend
b.year <- mod.jags$sims.list$b.year
int.temp <- mod.jags$sims.list$int.temp
slope.temp <- mod.jags$sims.list$slope.temp

# get median values
overall.int <- median(overall.int)
spp.int <- apply(spp.int, 2, median)
slope <- apply(slope, 2, median)
b.year.trend <- median(b.year.trend)
b.year <- apply(b.year, 2, median)
int.temp <- median(int.temp)
slope.temp <- median(slope.temp)

# Find median year for prediction
med.yr <- median(rep(all.dat$year, all.dat$n.ind))
level.yr <- which(as.numeric(levels(factor(all.dat$year))) == med.yr)

# Convert to standardised year value
yr <- unique(year[all.dat$year == med.yr])

# Range of standardised temperature values
tmp <- seq(0, 32, 0.01)
# rescale so temp values are on the scale used in the model
s <- sd(temp)
m <- mean(temp)
tt <- (tmp - m) / s 

# estimate prevalence (prev) for each species at temperatures tt from the fitted model
prev <- matrix(nrow = length(tt), ncol = length(slope))
m.prev <- numeric()
for(i in 1:length(sp.list)) {
  prev[, i] <- overall.int + spp.int[i] + slope[i] * tt + b.year.trend * yr + b.year[level.yr]
  
  # subset estimated prevalence and temperature to temperature range for species
  mp <- prev[tmp >= frog$min.temp[i] & tmp <= frog$max.temp[i], i]
  stmp <- tmp[tmp >= frog$min.temp[i] & tmp <= frog$max.temp[i]]
  # get the fitted occupancy distribution for that species across the temperature range
  yy <- genlog(stmp, as.numeric(frog[i, 1:3]))
  # get mean estimated prevalence, weighted by occupancy 
  m.prev[i] <- weighted.mean(mp, w = yy)
}

# add species-level estimated prevalence to species file
frog$est.prev <- m.prev
# add estimated within-species prevalence-thermal affinity slope for each species
frog$est.slope <- slope

#-------------------------------------------------------------------------------
# Visualise Results
#-------------------------------------------------------------------------------
pdf("./Figures/Figure 3.pdf", width = 7, height = 8)

par(mfrow = c(3, 2), mar = c(4, 4, 2, 1))

#-------------------------------------------------------------------------------
# Plot raw prevalence:MAT relationship
# estimate prevalence in 1 degree bins
sp <- all.dat |>
  mutate(grp = cut(temp, breaks = seq(0, 32, 1))) |>
  group_by(grp) |>
  summarise(n.pos = sum(n.pos),
            n.ind = sum(n.ind),
            temp = mean(temp)) |>
  mutate(prop = n.pos / n.ind) |>
  glimpse()

# plot with point size proportional to number of individuals sampled
sx <- sc.range(sp$n.ind, old.min = min(sp$n.ind), old.max = max(sp$n.ind), new.min = 0.5, new.max = 2.5)

plot(prop ~ temp, data = sp, pch = 19, bty = "l", cex = sx, 
     ylab = "Mean prevalence", xlab = "Mean annual temperature (\u00B0C)")
p <- c(20, 200, 2000)
lp <- sc.range(p, old.min = min(sp$n.ind), old.max = max(sp$n.ind), new.min = 0.5, new.max = 2.5)
legend(18, 1.1, legend = c("20", "200", "2000"), pch = 19, pt.cex = lp, cex = 0.8,
       title = "Number of individuals", xpd = NA)
mtext("A", adj = 0)

#-------------------------------------------------------------------------------
# Plot estimated prevalence across the temperature range for each species
# on the probability scale

# sens is esitmated prevalence at median temp (thermal affinity)
sens <- numeric()

# for first species
# temperature range of species
stmp <- tmp[tmp >= frog$min.temp[1] & tmp <= frog$max.temp[1]]
# prevalence values across that range on the probability scale
sprev <- unlogit(prev[, 1][tmp > frog$min.temp[1] & tmp < frog$max.temp[1]])

plot(sprev ~ stmp, type = "l", ylim = c(0, 1), xlim = c(3, 30),
     bty = "l", xlab = "", ylab = "", col = "grey")
sens[1] <- sprev[which.min(abs(frog$med.temp[1] - stmp))]
points(frog$med.temp[1], sens[1], pch = 19, col = "blue", xpd = NA)
title(ylab = "Prevalence", xlab = "Mean annual temperature (\u00B0C)", line = 2.5)
mtext("B", adj = 0)

# subsequent species
for(i in 2:nrow(frog)) {
  stmp <- tmp[tmp >= frog$min.temp[i] & tmp <= frog$max.temp[i]]
  sprev <- unlogit(prev[, i][tmp >= frog$min.temp[i] & tmp <= frog$max.temp[i]])
  lines(sprev ~ stmp, col = "grey")
  sens[i] <- sprev[which.min(abs(frog$med.temp[i] - stmp))]
  points(frog$med.temp[i], sens[i], pch = 19, col = "blue", xpd = NA)
}

#-------------------------------------------------------------------------------
# slope as a function of mean temperature

# line fitted in the model
pred.pt.slope <- int.temp + slope.temp * sc.ta_dif

# r2 of fit
rss <- sum((frog$est.slope - pred.pt.slope)^2)
ss <- sum((frog$est.slope - mean(frog$est.slope))^2)
r2 <- round(1 - rss/ss, 2)
eq <- substitute(italic(r)^2~"="~r2~", "~italic(P)~" <0.001", 
         list(r2 = r2))

plot(frog$est.slope ~ sc.ta_dif, pch = 16, bty = "l",
     ylab = "", xlab = "", xaxt = "n", cex = 1.3, col = "blue", main = eq)
abline(int.temp, slope.temp, lwd = 2)
abline(h = 0, lty = 3)
abline(v = (0 - m)/s, col = "red", lwd = 2)

# add x axis
locx <- seq(-5, 15, 5)
loc.sc <- (locx - m) / s
axis(1, at = loc.sc, labels = locx)
#title(eq, line = 0, cex = 0.9)
title(ylab = "Slope prevalence:MAT", xlab = "Thermal affinity difference (\u00B0C)", line = 2.5)
mtext("C", adj = 0)

#-------------------------------------------------------------------------------
# mean prevalence versus absolute thermal affinity
# fitted relationship
m1 <- lm(frog$est.prev ~ abs(frog$ta_dif))
r2 <- round(summary(m1)$r.squared, 2)
eq <- substitute(italic(r)^2~"="~r2~", "~italic(P)~" <0.001", 
         list(r2 = r2))

plot(frog$est.prev ~ abs(frog$ta_dif), pch = 19, cex = 1.3, bty = "l", yaxt = "n", 
     ylab = "Mean prevalence", xlab = "Absolute thermal affinity difference (\u00B0C)",
     col = "blue", main = eq)
axis(2, at = logit(c(0.001, 0.01, 0.1, 0.4, 0.8)), labels = c("0.001", "0.01", "0.1", "0.4", "0.8"))
abline(m1, lwd = 2)
abline(v = 0, col = "red", lwd = 2)
mtext("D", adj = 0)

#-------------------------------------------------------------------------------
# observed versus predicted

par(mar = c(6, 6, 2, 1))

m1 <- lm(frog$est.slope ~ frog$tmr_slope)
r2 <- round(summary(m1)$r.squared, 2)
eq <- substitute(italic(r)^2~"="~r2~", "~italic(P)~" <0.05", 
         list(r2 = r2))
plot(frog$est.slope ~ frog$tmr_slope, pch = 19, cex = 1.3, bty = "l", 
     ylab = "Observed\n Slope Prevalence:MAT", 
     xlab = "Predicted\n Slope LTMR:MAT", main = eq, col = "blue")
abline(m1, lwd = 2)
mtext("E", adj = 0)

# estimated mean prevalence versus predicted
m1 <- lm(frog$est.prev ~ frog$tmr)
r2 <- round(summary(m1)$r.squared, 2)
eq <- substitute(italic(r)^2~"="~r2~", "~italic(P)~" <0.001", 
         list(r2 = r2))

plot(frog$est.prev ~ frog$tmr, pch = 19, cex = 1.3, bty = "l", 
     ylab = "Osbserved\n Mean prevalence", xlab = "Predicted\n Mean LTMR", main = eq, col = "blue")
mtext("F", adj = 0)
abline(m1, lwd = 2)

dev.off()

#-------------------------------------------------------------------------------
# Trait Data Analysis
#-------------------------------------------------------------------------------
# Read in trait data
trait.dat <- read.csv("./data/Trait data.csv") |>
  glimpse()

# combine with frog data
frog.trait <- full_join(frog, trait.dat) |>
  glimpse()

#-------------------------------------------------------------------------------
# read in Australian frog phylogenetic tree 
# derived from Jetz and Pyron amphibian consensus tree: amph_shl_new_Consensus_7238.tre

aus_tree <- read.tree("./data/Australian frog tree.tre")
# adjust tip labels
aus_tree$tip.label <- gsub("_", " ", aus_tree$tip.label)

#-------------------------------------------------------------------------------
# model to assess factors that predict prevalence 
# accounting for phylogenetic relationships 

# construct terrestriality index by summing terretriality values for eggs, tadpoles and adults
eggs <- ifelse(frog.trait$eggs == "aquatic", 0, 1)
tadpoles <- ifelse(frog.trait$tadpoles == "aquatic", 0, 1)
guild <- as.numeric(factor(frog.trait$ecological_guild, levels = c("wetland", "terrestrial", "arboreal"))) - 1

frog.trait$terr_index <- eggs + tadpoles + guild
table(frog.trait$terr_index)

# convert to format for use in caper package
cd <- comparative.data(aus_tree, frog.trait, names.col = "species")

#-------------------------------------------------------------------------------
# Figures

pdf("./Figures/Figure 4.pdf")

layout(matrix(c(1, 1, 2, 3, 3, 2), ncol = 3, byrow = T), widths = c(0.7, 0.7, 1))
par(mar = c(5, 5, 2, 1))

# fit phylogenetic regression model with slope of prevalence:MAT as a function of covariates
m1 <- pgls(est.slope ~ sc(ta_dif) + sc(terr_index) + sc(female_max_svl) + sc(precip), data = cd, lambda = "ML")
summary(m1)
r2 <- round(summary(m1)$r.squared, 2)
eq <- substitute("Slope prevalence:MAT   "~italic(r)^2~"="~r2, 
         list(r2 = r2))

out <- summary(m1)$coef
pa <- data.frame(est = out[2:5, 1],
                 lcl = out[2:5, 1] - 2*out[2:5, 2],
                 ucl = out[2:5, 1] + 2*out[2:5, 2])

xx <- 1:4
yl <- c(min(pa$lcl), max(pa$ucl))
plot(pa$est ~ xx, bty = "l", pch = 19, cex = 2, ylim = yl, xaxt = "n", xlim = c(0.5, 4.5),
     xlab = "", ylab = "Effect size", main = eq)
  arrows(xx, pa$lcl, xx, pa$ucl, length = 0, lwd = 2)
  axis(1, at = 1:4, labels = c("Thermal affinity", "Terrestrial", "SVL", "Precipitation"))
  abline(h = 0, lty = 3)
  mtext("A", adj = 0)

#-------------------------------------------------------------------------------  
# Decline due to Bd or not
cl <- ifelse(frog$decline == "yes", rgb(1, 0, 0, 0.7), rgb(0.1, 0.1, 0.1, 0.7))
cx <- sc.range(frog$n.ind, old.min = min(frog$n.ind), old.max = max(frog$n.ind), new.min = 0.5, new.max = 2.5)

# fit phylogenetic ANOVA model
m1 <- pgls(unlogit(est.prev) ~ decline, data = cd, lambda = "ML")
summary(m1)
pval <- round(summary(m1)$coef[2, 4], 3)
eq <- substitute(italic(P)~" = "~pval, list(pval = pval))

par(bty = "l")
plot(unlogit(est.prev) ~ factor(decline), data = frog, col = "white",
     xlab = "Decline", ylab = "Mean prevalence", main = eq)
  points(unlogit(est.prev) ~ ifelse(decline == "no", 1, 2), data = frog, 
         pch = 19, col = cl, cex = cx)
  mtext("C", adj = 0)

#-------------------------------------------------------------------------------
# mean prevalence model
# fit phylogenetic regression model with mean prevalence as a function of covariates
m1 <- pgls(est.prev ~ sc(abs(ta_dif)) + sc(terr_index) + sc(female_max_svl) + sc(precip), data = cd, lambda = "ML")
summary(m1)
r2 <- round(summary(m1)$r.squared, 2)
eq <- substitute("Mean prevalence    "~italic(r)^2~"="~r2, 
         list(r2 = r2))

out <- summary(m1)$coef
pa <- data.frame(est = out[2:5, 1],
                 lcl = out[2:5, 1] - 2*out[2:5, 2],
                 ucl = out[2:5, 1] + 2*out[2:5, 2])

xx <- 1:4
yl <- c(min(pa$lcl), max(pa$ucl))
plot(pa$est ~ xx, bty = "l", pch = 19, cex = 2, ylim = yl, xaxt = "n", xlim = c(0.5, 4.5),
     xlab = "", ylab = "Effect size", main = eq)
  arrows(xx, pa$lcl, xx, pa$ucl, length = 0, lwd = 2)
  axis(1, at = 1:4, labels = c("Thermal affinity", "Terrestrial", "SVL", "Precipitation"))
  abline(h = 0, lty = 3)
  mtext("B", adj = 0)
  
dev.off()

