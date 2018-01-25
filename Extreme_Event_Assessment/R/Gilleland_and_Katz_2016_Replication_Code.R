library("extRemes")

# Replicate Fig. 1.

data("PORTw", package = "extRemes")
plot(PORTw$TMX1, type = "l", xlab = "Year",
     ylab = "Maximum winter temperature",
     col = "darkblue", lwd = 1.5, cex.lab = 1.25)

## Fitting the GEV to block maxima data.

# maximum winter temperatures at Port Jervis.
fit1 <- fevd(TMX1, PORTw, units = "deg C")

# printout of the estimated parameters (cf. Table 1).
fit1

# to get a named vector with just the parameters, 
# negative log-likelihood value,
# and parameter covariance.
distill(fit1)

# Some useful plots.

# Fig. 2
plot(fit1)

# Fig. 3
plot(fit1, "trace")

return.level(fit1)
return.level(fit1, do.ci = TRUE)

ci(fit1, return.period = c(2, 20, 100))
ci(fit1, type = "parameter")
ci(fit1, type = "parameter", which.par = 3, xrange = c(-0.4, 0.01), 
   nint = 100, method = "proflik", verbose = TRUE)

# Fig. 4
ci(fit1, method = "proflik", xrange = c(22, 28), verbose = TRUE)

# Fitting the Gumbel to the maximum winter temperatures.

fit0 <- fevd(TMX1, PORTw, type = "Gumbel", units = "deg C")
fit0

plot(fit0)
plot(fit0, "trace")

# Test of null hypothesis that the shape parameter is zero.
# The likelihood-ratio test and checking the confidence intervals
# for fit1.

lr.test(fit0, fit1)

ci(fit1, type = "parameter", which.par = 3, xrange = c(-0.4, 0.01),
   nint = 100, method = "proflik", verbose = TRUE)
ci(fit1, method = "proflik", xrange = c(22, 28), verbose = TRUE)

pextRemes(fit1, q = c(20, 25, 28.8, 28.9), lower.tail = FALSE)

# Drawing a random sample from the fitted model.

z <- rextRemes(fit1, 100)

fitz <- fevd(z)

# cf. with parameter estimates in Table 1.
fitz 

plot(fitz)

# Using Bayesian estimation.

fB <- fevd(TMX1, PORTw, method = "Bayesian")
fB

postmode(fB)

# cf. with Fig. 2 for the MLE fit.
plot(fB)

# Fig. 5
plot(fB, "trace")

# 100-year return level with 95% credible interval.
ci(fB) 

# 95% credible intervals for the parameters.
ci(fB, type = "parameter")


# Using L-moments estimation.
fitLM <- fevd(TMX1, PORTw, method = "Lmoments")
fitLM

# cf. with Fig. 2 for the MLE fit.  
plot(fitLM) 

# the following employ the parametric bootstrap procedure.
ci(fitLM)
ci(fitLM, type = "parameter")

fit2 <- fevd(TMX1, PORTw, location.fun = ~AOindex, units = "deg C")
fit2

# Fig. 6
plot(fit2)
plot(fit2, "trace")

# Test whether inclusion of the AO index is statistically significant
# using the likelihood-ratio test.

lr.test(fit1, fit2)

# Estimate "effective" return levels for values of AO index equal to -1 and 1.

v <- make.qcov(fit2, vals = list(mu1 = c(-1, 1)))
return.level(fit2, return.period = c(2, 20, 100), qcov = v)


# Add covariate to the scale parameter.
fit3 <- fevd(TMX1, PORTw, location.fun = ~AOindex, scale.fun = ~AOindex,
             units = "deg C")
lr.test(fit2, fit3)

# same as above, but use a log-linke function for the scale parameter.
fit4 <- fevd(TMX1, PORTw, location.fun = ~AOindex, scale.fun = ~AOindex,
    use.phi = TRUE, units = "deg C")

# Sometimes, it might be helpful to run the optimization using the log(scale)
# instead of the raw scale.
fevd(TMX1, PORTw, use.phi = TRUE)

## Fitting the GP df to data.

# Hurricane damage data
data("damage", package = "extRemes")

# Fig. 7
par(mfrow = c(2, 2))
plot(damage$Year, damage$Dam, xlab = "Year",
     ylab = "U.S. Hurricane Damage (billion USD)",
     cex = 1.25, cex.lab = 1.25,
     col = "darkblue", bg = "lightblue", pch = 21)
plot(damage[, "Year"], log(damage[, "Dam"]), xlab = "Year",
     ylab = "ln(Damage)", ylim = c(-10, 5), cex.lab = 1.25,
     col = "darkblue", bg = "lightblue", pch = 21)
qqnorm(log(damage[, "Dam"]), ylim = c(-10, 5), cex.lab = 1.25)

# Threshold selection.
# Fig. 8

threshrange.plot(damage$Dam, r = c(1, 15), nint = 20)
mrlplot(damage$Dam, xlim = c(0, 12))

# Estimate the average number of hurricanes per year.
range(damage$Year)
1995 - 1926 + 1 # 70

# 144 hurricanes in 70 years.
dim(damage)

144 / 70 
# about 2.06 per year.

fitD <- fevd(Dam, damage, threshold = 6, type = "GP", time.units = "2.06/year")
fitD

plot(fitD)

# Fig. 9
plot( fitD, type = "qq")

pextRemes(fitD, c(20, 40, 60, 100), lower.tail = FALSE)

# Fort Collins Precipitation Data
data("Fort", package = "extRemes")

# Fig. 10
par(mfrow = c(1, 2))
plot(Prec ~ obs, data = Fort, type = "h", 
     xlab = "observation", ylab = "Precipitation (inches)")

plot(Prec ~ month, data = Fort, pch = 22,
     col = "darkblue", bg = "lightblue", 
     xlab = "Month", ylab = "")

id <- Fort$day == 29 & Fort$month == 7 & Fort$year == 1997
points(Fort$month[id], Fort$Prec[id], pch = 22, 
       col = "darkred", bg = "red")

# Fig. 11
threshrange.plot(Fort$Prec, r = c(0.01, 2), nint = 30)
fitFC <- fevd(Prec, Fort, threshold = 0.395, type = "GP")
fitFC
plot(fitFC)

# With the density plot obtained from the above line, note that
# conventional nonparametric density function estimators do not
# work well near a boundary.

fitFC2 <- fevd(Prec, Fort, threshold = 0.395, scale.fun = 
  ~ cos(2 * pi * tobs / 365.25) + sin(2 * pi * tobs / 365.25), 
  type = "GP", use.phi = TRUE, units = "inches")
fitFC2
plot(fitFC2)
lr.test(fitFC, fitFC2)

# Fig. 12
ci(fitFC2, type = "parameter", which.par = 2, method = "proflik",
   xrange = c(-0.5, 0.01), verbose = TRUE)

ci(fitFC2, type = "parameter", which.par = 3, method = "proflik",
   xrange = c(-0.1, 0.2), verbose = TRUE)

# 95% normal approximation confidence intervals for all parameters.
ci(fitFC2, type = "parameter")

v <- matrix(1, 730, 5)

v[, 2] <- cos(2 * pi * rep(1:365 / 365.25, 2))
v[, 3] <- sin(2 * pi * rep(1:365 / 365.25, 2))
v[, 5] <- 0.395 # Threshold.

v <- make.qcov(fitFC2, vals = v, nr = 730)
v[1:10, ]

FCprobs <- pextRemes(fitFC2, q = c(rep(1.54, 365), rep(4.63, 365)),
                     qcov = v, lower.tail = FALSE)

# Compute scale parameter for each day of the year (based on model).
phi <- -1.24022241 - 0.31235985 * cos(2 * pi * 1:365 / 365.25) +
    0.07733595 * sin(2 * pi * 1:365 / 365.25)

sigma <- exp(phi)

pevd(c(1.54, 4.63), scale = sigma[210], shape = 0.177584305, threshold = 0.395,
     type = "GP", lower.tail = FALSE)

u <- numeric(dim(Fort)[1])
u[is.element(Fort$month, c(1, 2, 11, 12))] <- 0.395
u[is.element(Fort$month, c(3, 4, 9, 10))] <- 0.25
u[is.element(Fort$month, 5:8)] <- 0.895

fitFC3 <- fevd(Prec, Fort, threshold = u,
               scale.fun = ~cos(2 * pi * tobs / 365.25) + sin(2 * pi * tobs / 365.25),
               type = "GP", verbose = TRUE)
fitFC3
plot(fitFC3)

## PP model fitting

fit <- fevd(Prec, Fort, threshold = 0.395, type = "PP",
            units = "inches", verbose = TRUE)
fit

# Fig. 13
plot(fit)
plot(fit, "trace")

ci(fit, type = "parameter")

# Phoenix Sky Harbor temperature data.
data("Tphap", package = "extRemes")

# Fig. 14
plot(Tphap$MinT, pch = 22, col = "darkblue", bg = "lightblue", 
     xlab = "time", ylab = "Phoenix Summer Minimum Temperature (deg. F)")

# Fig. 15
threshrange.plot(-Tphap$MinT, r = c(-75.5,-68.5), nint = 20, type = "PP")

fit1 <- fevd(-MinT ~1, Tphap, threshold = -73, type = "PP",
             units = "deg F", time.units = "62/year", verbose = TRUE)
fit1

# Recall that the fit is for the negative of MinT.
plot(fit1)
plot(fit1, "trace")

fit2 <- fevd(-MinT ~1, Tphap, threshold = c(-68,-7),
             threshold.fun = ~I((Year - 48)/42),
             location.fun = ~I((Year - 48)/42),
             type = "PP", time.units = "62/year", verbose = TRUE)
fit2

plot(fit2)

# Fig. 16
plot(fit2, type = "rl")

# Estimate "effective" return levels for 1952 and 1988
# and take the difference with confidence intervals.

v1952 <- make.qcov(fit2, vals = list(mu1 = (5 - 48)/42))
v1988 <- make.qcov(fit2, vals = list(mu1 = (41 - 48)/42))

return.level(fit2, return.period = 100, do.ci = TRUE,
             qcov = v1988, qcov.base = v1952)

fitFCpp <- fevd(Prec, Fort, threshold = 0.395, type = "PP")
fitFCpp2 <- fevd(Prec, Fort, threshold = 0.395,
                 location.fun = ~cos(2 * pi * tobs / 365.25) + sin(2 * pi * tobs / 365.25),
                 type = "PP", units = "inches")

fitFCpp3 <- fevd(Prec, Fort, threshold = 0.395,
                 location.fun = ~cos(2 * pi * tobs / 365.25) + sin(2 * pi * tobs / 365.25),
                 scale.fun = ~cos(2 * pi * tobs / 365.25) + sin(2 * pi * tobs / 365.25),
                 type = "PP", use.phi = TRUE, units = "inches")

lr.test(fitFCpp, fitFCpp2)
lr.test(fitFCpp2, fitFCpp3)

v <- make.qcov(fitFCpp3, vals = list(mu1 = cos(2 * pi * 1:365 / 365.25),
                                     mu2 = sin(2 * pi * 1:365 / 365.25),
                                     phi1 = cos(2 * pi * 1:365 / 365.25),
                                     phi2 = sin(2 * pi * 1:365/365.25)))

v[1:10, ] # Look at first ten rows of v.

p1.54 <- pextRemes(fitFCpp3, rep(1.54, 365), lower.tail = FALSE, qcov = v)

ciEffRL100 <- ci(fitFCpp3, return.period = 100, qcov = v)

# Fig. 17
plot(Fort$tobs, Fort$Prec, ylim = c(0, 10), xlab = "day",
     ylab = "Fort Collins Precipitation (inches)")

lines(ciEffRL100[, 1], lty = 2, col = "darkblue", lwd = 1.25)
lines(ciEffRL100[, 2], col = "darkblue", lwd = 1.25)
lines(ciEffRL100[, 3], lty = 2, col = "darkblue", lwd = 1.25)

legend("topleft", legend = c("Effective 100-year return level",
                             "95% CI (normal approx)"),
       col = "darkblue", lty = c(1, 2), lwd = 1.25, bty = "n")

## Dependent sequences.

# Fig. 18
atdf(-Tphap$MinT, 0.99) 

extremalindex(-Tphap$MinT, -68 - 7 * (Tphap$Year - 48)/42)

y <- -Tphap$MinT

u <- -68 - 7 *(Tphap$Year - 48)/42

grp <- Tphap$Year - 47

look <- decluster(y, threshold = u, groups = grp)
look

plot(look)

look2 <- decluster(y, threshold = u, method = "intervals", groups = grp)
look2

# Fig. 19
plot(look2)

# Fig. 20
atdf(look2, 0.99) 

u <- (Tphap$Year - 48)/42
u <- -68 - 7 * u

fitDC <- fevd(y, data = data.frame(y = c(look2), time = 1:2666),
              threshold = u, location.fun = ~time,
              type = "PP", time.units = "62/year")
fitDC

plot(fitDC)

par(mfrow = c(2, 2))
tmp <- eeplot(fitDC)
plot(tmp)

# Dependence in extremes between two variables.

z <- matrix(rnorm(2000), ncol = 2)

rho <- cbind(c(1, 0.9), c(0.9, 1))
rho

rho <- chol(rho)

t(rho) %*% rho

z <- t(rho %*% t(z))
cor(z[, 1], z[, 2])

plot(z[, 1], z[, 2])

taildep(z[, 1], z[, 2], 0.99)

taildep.test(z[, 1], z[, 2]) # Recall that the null hypothesis is tail dependence!

# Changing the options in the optimization routine.

# Use finite differences instead of the actual gradients.

fit <- fevd(-MinT ~1, Tphap, threshold = -73, type = "PP",
            units = "deg F", time.units = "62/year", use.phi = TRUE,
            optim.args = list(method = "BFGS", gr = NULL), verbose = TRUE)

# Super heavy tails

z <- revd(1000, loc = 2, scale = 1.5, shape = 0.2)
z <- exp(z)

hist(z)
hist(z[z < 1000])

threshrange.plot(z) # Lots of errors!  Will use u = exp(8).

fitZ <- fevd(z, threshold = exp(8), type = "GP")

threshrange.plot(log(z))

fitlZ <- fevd(log(z), threshold = 8, type = "GP")

par(mfrow = c(1, 2))
plot(fitZ, "qq2") # cf. Fig. 21 (left).
plot(fitlZ, "qq2") # cf. Fig. 21 (right).

