
#-------------------------------------------------------------------------------
# Custom helper functions
#-------------------------------------------------------------------------------
# Logit transformation (log-odds)
logit <- function(x) log(x / (1-x))

#-------------------------------------------------------------------------------
# Inverse logit transformation (probability)
unlogit <- function(x) exp(x) / (1 + exp(x))

#-------------------------------------------------------------------------------
# Scale function (center and standardize)
sc <- function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)

#-------------------------------------------------------------------------------
# Function to scale an old range to a new range
sc.range <- function(x, old.min, old.max, new.min, new.max) {
  old.per <- (x - old.min) / (old.max - old.min)
  new.x <- ((new.max - new.min) * old.per) + new.min
  return(new.x)
}

#-------------------------------------------------------------------------------
# Generalized logistic distribution function
genlog <- function(x, b) {
  dx <- (b[3]/b[2]) * ((exp((b[1] - x) / b[2])) / (1 + exp((b[1] - x) / b[2]))^(b[3] + 1)) 
  return(dx)
}

#-------------------------------------------------------------------------------
# Likelihood function for generalized logistic distribution
genlog.lik <- function(b, x) {
  logl <- sum(log((b[3]/b[2]) * ((exp((b[1] - x) / b[2])) / (1 + exp((b[1] - x) / b[2]))^(b[3] + 1))))
  return(-logl)
}

#-------------------------------------------------------------------------------
# function to calculate TMR and plot from fitted temperature curves
tmr_func <- function(q1, q9, a, b, c, med.temp, plot = T, yl = F) {
  # range of temp values
  xx <- seq(-5, 32, 0.1)
  # parameters of distribution and curve fit
  yy <- genlog(xx, c(a, b, c))
  # scale to 1
  yy <- yy / max(yy)
  # Bd values
  bd <- c(b.curve$a, b.curve$b, b.curve$c)
  bb <- genlog(xx, bd)
  bb <- bb / max(bb)
  # restrict to quantiles
  yy <- yy[xx >= q1 & xx <= q9]
  bb <- bb[xx >= q1 & xx <= q9]
  temp <- xx[xx >= q1 & xx <= q9]
  # calculate TMR as weighted sum of log ratio
  tmr_ratio <- log(bb / yy)
  # fit line
  m1 <- lm(tmr_ratio ~ temp)
  # scale point size by frequency
  sz <- (yy / max(yy)) * 2
  # plot
  if(plot == T) {
    plot(tmr_ratio ~ temp, pch = 19, cex = sz, bty = "l", ylim = c(-2.5, 1.5), col = "grey",
         ylab = "", cex.lab = 1.5, xlab = "")
    if(yl == T) title(ylab = "Thermal mismatch ratio (LTMR)", line = 2, xpd = NA, cex.lab = 1.3)
    tyy <- coef(m1)[1] + coef(m1)[2] * temp
    lines(tyy ~ temp, col = "blue", lwd = 3)
    # median temp
    my <- coef(m1)[1] + coef(m1)[2] * med.temp
    points(my ~ med.temp, pch = 19, cex = 2, col = "blue")
  }

  # mean TMR
  tmr <- weighted.mean(tmr_ratio, w = yy)
  # calculate TMR-temp slope
  tmr_slope <- coef(m1)[2]
  return(c(tmr, tmr_slope))
}


