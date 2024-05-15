####################
# Price shock identification
####################
library(ggplot2)
library(zoo)
library(lmtest)

# Read the data from xlsx file
urea <- readxl::read_xlsx('Urea_price.xlsx', sheet = 'Monthly')
#urea <- urea[36:64,] #for yearly
urea <- urea[529:732,] #for monthly
urea$index <- rownames(urea)
urea$month <- as.Date(as.yearmon(urea$month))

# Using k-fold approach to figure out the best span parameter, basing on the MSE value
# From source httPS_15_15://stats.stackexchange.com/questions/2002/how-do-i-decide-what-span-to-use-in-loess-regression-in-r
k = 5
span.seq <- seq(from = 0.02, to = 0.98, by = 0.02)
cv.error.mtrx <- matrix(rep(x = NA, times = k * length(span.seq)), 
                        nrow = length(span.seq), ncol = k)

set.seed(42)
folds <- sample(x = 1:k, size = length(urea$month), replace = TRUE)

for(i in 1:length(span.seq)) {
  for(j in 1:k) {
    loess.fit <- loess(formula = price ~ as.numeric(month), data = urea[folds != j, ], span = span.seq[i])
    preds <- predict(object = loess.fit, newdata = urea[folds == j, ])
    cv.error.mtrx[i, j] <- mean((urea$price[folds == j] - preds)^2, na.rm = TRUE)
    # some predictions result in `NA` because of the `x` ranges in each fold
  }
}

cv.errors <- rowMeans(cv.error.mtrx)
best.span.i <- which.min(cv.errors)
best.span.i
span.seq[best.span.i]

# Setting time span q
q = .60
#q = span.seq[best.span.i]

# LOESS regression
PS_60 <- loess(price ~ as.numeric(month), data = urea, span = q, degree = 2)

# Autocorrelation detection
mod.err <- lm(PS_60$residuals[2:204] ~ PS_60$residuals[1:203])
bptest(mod.err) # for heteroskedasticity
summary(mod.err)

# Compute Cook's distance
D_60 <- data.frame(cooks.distance(mod.err))
D_60$index <- rownames(D_60) # for merge function

shock_60 <- merge(urea,D_60, all = T, sort = F)
shock_60$shock <- NA
sumofshock <- NULL

for (j in seq(0.95,0.99,0.01)) {
  d = quantile(D_60$cooks.distance.mod.err.,j)
  shock_60$shock <- NA
  for (i in 6:nrow(shock_60)-1) {
    if(shock_60[i,"cooks.distance.mod.err."] > d){
      if(shock_60$price[i] > sum(shock_60$price[i-5],shock_60$price[i-4],shock_60$price[i-3],shock_60$price[i-2],shock_60$price[i-1])/5){
        shock_60$shock[i] = 1
      }
    }
  }
  sumofshock <- rbind(sumofshock, data.frame(CI = j, N = sum(shock_60$shock, na.rm = T), D = d))
}

sumofshock.60 <- sumofshock

resid_60 <- mod.err$model
resid_60$index <- rownames(resid_60)
resid_60 <- merge(shock_60, resid_60)

# Setting cut-off line for Cook's distance
#d = 4/nrow(urea)
#d = 3*mean(D_60$cooks.distance.mod.err.)
d = quantile(D_60$cooks.distance.mod.err.,.98)

# Plot whole figure
p1_60 <- ggplot(data = shock_60, aes(x = month, y = price)) +
  geom_smooth(method = 'loess', span = q, col = 'black', lwd = .5)+
  geom_point(data = shock_60, aes(x = month, y = price), col = 'black', size = .7) +
  theme_bw() + theme(legend.position="bottom")+
  xlab('Year') +
  ylab("Price ($/mt)")+
  ggtitle('span = 0.60')

# Plot error
p2_60 <- ggplot(data = resid_60, aes(x = PS_60$residuals[1:203], y = PS_60$residuals[2:204]))+
  geom_smooth(method = 'lm', col = 'black', lwd = .5)+
  geom_point(data = resid_60, aes(x = PS_60$residuals[1:203], y = PS_60$residuals[2:204]), col = 'black', size = .7)+
  theme_bw() + theme(legend.position="bottom")+
  xlab('Residual at i')+
  ylab('Residual at i+1')
  
# Plot cook's distance
p3_60 <- ggplot(merge(urea,D_60), aes(x = month, y = cooks.distance.mod.err.))+
  geom_line()+
  geom_hline(yintercept = d, lty = 2)+
  theme_bw()+
  xlab('Year')+
  ylab("Cook's distance")

# Combine figures together to a jpeg form
jpeg('shock identification.jpeg', width = 7.5, height = 5, units = 'in', res = 200)
ggpubr::ggarrange(p1_15,p1_60,p2_15,p2_60,p3_15,p3_60, ncol = 2, nrow = 3, common.legend = T, labels = 'auto', vjust = c(1.5,1.5,-.5,-.5,-.5,-.5), font.label = list(size = 10))
dev.off()
