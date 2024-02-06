dfSim <- function (c0 = 1, c2 = 1, n = 400, scale = 2, verbose = TRUE, add = FALSE){
	if (add  == FALSE) {
	        if (verbose) 
            cat("2 term additive")
        x0 <- runif(n, 0, 1)
        #x1 <- runif(n, 0, 1)
        x2 <- runif(n, 0, 1)
        x3 <- runif(n, 0, 1)
        f0 <- function(x) x + 10 * sin(pi * c0 * x)
        #f1 <- function(x) exp(2 * x)
        f2 <- function(x) 0.2 * x^11 * (10 * (1 - x))^6 + 10 * 
            (10 * x)^3 * (1 - x)^10
        #f3 <- function(x) 0 * x
        f <- f0(x0) + c2*f2(x2)
		#f <- c2*f2(x2)
		
		e <- rnorm(n, 0, scale)
        y <- f + e
        
        data <- data.frame(y = y, x0 = x0, x2 = x2, x3 = x3,
            f = f, f0 = f0(x0))
        return(data)
    } else {
        if (verbose) 
            cat("2 term additive + random effect")
			        x0 <- runif(n, 0, 1)
        #x1 <- runif(n, 0, 1)
        x2 <- runif(n, 0, 1)
        x3 <- runif(n, -10, 10)
        f0 <- function(x) x + 2 * sin(pi * c0 * x)
        #f1 <- function(x) exp(2 * x)
        f2 <- function(x) 0.2 * x^11 * (10 * (1 - x))^6 + 10 * 
            (10 * x)^3 * (1 - x)^10
        #f3 <- function(x) 0 * x
        f <- f0(x0) + c2*f2(x2)
		#f <- c2 * f2(x2)
		
		e <- rnorm(n, 0, scale)
        y <- f + e
        
        fac <- rep(1:10, n/10)
        data <- data.frame(y = y, x0 = x0, x2 = x2, x3 = x3, f = f, f0 = f0(x0), fac = as.factor(fac))
        return(data)
    }
}


likeli <- function(b) {
  if (class(b)[1]=="lm") {
    return(sum(dnorm(b$model$y, mean = b$fitted.values, sd = sigma(b), log = TRUE)))
	} else {
    return(sum(dnorm(b$y, mean = b$fitted.values, sd = sqrt(b$sig2), log = TRUE)))
    }
}

source("fcts.r")

# Whether a server running a unix system is used. If not automatically a windows system is assumed
server <- TRUE
	
# How many cores in the parallized part are used.
NumberOfCores <- 100

if(server == TRUE){
require("doMC")
require("foreach")
registerDoMC(cores = NumberOfCores)
}

require(mgcv)
require(splines)

set.seed(1234)

n      <- 20  
size1  <- c(0, .1, .2, .4, .6, .8, 1, 1.2)#seq(0,2,by=.2)
size2  <- 0
sig2   <- .2
NoSets <- 1000
crit   <- "GCV.Cp"

nmrr0 <- nmrrdf <- biasinc <- numeric(length(size1))

selectRates0 <- selectRatesdf <- matrix(0, nrow = 3, ncol = length(size1)) 

postSelectBias <- freqAICdf <- freqAIC0 <- postSelectBiasAIC0 <- postSelectBiasAICdf <- numeric(NoSets)

for(m in 1:length(size1)){
for(j in 1:NoSets){

#datas <- dfSim(c0 = size1[m], c2 = size2, n = n, scale = sig2, verbose = FALSE)
datas <- dfSim(c0 = size1[m], c2 = size2, n = n, scale = sig2, verbose = FALSE)

#estimate over-simple complex, over-complex
	m1 <- lm(y ~ x0, data = datas)
	m2 <- gam(y ~ s(x0, bs="ps"), data = datas, method = crit)
	#m3 <- gam(y ~ s(x2, k= 25), data = datas, method = crit)
	m3 <- gam(y ~ s(x0, bs="ps") + s(x3, bs="ps", k = 11), data = datas, method = crit)
	#m3 <- gam(y ~ s(x0, bs="ps") + s(fac, bs="re"), data = datas, method = crit)
	
#simulate evaluation data
evalSize <- 1000
#evalData <- dfSim(c0 = size1[m], c2 = size2, n = evalSize, scale = sig2, verbose = FALSE)
evalData <- dfSim(c0 = size1[m], c2 = size2, n = evalSize, scale = sig2, verbose = FALSE)

#calculate post selection bias
biases <- c(
sum((predict(m1, evalData)-evalData$y)^2)/evalSize,
sum((predict(m2, evalData)-evalData$y)^2)/evalSize,
sum((predict(m3, evalData)-evalData$y)^2)/evalSize
)

freqAIC0[j] <- which.min(c(
-likeli(m1) + sum(!is.na(coef(m1))),
-likeli(m2) + sum(m2$edf),
-likeli(m3) + sum(m3$edf)
))

postSelectBiasAIC0[j] <- biases[freqAIC0[j]]

freqAICdf[j] <- which.min(c(
-likeli(m1) + sum(!is.na(coef(m1))),
-likeli(m2) + sum(m2$edf) + sum(cAIC(m2)),
-likeli(m3) + sum(m3$edf) + sum(cAIC(m3))
))

postSelectBiasAICdf[j] <- biases[freqAICdf[j]]

postSelectBias[j] <- postSelectBiasAIC0[j] - postSelectBiasAICdf[j]

}

selectRates0[1, m] <- sum(freqAIC0==1)/length(freqAIC0)
selectRates0[2, m] <- sum(freqAIC0==2)/length(freqAIC0)
selectRates0[3, m] <- sum(freqAIC0==3)/length(freqAIC0)

selectRatesdf[1, m] <- sum(freqAICdf==1)/length(freqAICdf)
selectRatesdf[2, m] <- sum(freqAICdf==2)/length(freqAICdf)
selectRatesdf[3, m] <- sum(freqAICdf==3)/length(freqAICdf)

nmrr0[m] <- 1- sum(freqAIC0==1)/length(freqAIC0)

nmrrdf[m] <- 1- sum(freqAICdf==1)/length(freqAICdf)

biasinc[m] <- round(mean(postSelectBiasAIC0)/mean(postSelectBiasAICdf)*100,2)

cat("for k = ", m, "null model rejection rate AIC0 ", nmrr0[m], "\n")

cat("null model rejection rate AICdf ", nmrrdf[m], "\n")

cat("bias = ", sum(postSelectBias), "percentage = ", biasinc[m], "\n")

#boxplot(postSelectBiasAIC0, postSelectBiasAICdf, main=paste("s = ", s, ", l = ", m))
}


#par(mfrow = c(3, 1))
#plot(size1, biasinc-100, type ="l", ylim=c(0,5), ylab = "bias increase in %", xlab = "effect size")
## Create the input vectors.
#colors = c("lightgreen","plum1","royalblue")
#models <- c("m1","m2","m3")
## Create the bar chart
#barplot(selectRatesdf, main = "Selection rates AICdf", names.arg = size1, xlab = "effect size", ylab = "rate", col = colors)
## Add the legend to the chart
#legend("topleft", models, cex = 1, fill = colors)
#barplot(selectRates0, main = "Selection rates AIC0", names.arg = size1, xlab = "effect size", ylab = "rate", col = colors)
## Add the legend to the chart
#legend("topleft", models, cex = 1, fill = colors)

filename <- paste("sin_1234_REML_bar.eps")

#setEPS()
#postscript(file = filename, width = 16, height = 5)
pdf(file = filename, width = 16, height = 5)

par(mfrow = c(1, 3), cex=1.2)
plot(size1, biasinc-100, type ="l", ylim=c(0,6), ylab = "bias increase in %", xlab = "effect size")
# Create the input vectors.
colors = c("white","gray","black")
models <- c("m1","m2","m3")
# Create the bar chart
barplot(selectRatesdf, main = expression(paste("Including df(", lambda, ")")), names.arg = size1, xlab = "effect size", ylab = "rate", col = colors)
# Add the legend to the chart
#legend("topleft", models, cex = 1, fill = colors)
barplot(selectRates0, main = expression(paste("Ingnoring df(", lambda, ")")), names.arg = size1, xlab = "effect size", ylab = "rate", col = colors)

dev.off()
# Add the legend to the chart
#legend("topleft", models, cex = 1, fill = colors)
#plot the effect sizes
#f0 <- function(x,c0) x + 2 * sin(pi * c0 * x)
#x<-seq(0,1,by=.01)
#par(mfrow = c(3, 3))
#for(i in size1){
#plot(x,f0(x,i), type ="l")
#}
#
#f2 <- function(x,c0) x + 0.2 * c0 * (10 * x)^3 * (1 - x)^10 + 2 * sin(pi * c0 * x)
#x<-seq(0,1,by=.01)
#plot(x,f2(x,2.5), type ="l")
##boxplot(postSelectBiasAIC0, postSelectBiasAICdf, main=paste("mean0 = ", round(mean(postSelectBiasAIC0),2), "meandf = ", round(mean(postSelectBiasAICdf), 2)))