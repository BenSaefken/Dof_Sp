
likeli <- function(b) {
  if (class(b)[1]=="lm") {
    return(sum(dnorm(b$model$y, mean = b$fitted.values, sd = sigma(b), log = TRUE)))
	} else {
    return(sum(dnorm(b$y, mean = b$fitted.values, sd = sqrt(b$sig2), log = TRUE)))
    }
}

source("fcts.r")

require(mgcv)
require(splines)

set.seed(101)


sambiaData <- read.table(paste("sambia.raw", sep = ""), header = TRUE)
names(sambiaData) <- c("zscore", "csex", "cfeed", "cage", "mage", "mheight", "mbmi", "medu", "mwork", "district", "region", "time")
#sambiaData$zscore <- (sambiaData$zscore-mean(sambiaData$zscore))/sqrt(var(sambiaData$zscore))
sambiaData$district <- as.factor(sambiaData$district)

for (i in unique(sambiaData$cage)){
  correctAge <- sambiaData$cage == i
  zscores <- sambiaData[correctAge,]$zscore
  sambiaData$zscore[correctAge] <- (zscores-median(zscores))/sqrt(var(zscores))
}

#sambiaData <- sambiaData[1:100,]

modelFormulas <- list()


modelFormulas[[1]] <- zscore ~ csex + medu + mwork + s(cfeed) + s(cage) + s(mage) + s(mheight) + s(mbmi) + s(district, bs = "re")#
 
modelFormulas[[2]] <- zscore ~ csex + medu + mwork + s(cfeed) + s(cage) + s(mage) + s(mheight) + mbmi + s(district, bs = "re")#
modelFormulas[[3]] <- zscore ~ csex + medu + mwork + s(cfeed) + s(cage) + s(mage) + mheight + s(mbmi) + s(district, bs = "re")#
modelFormulas[[4]] <- zscore ~ csex + medu + mwork + s(cfeed) + s(cage) + mage + s(mheight) + s(mbmi) + s(district, bs = "re")#
modelFormulas[[5]] <- zscore ~ csex + medu + mwork + s(cfeed) + cage + s(mage) + s(mheight) + s(mbmi) + s(district, bs = "re")
modelFormulas[[6]] <- zscore ~ csex + medu + mwork + cfeed + s(cage) + s(mage) + s(mheight) + s(mbmi) + s(district, bs = "re")

modelFormulas[[7]] <- zscore ~ csex + medu + mwork + s(cfeed) + s(cage) + s(mage) + mheight + mbmi + s(district, bs = "re")#
modelFormulas[[8]] <- zscore ~ csex + medu + mwork + s(cfeed) + s(cage) + mage + s(mheight) + mbmi + s(district, bs = "re")#
modelFormulas[[9]] <- zscore ~ csex + medu + mwork + s(cfeed) + cage + s(mage) + s(mheight) + mbmi + s(district, bs = "re")
modelFormulas[[10]] <- zscore ~ csex + medu + mwork + cfeed + s(cage) + s(mage) + s(mheight) + mbmi + s(district, bs = "re")

modelFormulas[[11]] <- zscore ~ csex + medu + mwork + s(cfeed) + s(cage) + mage + mheight + s(mbmi) + s(district, bs = "re")#
modelFormulas[[12]] <- zscore ~ csex + medu + mwork + s(cfeed) + cage + s(mage) + mheight + s(mbmi) + s(district, bs = "re")
modelFormulas[[13]] <- zscore ~ csex + medu + mwork + cfeed + s(cage) + s(mage) + mheight + s(mbmi) + s(district, bs = "re")
  
modelFormulas[[14]] <- zscore ~ csex + medu + mwork + s(cfeed) + cage + mage + s(mheight) + s(mbmi) + s(district, bs = "re")
modelFormulas[[15]] <- zscore ~ csex + medu + mwork + cfeed + s(cage) + mage + s(mheight) + s(mbmi) + s(district, bs = "re") 

modelFormulas[[16]] <- zscore ~ csex + medu + mwork + cfeed + cage + s(mage) + s(mheight) + s(mbmi) + s(district, bs = "re") 

modelFormulas[[17]] <- zscore ~ csex + medu + mwork + s(cfeed) + s(cage) + mage + mheight + mbmi + s(district, bs = "re")#
modelFormulas[[18]] <- zscore ~ csex + medu + mwork + s(cfeed) + cage + s(mage) + mheight + mbmi + s(district, bs = "re")
modelFormulas[[19]] <- zscore ~ csex + medu + mwork + cfeed + s(cage) + s(mage) + mheight + mbmi + s(district, bs = "re") 

modelFormulas[[20]] <- zscore ~ csex + medu + mwork + s(cfeed) + cage + mage + s(mheight) + mbmi + s(district, bs = "re")
modelFormulas[[21]] <- zscore ~ csex + medu + mwork + cfeed + s(cage) + mage + s(mheight) + mbmi + s(district, bs = "re") 

modelFormulas[[22]] <- zscore ~ csex + medu + mwork + cfeed + cage + s(mage) + s(mheight) + mbmi + s(district, bs = "re") 

modelFormulas[[23]] <- zscore ~ csex + medu + mwork + s(cfeed) + cage + mage + mheight + s(mbmi) + s(district, bs = "re")
modelFormulas[[24]] <- zscore ~ csex + medu + mwork + cfeed + s(cage) + mage + mheight + s(mbmi) + s(district, bs = "re") 

modelFormulas[[25]] <- zscore ~ csex + medu + mwork + cfeed + cage + s(mage) + mheight + s(mbmi) + s(district, bs = "re") 

modelFormulas[[26]] <- zscore ~ csex + medu + mwork + cfeed + cage + mage + s(mheight) + s(mbmi) + s(district, bs = "re") 

modelFormulas[[27]] <- zscore ~ csex + medu + mwork + s(cfeed) + cage + mage + mheight + mbmi + s(district, bs = "re")
modelFormulas[[28]] <- zscore ~ csex + medu + mwork + cfeed + s(cage) + mage + mheight + mbmi + s(district, bs = "re")
modelFormulas[[29]] <- zscore ~ csex + medu + mwork + cfeed + cage + s(mage) + mheight + mbmi + s(district, bs = "re")
modelFormulas[[30]] <- zscore ~ csex + medu + mwork + cfeed + cage + mage + s(mheight) + mbmi + s(district, bs = "re")
modelFormulas[[31]] <- zscore ~ csex + medu + mwork + cfeed + cage + mage + mheight + s(mbmi) + s(district, bs = "re")

modelFormulas[[32]] <- zscore ~ csex + medu + mwork + cfeed + cage + mage + mheight + mbmi + s(district, bs = "re")#

df0 <- dfl <- numeric(length(modelFormulas))

crit <- "GCV.Cp"

crit <- "REML"

for(k in 1:length(modelFormulas)){
	model <- gam(modelFormulas[[k]], data = sambiaData, method = crit)
	df0[k] <-  -likeli(model) + sum(model$edf)
	cat("df0_", k, "\n")
	dfl[k] <- -likeli(model) + sum(model$edf) + sum(cAIC(model))
	cat("dfl_", k, "\n")

}

#tableString <- "results/mixedModelFull_missing.txt"
               
results <- round(data.frame(df0, dfl), 2)

which(df0 == min(df0, na.rm = TRUE))
which(dfl == min(dfl, na.rm = TRUE))

which(results[,1] == min(results[,1], na.rm = TRUE))
which(results[,2] == min(results[,2], na.rm = TRUE))

resultsREML[c(1,2,3,4,7,8,11,17,32),]

model <- gam(modelFormulas[[1]], data = sambiaData, method = crit)

par(mfrow = c(2, 3), cex = 1.2) 
plot(model, select = 1, main = "(a)", shade = TRUE)
plot(model, select = 2, main = "(b)", shade = TRUE)
plot(model, select = 3, main = "(c)", shade = TRUE)
plot(model, select = 4, main = "(d)", shade = TRUE)
plot(model, select = 5, main = "(e)", shade = TRUE)
plot(model, select = 6, main = "(f)", shade = TRUE)


#names(results) <- c("cv", "boot", "cs", "error", "cvError", "bootError", "csError")
              
#write.table(results, file = tableString, quote = FALSE, row.names = FALSE)