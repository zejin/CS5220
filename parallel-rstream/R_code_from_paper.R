# R_code_from_paper.R
# Using RngStreams for Parallel Random Number Generation in C++ and R
# Computational Statistics Andrew Karl, Randy Eubank, Jelena Milovanovic, Mark
# Reiser, Dennis Young akarl@asu.edu
x <- NULL
for (i in 1:200) {
  set.seed(123)
  .Random.seed[3:626] <- as.integer(.Random.seed[3:626] + 
                                    (i - 1) * 100)
  zTemp <- rnorm(1000)
  x <- c(x, sqrt(1000) * mean(zTemp))
}
set.seed(123)
y <- NULL
for (i in 1:200) {
  zTemp <- rnorm(1000)
  y <- c(y, sqrt(1000) * mean(zTemp))
}
cat(shapiro.test(x)$p.value, shapiro.test(y)$p.value, "\n")
#################################### 
RNGkind("L'Ecuyer-CMRG")
set.seed(123)
cat(.Random.seed[2:7], "\n")
.Random.seed[3] + 2^{32}
runif(1)
#################################### 
library(rstream)
rngList <- c(new("rstream.mrg32k3a", seed = c(1806547166, 
                 3311292359, 643431772, 1162448557, 3335719306, 
                 4161054083), force.seed = TRUE), 
             replicate(3, new("rstream.mrg32k3a")))
sapply(rngList, rstream.sample)
#################################### 
library(parallel)
cl <- makeCluster(4)
clusterSetRNGStream(cl, 123)
parSapply(cl, rep(1, 4), runif)
stopCluster(cl)
#################################### 
library(parallel)
RNGkind("L'Ecuyer-CMRG")
set.seed(123)
tryCatch(
    {
	unlist(mclapply(rep(1, 3), runif), mc.cores = 3)
	# the above line will generate an error on Windows machines
    },
    error = function(e) {
    	print(e)
    }
)
runif(1)
#################################### 
library(rstream)
library(parallel)
rngList <- c(new("rstream.mrg32k3a", 
                 seed = c(1806547166, 3311292359, 643431772,
                          1162448557, 3335719306, 4161054083), 
                 force.seed = TRUE), 
             replicate(3, new("rstream.mrg32k3a")))
unlist(mclapply(rngList, rstream.sample, mc.cores = 4))
#################################### 