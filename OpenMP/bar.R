dyn.load('bar.so')

a <- rnorm(10)
b <- rnorm(10)

.Call('dotProd', a, b)