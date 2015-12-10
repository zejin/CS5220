dyn.load('bar.so')

a <- rnorm(8)

.Call('noise', a)