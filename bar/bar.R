dyn.load("bar.so")

bar <- function(x) {
   if (!is.numeric(x))
      stop("argument x must be numeric")
   out <- .C("bar",
             n=as.integer(length(x)),
             x=as.double(x))
   return(out$x)
}

x <- rnorm(10)

x2 <- bar(x)

identical(x^2, x2)