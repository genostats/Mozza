
triangle.plot <- function(a, b, c) {
  S <- a+b+c
  list(x = (b + 0.5*c)/S, y = sqrt(3)/2*c/S)
}

triangle.tiling <- function(n) {
  a <-  unlist( sapply( seq(n, 1), function(x) rep( (n-x)/(n-1), x) ))
  b <- unlist( sapply( seq(0, 1, by = 1/(n-1)), function(x) seq(0, 1-x, by = 1/(n-1)) ) )
  c <- 1 - a - b
  data.frame(a,b,c) 
}

# plot( do.call(triangle.plot, tiling(21)), pch = 16 )

