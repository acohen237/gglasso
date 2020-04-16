# short script to play around with the SGL package (not included in this package)
# and to compare with our package. Ultimate goal is to confirm our routines--in particular the ls_sparse subroutine-- are working.

# To start, let's look at the example given in the SGL package documentation
# This is from: https://cran.r-project.org/web/packages/SGL/SGL.pdf

library(SGL)

n = 50; p = 100; size.groups = 10
index <- ceiling(1:p / size.groups)
X = matrix(rnorm(n*p), ncol = p, nrow = n)
beta = (-2:2)
y = X[,1:5] %*% beta + 0.1*rnorm(n)
data = list(x=X, y=y)
fit = SGL(data, index, type = "linear")

# Examine the beta matrix; it seems like this algorithm is NOT recovering the proper coefficients...
