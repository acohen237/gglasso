# For testing the lasso, group lasso, sparse group lasso
# In other words, the gglasso DA project

# Set the n and p values
set.seed(1010)
n <- 1000
p <- 100

# There are different ways to generate a data matrix, I'm gonna be a dummy

X <- matrix(data = rnorm(n*p, mean=0, sd=1), nrow = n, ncol = p)

eps <- rnorm(n, mean = 0, sd=1)

beta_star <- c(rep(5,5), c(5,-5,2,0,0), rep(-5,5), c(2,-3,8,0,0), rep(0,(p-20)))


# we generate the data so that there are groups of 5 vars each
grp <- rep(1:20, each=5)

y <- X%*%beta_star + eps

# Now we need to try the lassos
out <- gglasso(X,y,group=grp, loss='ls')

out$beta[,50]