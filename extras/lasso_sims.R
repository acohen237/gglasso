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

y <- X%*%beta_star + eps

groups <- rep(1:20, each=5)

# Now we need to try the lassos

out <- gglasso(X, y, group = groups, loss = 'ls')

out_sp <- gglasso(X,y, group = groups, loss="ls_sparse")


grp_norm <- function(x, grp, pf = sqrt(table(grp))){
  out = by(x,grp,function(x) sqrt(sum(x^2)))
  sum(out*pf)
}

sparse_grp_norm <- function(x, grp, alpha = .05, pf = sqrt(table(grp))){
  out1 = grp_norm(x,grp,pf)
  (1-alpha)*out1 + alpha*sum(abs(x))
}

matplot(apply(out$beta, 2, grp_norm, grp=groups), t(out$beta), ty='l',lty=1,col=2)
matlines(apply(out_sp$beta, 2, sparse_grp_norm, grp=groups), t(out_sp$beta), lty=1,col=4)
