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

beta50 <- out$beta[,50]

# Next we try ls_sparse
out_sp <- gglasso(X,y,group=grp, loss = 'ls_sparse')

beta50_sparse <- out$beta[,50]


###############################################################################
## A much smaller data set for debugging


# Set the n and p values
set.seed(1010)
n <- 100
p <- 9

# There are different ways to generate a data matrix, I'm gonna be a dummy

X <- matrix(data = rnorm(n*p, mean=0, sd=1), nrow = n, ncol = p)

eps <- rnorm(n, mean = 0, sd=1)

beta_star <- c(rep(5,3), c(-5,0,-3),rep(0,3))


# we generate the data so that there are groups of 5 vars each
grp <- rep(1:3, each=3)

y <- X%*%beta_star + eps

# Now we need to try the lassos
out <- gglasso(X,y,group=grp, loss='ls')

#beta50 <- out$beta[,50]

# Next we try ls_sparse
out_sp <- gglasso(X,y,group=grp, loss = 'ls_sparse')

#beta50_sparse <- out$beta[,50]

group_norm <- function(x, grp) sum(by(x, grp, function(x) sqrt(sum(x^2))))
sp_group_norm <- function(x, grp, alp=.05) group_norm(x,grp)*(1-alp) + alp*sum(abs(x))




# Dan's code
par(mfrow=c(2,2))
plot(out)
plot(out_sp)
b1 = apply(out$beta, 2, sd)
b2 = apply(out_sp$beta, 2, sd)
matplot(b1, t(out$beta), ty='l', lty=1, col = grp)
matplot(b2, t(out_sp$beta), ty='l', lty=1, col = grp)




###########################################################################
# Dan's code, redux (using group_norm instead of sd)
group_norm <- function(x) sum(by(x, grp, function(x) sqrt(sum(x^2))))
sp_group_norm <- function(x, alp=.05) group_norm(x)*(1-alp) + alp*sum(abs(x))



par(mfrow=c(1,2))
b1 = apply(out$beta, 2, group_norm)
b2 = apply(out_sp$beta, 2, sp_group_norm)
matplot(b1, t(out$beta), ty='l', lty=1, col = grp)
matplot(b2, t(out_sp$beta), ty='l', lty=1, col = grp)


###########################################################################
###############################################################################
###########################################################################
# OK, let's play again with some data, for both the n<p and p<n cases
# Recall Dan said there are 4 cases: (sparsity/no sparcity within groups) X (p < n vs p > n)

set.seed(1010)
n <- 10000
p <- 1000



X <- matrix(data = rnorm(n*p, mean=0, sd=1), nrow = n, ncol = p)

eps <- rnorm(n, mean = 0, sd=1)


beta_star <- rep(0,p)
for (i in 0:(20-1)) {
  beta_star[(50*i + 1):(50*(i+1))] <- rep((-1)^i*i, 50)
  beta_star[(50*i + 20):(50*i + 29)] <- rep(0,10)
}


grp <- rep(1:20, each=50)

y <- X%*%beta_star + eps

# Now we need to try the lassos
out <- gglasso(X,y,group=grp, loss='ls')

# Next we try ls_sparse
out_sp <- gglasso(X,y,group=grp, loss = 'ls_sparse')


group_norm <- function(x) sum(by(x, grp, function(x) sqrt(sum(x^2))))
sp_group_norm <- function(x, alp=.05) group_norm(x)*(1-alp) + alp*sum(abs(x))


matplot(b1, t(out$beta), ty='l', lty=1, col = grp)
matplot(b2, t(out_sp$beta), ty='l', lty=1, col = grp)

