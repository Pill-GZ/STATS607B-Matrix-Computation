m <- 100
n <- 15
i <- 1:m
j <- 1:n
X <- outer(i, j, function(i,j){((i-1)/(m-1))^(j-1)})
Y = exp(sin(4*(i-1)/(m-1)))/2006.787453080206

lm(Y~X)

# QRstdGS
Q <- matrix(0, m, n)
R <- matrix(0, n, n)
for (j in 1:n) {
  #cat("j = ", j)
  v <- X[,j]
  if (j>1) {
    for (i in 1:(j-1)) {
      #cat("i = ",i)
      R[i,j] <- sum(X[,j]*Q[,i])
      v <- v - R[i,j]*Q[,i]
    }
  }
  R[j,j] <- sqrt(sum(v^2))
  Q[,j] <- v/R[j,j]
}
max(abs(Q%*%R-X))
# solve Rb=Q'Y=y
b <- numeric(n)
y <- t(Q)%*%Y
b[n] <- y[n]/R[n,n]
for (j in (n-1):1) {
  b[j] <- (y[j] - sum(b[(j+1):n]*R[j,(j+1):n]))/R[j,j]
}
(b15QRstdGE <- b[15])

# QRmodGS
Q <- matrix(0, m, n)
R <- matrix(0, n, n)
for (j in 1:n) {
  #cat("j = ", j)
  v <- X[,j]
  if (j>1) {
    for (i in 1:(j-1)) {
      #cat("i = ",i)
      R[i,j] <- sum(v*Q[,i])
      v <- v - R[i,j]*Q[,i]
    }
  }
  R[j,j] <- sqrt(sum(v^2))
  Q[,j] <- v/R[j,j]
}
max(abs(Q%*%R-X))
# solve Rb=Q'Y=y
b <- numeric(n)
y <- t(Q)%*%Y
b[n] <- y[n]/R[n,n]
for (j in (n-1):1) {
  b[j] <- (y[j] - sum(b[(j+1):n]*R[j,(j+1):n]))/R[j,j]
}
(b15QRmodGE <- b[15])

# QRHH
#Q <- diag(m)
R <- X
y <- Y
for (j in 1:n) {
  x <- R[j:m,j]
  v <- x + sign(x[1]) * sqrt(sum(x^2)) * c(1,rep(0,m-j))
  R[j:m,j:n] <- (diag(m-j+1) - 2*outer(v,v)/sum(v^2)) %*% R[j:m,j:n]
  #  Q[j:m,j:n] <- (diag(m-j+1) - 2*outer(v,v)/sum(v^2)) %*% Q[j:m,j:n]
  y[j:m] <- (diag(m-j+1) - 2*outer(v,v)/sum(v^2)) %*% y[j:m]
}
#max(Q%*%R - X)
# solve Rb=Q'Y=y
b[n] <- y[n]/R[n,n]
for (j in (n-1):1) {
  b[j] <- (y[j] - sum(b[(j+1):n]*R[j,(j+1):n]))/R[j,j]
}
(b15QRHouseHolder <- b[15])

# Chol and solve normal equation (X'X)b=X'Y=w
XX <- t(X)%*%X
w <- t(X)%*%Y
L <- matrix(0,n,n)
for (j in 1:n) {
  L[j,j] <- sqrt(XX[j,j] - sum(L[j,1:(j-1)]^2))
  for (i in (j+1):n) {
    L[i,j] <- (XX[i,j] - sum(L[j,1:(j-1)]*L[i,1:(j-1)])) / L[j,j]
  }
}
L
# the Cholesky decomposition breaks down due to rounding in R!
# although LL' is close to X'X 
max(abs(L%*%t(L)-XX),na.rm = T)

# perhaps I should do it in MATLAB

## Question 2



