# ---------------- SCORE FUNCTION ----------------
score <- function(theta, x) {
alpha <- theta[1]
lambda <- theta[2]
  
if(alpha <= 0 || lambda <= 0) return(c(NA, NA))
  
n <- length(x)
A <- 1 - exp(-lambda * x)
  
A[A <= 1e-12] <- 1e-12
  
d_alpha <- n/alpha + sum(log(A))
d_lambda <- n/lambda - sum(x) +(alpha - 1)*sum(x * exp(-lambda*x)/A)

return(c(d_alpha, d_lambda))
}

# ---------------- JACOBIAN ----------------------
jacobian <- function(theta, x) {
alpha <- theta[1]
lambda <- theta[2]
  
if(alpha <= 0 || lambda <= 0)
return(matrix(NA,2,2))
  
n <- length(x)
A <- 1 - exp(-lambda*x)
  
A[A <= 1e-12] <- 1e-12
  
d11 <- -n/(alpha^2)
d12 <- sum(x * exp(-lambda*x)/A)
  
d22 <- -n/(lambda^2) - (alpha - 1)*sum(x^2 * exp(-lambda*x)/(A^2))

return(matrix(c(d11, d12,
                  d12, d22), 2, 2))
}

# ---------------- NEWTON-RAPHSON ----------------
NR <- function(x, init, tol=1e-6, max_iter=100) {
  
theta <- init_guess(x)
  
for(i in 1:max_iter) {
    
S <- score(theta, x)
J <- jacobian(theta, x)
    
if(any(is.na(S)) || any(is.na(J)))
return(c(NA, NA))
    
step <- tryCatch(solve(J, S), error=function(e) c(NA, NA))
    
if(any(is.na(step))) return(c(NA, NA))
    
theta_new <- theta - step
    
# enforce positivity
if(any(theta_new <= 0)) return(c(NA, NA))
    
# stopping rule
if(max(abs(theta_new - theta)) < tol) {
cat("Converged in", i, "iterations\n")
return(theta_new)
 }
    
theta <- theta_new
 }
  
return(c(NA, NA))
}

# ---------------- INITIAL VALUE ----------------
lambda_init <- 1 / mean(data)

A <- 1 - exp(-lambda_init * data)
alpha_init <- mean(A) / (1 - mean(A))

init <- c(alpha_init, lambda_init)
init

# ---------------- MLE --------------------------
mle <- NR(data, init)

alpha_hat <- mle[1]
lambda_hat <- mle[2]

cat("MLE alpha =", alpha_hat, "\n")
cat("MLE lambda =", lambda_hat, "\n")

# ---------------- INVERSE TRANSFORM ----------------
r_sample <- function(n, alpha, lambda) {
u <- runif(n)
x <- - (1/lambda) * log(1 - u^(1/alpha))
return(x)
}

# ---------------- PARAMETRIC BOOTSTRAP ----------
B <- 500

alpha_boot <- c()
lambda_boot <- c()

for(b in 1:B) {
x_star <- r_sample(n, alpha_hat, lambda_hat)
est <- NR(x_star, init)
  
if(!any(is.na(est))) {
alpha_boot <- c(alpha_boot, est[1])
lambda_boot <- c(lambda_boot, est[2])
  }
}

cat("Parametric CI alpha:",
    quantile(alpha_boot, c(0.025, 0.975)), "\n")

cat("Parametric CI lambda:",
    quantile(lambda_boot, c(0.025, 0.975)), "\n")

# ---------------- NON-PARAMETRIC BOOTSTRAP -------
alpha_boot_np <- c()
lambda_boot_np <- c()

for(b in 1:B) {
x_star <- sample(data, n, replace=TRUE)
est <- NR(x_star, init)
  
if(!any(is.na(est))) {
alpha_boot_np <- c(alpha_boot_np, est[1])
lambda_boot_np <- c(lambda_boot_np, est[2])
  }
}

cat("Non-parametric CI alpha:",
    quantile(alpha_boot_np, c(0.025, 0.975)), "\n")

cat("Non-parametric CI lambda:",
    quantile(lambda_boot_np, c(0.025, 0.975)), "\n")

