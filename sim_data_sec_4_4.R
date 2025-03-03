library(data.table)
library(copula)

expit <- function(x) {
  1 / (1 + exp(-x))
}

sim_func <- function(n=500,beta = -0.3){
  
  # Baseline ----------------------------------------------------------------
  
  # obesity
  B <- rbinom(n,1,0.1)
  # age
  C <- runif(n,25,35)
  # L_0
  
  w0 <- runif(n)
  L_0 <- qnorm(w0, mean  = 11-0.05*B -0.02*C, sd = 0.5)
  A_0 <- rbinom(n,1,0.5)
  
  # Later time points -------------------------------------------------------
  w1 <- runif(n)
  L_1 <- qnorm(w1, mean  = L_0 + 0.5*A_0, sd = 0.1)
  A_1 <- A_0
  
  w2 <- runif(n)
  L_2 <- qnorm(w2, mean  = L_1 + 0.5*A_1, sd = 0.1)
  A_2 <- A_0
  w3 <- runif(n)
  L_3 <- qnorm(w3, mean  = L_2 + 0.5*A_2, sd = 0.1)
  A_3 <- rbinom(n,1,expit(20 - 2*L_3 + 0.05*B + 0.01 *C))
  
  w4 <- runif(n)
  L_4 <- qnorm(w4, mean  = L_3 + 0.5*A_3, sd = 0.1)
  A_4 <- rbinom(n,1,expit(20 - 2*L_4  + 0.05*B + 0.01 *C))
  
  
  # Simulate Y --------------------------------------------------------------
  
  cop <- tCopula(param = -0.5, df = 5, dim = 2, dispstr = "un")
  
  wy <- runif(n)
  
  nu4 <- cCopula(as.matrix(cbind(w4,wy)),copula = cop,inverse = T)[,2]
  nu3 <- cCopula(as.matrix(cbind(w3,nu4)),copula = cop,inverse = T)[,2]
  nu2 <- cCopula(as.matrix(cbind(w2,nu3)),copula = cop,inverse = T)[,2]
  nu1 <- cCopula(as.matrix(cbind(w1,nu2)),copula = cop,inverse = T)[,2]
  nu0 <- cCopula(as.matrix(cbind(w0,nu1)),copula = cop,inverse = T)[,2]
  
  Y  <- qbinom(nu0,size = 1, prob = expit(-2 + 0.1*B + 0.02*C + beta *(A_0 + A_1 +A_2 + A_3 + A_4)))
  
  
  dat <- data.table(B = B, C= C,
                    A_0 = A_0, L_0 = L_0, Y_0 = 0,
                    A_1 = A_1, L_1 = L_1, Y_1 = 0,
                    A_2 = A_2, L_2 = L_2, Y_2 = 0,
                    A_3 = A_3, L_3 = L_3, Y_3 = 0,
                    A_4 = A_4, L_4 = L_4, Y_4 = Y)
  
  
  return(dat)
}

sim_func()
