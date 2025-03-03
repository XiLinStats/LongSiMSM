library(data.table)
library(copula)

expit <- function(x) {
  1 / (1 + exp(-x))
}

sim_func <- function(n = 500, 
                     rho1 = -0.1, rho2 = 0.2,
                     K = 5, betaS = -0.5,
                     betaB2 = -0.5, betaSB2 = 0.3){
  # copulas
  copL <- tCopula(param = rho1, df = 5, dim = 2, dispstr = "un")
  copW <- normalCopula(param = rho2, dim = 2, dispstr = "un")
  dat <- data.table()
  dat$surv <- rep(1,n)
  
  # baseline vars
  dat[, B1 := rexp(n, rate = 2)]
  dat[, B2 := rbinom(n, size = 1, prob = 0.7)]
  baseline_vars <- c("B1","B2")
  
  # a data table to store calculated quantiles to be carried over to the next time step
  comp_dat <- copy(dat[,baseline_vars, with = F])
  
  for (k in 0:K){
    n_surv <- sum(dat$surv, na.rm = T)
    if (k>0){
      comp_dat <- comp_dat[Y == 0,]
      if (k > 1){
        comp_dat[,grep("_prev$", colnames(comp_dat), value = TRUE):= NULL]
      }
      # reset names
      setnames(comp_dat,
               old = setdiff(names(comp_dat), baseline_vars),
               new = paste0(setdiff(names(comp_dat), baseline_vars), "_prev"))
    }
    
    # random quantiles
    comp_dat[, paste0("F_L",k,"|surv"):= runif(n_surv)]
    comp_dat[, paste0("F_W",k,"|L",k,"_surv"):= runif(n_surv)]
    comp_dat[, paste0("F_S",k,"|L",k,"_surv"):= runif(n_surv)]
    comp_dat[, paste0("F_Y|LW",k,"_surv"):= runif(n_surv)]
    
    if (k == 0){
      comp_dat[,L := qgamma(get(paste0("F_L",k,"|surv")), shape = 1 + 0.5*B1 + 0.5*B2, scale = 1)]
      comp_dat[,W := qbinom(get(paste0("F_W",k,"|L",k,"_surv")), size = 1, prob = expit(-0.2 +0.5 * B1 + 0.5*B2))]
    } else{
      comp_dat[,L := qgamma(get(paste0("F_L",k,"|surv")), shape = 1 + 0.5*B1 + 0.5*B2 +  0.1 * L_prev  - 0.5 * S_prev, scale = 1)]
      comp_dat[,W := qbinom(get(paste0("F_W",k,"|L",k,"_surv")), size = 1, prob = expit(-0.2 + 0.5 * B1 + 0.5*B2 + W_prev - 0.6* S_prev))]
    }
    
    comp_dat[,S := qbinom(get(paste0("F_S",k,"|L",k,"_surv")),size = 1, prob = expit(-1.5 + 0.5*B1 + 0.5*B2 + 0.5* L + 0.5 * W))]
    
    if (k >0){
      #get the distributions of Ls and Ws in the survivors
      for (j in (1:k)){
        if (j<k){
          comp_dat[,paste0("F_W",k-j,"|L",k-j,"_surv"):= cCopula(as.matrix(cbind(get(paste0("F_Y|L",k-j,"_W", k-j-1,"_surv_prev")),
                                                                                 get(paste0("F_W",k-j,"|L",k-j,"_surv_prev")))),copula = copW)[,2]]
          comp_dat[,paste0("F_L",k-j,"|surv"):= cCopula(as.matrix(cbind(get(paste0("F_Y|LW",k-j-1,"_surv_prev")),
                                                                        get(paste0("F_L",k-j,"|surv_prev")))),copula = copL)[,2]]
        } else{
          comp_dat[,paste0("F_W",k-j,"|L",k-j,"_surv"):= cCopula(as.matrix(cbind(get(paste0("F_Y|L",k-j,"_surv_prev")),
                                                                                 get(paste0("F_W",k-j,"|L",k-j,"_surv_prev")))),copula = copW)[,2]]
          comp_dat[,paste0("F_L",k-j,"|surv"):= cCopula(as.matrix(cbind(F_Y_surv_prev,
                                                                        get(paste0("F_L",k-j,"|surv_prev")))),copula = copL)[,2]]
        }
      }
    }
    
    for (j in (0:k)){
      if (j<k){
        comp_dat[,paste0("F_Y|L",k-j,"_W", k-j-1,"_surv"):= cCopula(as.matrix(cbind(get(paste0("F_W",k-j,"|L",k-j,"_surv")),
                                                                                    get(paste0("F_Y|LW",k-j,"_surv")))),copula = copW,inverse = T)[,2]]
        comp_dat[,paste0("F_Y|LW",k-j-1,"_surv"):= cCopula(as.matrix(cbind(get(paste0("F_L",k-j,"|surv")),
                                                                           get(paste0("F_Y|L",k-j,"_W", k-j-1,"_surv")))),copula = copL,inverse = T)[,2]]
      }else{
        comp_dat[,paste0("F_Y|L",k-j,"_surv"):= cCopula(as.matrix(cbind(get(paste0("F_W",k-j,"|L",k-j,"_surv")),
                                                                        get(paste0("F_Y|LW",k-j,"_surv")))),copula = copW,inverse = T)[,2]]
        comp_dat[,F_Y_surv:= cCopula(as.matrix(cbind(get(paste0("F_L",k-j,"|surv")),
                                                     get(paste0("F_Y|L",k-j,"_surv")))),copula = copL,inverse = T)[,2]]
      }
    }
    
    # hazard
    comp_dat[, lambda := exp(-2) *exp(betaS * S + betaB2 * B2 + betaSB2*B2*S)]
    comp_dat[, T := qexp(F_Y_surv, rate = lambda)]
    comp_dat[, Y := 1*(T<1)]
    
    dat[surv == 1, paste0("L","_",k) := comp_dat$L]
    dat[surv == 1, paste0("W","_",k) := comp_dat$W]
    dat[surv == 1, paste0("S","_",k) := comp_dat$S]
    dat[surv == 1, paste0("Y","_",k) := comp_dat$Y]
    dat[surv == 1, T := (comp_dat$T + k)]
    dat[,surv := 1*(get(paste0("Y_",k)) == 0)]
  }
  return(dat)
}
sim_func()
