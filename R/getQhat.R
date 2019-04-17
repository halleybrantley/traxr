
getQhat <- function(sims, grid, bases, design){
  for (i in sims){
    if(file.exists(sprintf("../Results/simulation/D%d/bases%d/param_%dm_%d.RData", 
                           design, bases, grid, i))){
      load(sprintf("../Results/simulation/D%d/bases%d/param_%dm_%d.RData", 
                   design, bases, grid, i))
      load(sprintf("../Data/simulated/D%d/Sim_%dm_%d.RData", 
                   design, grid, i))
      M <- dim(betaHat)[1]
      n <- length(rel_eff)
      B <- bs(seq(1,n,1), df = M, degree=3, intercept = TRUE)
      if (i == 1){
        Qhat <- array(NA, c(n, ncol(X), length(sims)))
      }
      Qhat[,,i] <- B%*%betaHat
    }
  }
  return(Qhat)
}