getQ_app <- function(grid, bases, days, node, tau){
  Q_all <- {}
  for (day in days){
    load(sprintf("../Results/application/param/June_%d_param_%dbases_%dm_%s.RData",
                 day, bases, grid, node))
    n <- length(y)
    betaHat <- apply(beta, c(2,3), quantile, tau)
    B <- bs(seq(1,n,1), df = m, degree=3, intercept = TRUE)
    Q <- B%*%betaHat
    colnames(Q) <- colnames(X)
    load(sprintf("../Data/processed/summary_10min_2018-06-%d_%s.RData",
                 day, node))
    Q <- as_tibble(Q)
    Q$time <- spod_stats$time
    Q_all <- bind_rows(Q_all, Q)
  }
  Q_all %>% gather("layer", "value", -time)
}
