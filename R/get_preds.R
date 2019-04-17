get_preds <- function(grid, bases, days, node){
  spod <- {}
  for (day in days){
    load(sprintf("../Results/application/param/June_%d_param_%dbases_%dm_%s.RData", 
                 day, bases, grid, node))
    muHat <- apply(mu, 2, mean)
    load(sprintf("../Data/processed/summary_10min_2018-06-%d_%s.RData", 
                 day, node))
    spod_stats$pred <- muHat
    spod <- bind_rows(spod, spod_stats)
  }
  spod
}