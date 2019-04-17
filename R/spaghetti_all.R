spaghetti_all <- function(Qhat, sourceInd, true_Q, design = 1){
  if (design == 1){
    colPal <- c("grey70", "#99d594")
  } else if (design == 2){
    colPal <- c("grey70", "#3288bd")

  } else {
    colPal <- c("grey70", "#3288bd", "#99d594")
  }
  n <- nrow(Qhat)
  Qsource <- {}
  for (i in 1:ncol(Qhat)){
    Qsource <- bind_rows(Qsource,
                         data.frame(Q = Qhat[,i,],
                                    cell = i,
                                    time =  seq(1, n, 1)))
  }
  nSim <- dim(Qhat)[3]
  Qsource$mean <- rowMeans(Qsource[,1:nSim], na.rm=T)
  Qsource$true <- as.numeric(true_Q)

  Qlong <- gather(Qsource, id, value, -c(time, cell))
  Qlong$col <- 0
  Qlong$col[which(Qlong$id == "mean")] <- 1
  Qlong$col[which(Qlong$id == "true")] <- 2
  Qlong$source <- 0
  Qlong$source[which(Qlong$cell == sourceInd[1])] <- 1
  Qlong$source[which(Qlong$cell == sourceInd[2])] <- 2
  Qlong$source <- factor(Qlong$source)
  Qlong$cellid <- str_c(Qlong$cell, Qlong$id)
  Qlong %>%
    ggplot(aes(x=time, y=value, group = cellid, col = source)) +
    geom_line(alpha = 0.2) +
    theme_bw() +
    geom_line(data=Qlong[which(Qlong$id == "mean"),],
              aes(group = cell), size = 1, linetype = 2) +
    geom_line(data=Qlong[which(Qlong$id == "true"),],
              size = 1)+
    scale_color_manual(values = colPal) +
    guides(col = "none")+
    ylim(c(0,5.5)) +
    labs(y="Source Strength", x = "Time")
}

