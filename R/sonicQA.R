#' Remove Spikes from Sonic Data
#' @description Used to remove sharp jumps in sonic data
#' @param spod data.frame with SPod data
#' @param nodes vector of letters with spod nodes
#' @export
sonicQA <- function(spod, nodes){
  spodQA <- spod
  for (node in nodes){
    spodQA[, paste0(node, ".Sonic.Temp..oC.")] <-
      as.numeric(as.character(spodQA[, paste0(node, ".Sonic.Temp..oC.")]))
    for (windVec in c("U", "V", "W")){
      spodQA[, paste0(node, ".Sonic.", windVec, "...")] <-
      as.numeric(as.character(spodQA[, paste0(node, ".Sonic.",
                                              windVec, "...")]))
    }
    spod[, paste0(node, ".Sonic.Temp..oC.")] <-
      zoo::na.locf( spod[,paste0(node, ".Sonic.Temp..oC.")])
    formula <- as.formula(paste0(node, ".Sonic.Temp..oC.~splines::ns(time,84)"))
    preds <- predict(lm(formula, spod))
    resid <- resid(lm(formula, spod))
    miss <- which(resid < -2.5)
    spodQA[miss, paste0(node, ".Sonic.Temp..oC.")] <- NA
    spodQA[miss, paste0(node, ".Sonic.", c("U", "V", "W"), "...")] <- NA
  }
  return(spodQA)
}
