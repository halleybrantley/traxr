#' @useDynLib traxr, .registration = TRUE
#' @importFrom Rcpp sourceCpp
coreModel <- function(u, v, w, zSens, ustar, L, Zo, bw, sigmaUustar,
                      sigmaVustar, kv, C0, alpha, MaxFetch, outSteps = 50){

	Linv <- 1/L

  if (Linv < 0) {
    out <- csFi(u, v, w, zSens, ustar, Linv, Zo, bw,
                sigmaUustar, sigmaVustar, kv, C0, alpha, MaxFetch, outSteps)
  } else {
    out <- csFs(u, v, w, zSens, ustar, Linv, Zo, bw,
                sigmaUustar, sigmaVustar, kv, C0, alpha, MaxFetch, outSteps)
  }

  return(out)
}


