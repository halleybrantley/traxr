#' Get Monin Obhukov Length
#' @export
getL <- function(ustar, flux, temp){
  temp <- unlist(temp)
  kv <- 0.4
  L <- -mean(temp+273.15, na.rm=T)*ustar^3/(kv*9.8*flux)
  return(L)
}

#' Get heat flux
#' @export
getFlux <- function(uvw, temp){
  uvw <- as.data.frame(uvw)
  wtheta <- cov(uvw$w, unlist(temp), use = "complete.obs")
  return(wtheta)
}

#' Get turbulence
#' @export
getUstar <- function(uvw){
  uvw <- as.data.frame(uvw)
  ustar <- sqrt(abs(cov(uvw$u, uvw$w, use = "complete.obs")))
  return(ustar)
}



