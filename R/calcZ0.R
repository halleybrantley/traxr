psi <- function(x, L){
  if (L > 0){
    return(-4.8*x)
  } else {
    alpha <- (1-16*x)^(.25)
    2*log((1+alpha)/2) + log((1+alpha^2)/2) - 2*atan(alpha) + pi/2
  }
}


z0_obj <- function(z0, ubar, ustar, L, z, d){
  z0 <- exp(z0)
  (ubar - ustar/0.4 * (log((z-d)/z0) - psi((z-d)/L, L) + psi(z0/L, L)))^2
}


