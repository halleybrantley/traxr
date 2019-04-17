#' Rotate wind vector to downwind, crosswind, and vertical
#' @export
rotateUVW <- function(u, v, w){
  theta <-  atan2(mean(v, na.rm=TRUE),mean(u, na.rm=TRUE))
  u_rot1 <- u*cos(theta) + v*sin(theta)
  vp <- v*cos(theta) - u*sin(theta)
  phi <-  atan2(mean(w, na.rm=TRUE), mean(u_rot1, na.rm=TRUE))
  up <- u_rot1*cos(phi) + w*sin(phi)
  wp <- w*cos(phi) - u_rot1*sin(phi);

  psi <- -atan2(2*cov(vp,wp, use="complete"),
                (var(vp, na.rm=T)-var(wp, na.rm=T)))/2
  u_f <- up
  v_f <- vp*cos(psi) - wp*sin(psi)
  w_f <- vp*sin(psi) + wp*cos(psi)
  # mean(v_f*w_f, na.rm=T)
  # mean(v_f, na.rm=T)
  # mean(w_f, na.rm=T)
  return(list(uvw = data.frame(u=u_f, v=v_f, w=w_f),
              rot_angles = c(theta, phi, psi)))
}

rotateUV <- function(u, v){
  theta <-  atan2(mean(v, na.rm=TRUE),mean(u, na.rm=TRUE))
  up <- u*cos(theta) + v*sin(theta)
  vp <- v*cos(theta) - u*sin(theta)
  return(list(uv = data.frame(u=up, v=vp),
              rot_angles = c(theta)))
}
