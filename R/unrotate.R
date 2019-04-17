#' @title Rotate back to orginal coordinate system
#' @param u_f along wind coordinate
#' @param v_f cross wind coordinate
#' @param w_f vertical coordianate
#' @param angles 3 rotation angles used to streamline coordinates
unrotate <- function(u_f, v_f, w_f, angles){
  angles <- -angles
  up <- u_f
  vp <- v_f*cos(angles[3]) - w_f*sin(angles[3])
  wp <- v_f*sin(angles[3]) + w_f*cos(angles[3])

  u_rot <- up*cos(angles[2]) + wp*sin(angles[2])
  w <- wp*cos(angles[2]) - up*sin(angles[2])

  u <- u_rot*cos(angles[1]) + vp*sin(angles[1])
  v <- vp*cos(angles[1]) - u_rot*sin(angles[1])

  return(uvw = data.frame(u=u, v=v, w=w))
}

unrotate_xy <- function(x, y, angle) {
  angle <- -angle
  x0 <- x*cos(angle) + y*sin(angle)
  y0 <- y*cos(angle) - x*sin(angle)
  return(xy = data.frame(x=x0, y=y0))
}
