#' Get backward trajectories
#' @export
getTraj <- function(sigmaU, sigmaV, sigmaW, ustar, ubar, flux, L, u, v, w, nTraj,
                    outSteps, height_sens, MaxFetch=3000){

  zSens <- height_sens
  ustar2 <- ustar^2
  sigmaUu <- sqrt(sigmaU/ustar2)
  sigmaVu <- sqrt(sigmaV/ustar2)
  sigmaWu <- sqrt(sigmaW/ustar2)
  # z_sWu = z_sens (heigth of wind measurements) defined in gen_interval
  bw <- calcbw(sigmaWu, zSens/L)
  C0 <- calcC0(bw, 0.4, A = 1)
  if (C0 > 9){
    C0 <- 9
  } else if (C0 < 2){
    C0 <- 2
  }
  alpha <- .02
  uv <- rotateUV(u, v)
  Zo <- exp(optimize(z0_obj, ubar=ubar, ustar=ustar, L=L, z = zSens, d=.1,
                     interval = c(-20, 20))$minimum)
  Sigma <- matrix(c(sigmaU, 0, -ustar2, 0, sigmaV, 0, -ustar2, 0, sigmaW),
                  nrow=3)
  Sigma_sqrt <- t(chol(Sigma))
  wind_mat <- t(Sigma_sqrt%*%matrix(rnorm(3*nTraj), nrow=3))

  output <- coreModel(wind_mat[,1]+ubar,
                      wind_mat[,2],
                      wind_mat[,3],
                      zSens, ustar, L, Zo,
                      bw, sigmaUu, sigmaVu,
                      kv=0.4, C0, alpha=alpha, MaxFetch, outSteps)


  traj <- data.frame(time = output$TimeOut, ID = output$Traj_IDOut,
                     x=output$xOut, y=output$yOut, z=output$zOut)
  # traj <- traj %>%
  #   mutate(dist = sqrt(x^2+y^2)) %>%
  #   filter(dist > 10)

  xyz <- unrotate_xy(traj$x, traj$y, uv$rot_angle)
  xyz$z <- traj$z
  xyz$ID <- traj$ID
  #plot(y~x, xyz)
  return(xyz)
}
