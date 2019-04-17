
#' Makes raster grid for aggregating trajecotries
#'
#' \code{make_grid}
#'
#' @param site_coords vector of lon, lat
#' @param grid_res resolution of grid cells
#' @param max_dist maximum distance from center point in cardinal directions
#' @export
grid_to_latlon <- function(rast_grid, deg = TRUE){
  grid <- rasterToPoints(rast_grid)
  grid_m <- SpatialPoints(cbind(grid[,1], grid[,2]),
                          proj4string=CRS("+init=epsg:26916"))

  if (deg){
    grid_deg <- as.data.frame(coordinates(spTransform(grid_m,
                                                    CRS("+proj=longlat"))))
  } else {
    grid_deg <- as.data.frame(coordinates(grid_m))
  }

  names(grid_deg) <- c("lon", "lat")
  grid_deg$layer <- grid[,3]
  return(grid_deg)
}
