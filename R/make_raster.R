#' Makes raster grid for aggregating trajecotries
#'
#' \code{make_grid}
#'
#' @param site_coords vector of lon, lat
#' @param grid_res resolution of grid cells
#' @param max_dist maximum distance from center point in cardinal directions
#' @export
make_raster <- function(site_coords, grid_res, max_dist){
  if (site_coords[2] < 0) {
    warning("Latitude is negative, longitude should be first element of site_coords.")
  }
  site_coords <- matrix(site_coords, nrow=1, ncol=2)
  colnames(site_coords) <- c("lon", "lat")
  coord.deg = SpatialPoints(site_coords,
                            proj4string=CRS("+proj=longlat"))
  
  center_m <- coordinates(spTransform(coord.deg, CRS("+init=epsg:26916")))
  rast_grid <- raster(xmn=center_m[1]-max_dist, xmx=center_m[1]+max_dist,
                      ymn=center_m[2]-max_dist, ymx=center_m[2]+max_dist,
                      res=grid_res)
  values(rast_grid) <- 1:ncell(rast_grid)
  return(rast_grid)
}
