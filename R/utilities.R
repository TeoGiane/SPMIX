#' Deserialize messages of the \code{SPMIX} package
#'
#' This funciton is a wrapper to deserialize raw vectors using the \code{.proto} files available in this package.
#' We rely on Google Protocol Buffers for the serialization procedure and on \code{\link{RProtoBuf}} package to
#' provide and easy-to-use interface for \code{R} users
#'
#' @param message_type A string containing the name of the Message in which the serialized message will be converted
#' @param raw_vector A vector of type \code{raw}. The message to be unserialized
#' @return An object of class \code{RProtoBuf::Message}, which stores the unserialized message and can be manipulated
#' using the \link{RProtoBuf} package.
#'
#' @export
DeserializeSPMIXProto <- function(message_type, raw_vector) {

  # Check Message Descriptor
  if (message_type == "EigenMatrix") {
    RProtoBuf::readProtoFiles(system.file("proto/eigen.proto", package = "SPMIX"))
  } else if (message_type == "spmix.SamplerParams") {
    RProtoBuf::readProtoFiles(system.file("proto/sampler_params.proto", package = "SPMIX"))
  } else if (message_type == "spmix.UnivariateState") {
    RProtoBuf::readProtoFiles(system.file("proto/univariate_mixture_state.proto", package = "SPMIX"))
  } else if (message_type == "spmix.OptimOptions") {
    RProtoBuf::readProtoFiles(system.file("proto/optimization_options.proto", package = "SPMIX"))
  } else {
    stop("Input 'message_type' is of uknown type")
  }

  # Read state from proper descriptor
  state <- RProtoBuf::read(get(message_type), raw_vector)
  return(state)
}


#' Compute the geometry of the boundaries between areal units
#'
#' This utility takes as input a list of boundaries between areal units and the
#' corresponding \code{sf} geometry object and computes the geometry of the boundaries
#' between the areal units.
#'
#' @param boundary_list A list of length \mjseqn{I}, where \mjseqn{I} is the number of areal units.
#' Element \mjseqn{i} of the list contains a vector of indices representing the areal units
#' that share a boundary with areal unit \mjseqn{i}.
#' @param sf_geometry An \code{sf} object representing the geometry of the areal units.
#' @return An \code{sf} object containing the geometry of the boundaries between the areal units.
#'
#' @export
BoundaryGeometry <- function(boundary_list, sf_geometry) {

  # Check inputs
  if (!is.list(boundary_list)) {
    stop("'boundary_list' must be a list")
  }
  if (!inherits(sf_geometry, "sf")) {
    stop("'sf_geometry' must be an sf object")
  }
  if (nrow(sf_geometry) != length(boundary_list)) {
    stop("Length of 'boundary_list' must match number of rows in 'sf_geometry'")
  }

  # Add id column if not present
  if(!("id" %in% names(sf_geometry))){
    sf_geometry$id <- 1:nrow(sf_geometry)
  }

  # Create empty list to store boundary geometries
  geom_bdd <- list()

  # Populate list
  for(i in 1:nrow(sf_geometry)) {
    # Get current area and its boundaries
    if (length(boundary_list[[i]]) > 0) {
      for (j in boundary_list[[i]]) {
        # Compute geometry of boundary
        sel_geom <- sf_geometry[c(i, j), "id"]
        bounds <- suppressWarnings(sf::st_intersection(sel_geom, sel_geom))
        bounds <- sf::st_union(sf::st_geometry(bounds[bounds$id != bounds$id.1, ]))
        # Append to list
        geom_bdd <- append(geom_bdd, list(sf::st_sf(geometry = bounds)))
      }
    }
  }
  # Bind all objects
  geom_bdd <- do.call(rbind, geom_bdd)

  # Drop points if present
  points <- which(attr(geom_bdd$geometry, "classes") == "POINT")
  if(length(points) > 0){
    geom_bdd <- geom_bdd[-points, ]
  }

  # Condense everything into a unique sf object and return
  return(geom_bdd)
}


#' Transform ggmap object to EPSG:3857
#'
#' This utility takes as input a \code{ggmap} object and transforms its bounding box
#' to EPSG:3857 coordinate reference system.
#'
#' @param map Output of \code{ggmap::get_map} function
#' @return The same map, but with a bounding box suitable to be plot in \code{ggplot}, i.e., in EPSG:3857
#' coordinate reference system.
#'
#' @export
sf_ggmap <- function(map) {
  if (!inherits(map, "ggmap")) stop("map must be a ggmap object")

  # Extract the bounding box (in lat/lon) from the ggmap to a numeric vector,
  # and set the names to what sf::st_bbox expects:
  map_bbox <- setNames(unlist(attr(map, "bb")),
                       c("ymin", "xmin", "ymax", "xmax"))

  # Coonvert the bbox to an sf polygon, transform it to 3857, and convert back to a bbox
  bbox_3857 <- sf::st_bbox(sf::st_transform(sf::st_as_sfc(sf::st_bbox(map_bbox, crs = 4326)), 3857))

  # Overwrite the bbox of the ggmap object with the transformed coordinates
  attr(map, "bb")$ll.lat <- bbox_3857["ymin"]
  attr(map, "bb")$ll.lon <- bbox_3857["xmin"]
  attr(map, "bb")$ur.lat <- bbox_3857["ymax"]
  attr(map, "bb")$ur.lon <- bbox_3857["xmax"]

  # Return
  return(map)
}
