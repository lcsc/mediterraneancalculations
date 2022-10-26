#' @title mediterraneancalculationsNews
#' 
#' @description Show the NEWS file of the \pkg{mediterraneancalculations} package.
#'
#' @details (See description)
#' 
#' @export
#' 
ClimIndNews <- function() {
    file <- file.path(system.file(package="mediterraneancalculations"), "NEWS.md")
    file.show(file)
}