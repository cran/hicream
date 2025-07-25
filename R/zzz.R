#' @importFrom reticulate py_require

sklearn <- NULL
kneebow <- NULL
pandas <- NULL
numpy <- NULL

.onLoad <- function(libname, pkgname) {
  reticulate::py_require("scikit-learn")
  sklearn <<- reticulate::import("scikit-learn", delay_load = TRUE)
  reticulate::py_require("kneebow")
  kneebow <<- reticulate::import("kneebow", delay_load = TRUE)
  reticulate::py_require("pandas")
  pandas <<- reticulate::import("pandas", delay_load = TRUE)
  reticulate::py_require("numpy")
  numpy <<- reticulate::import("numpy", delay_load = TRUE)
}
