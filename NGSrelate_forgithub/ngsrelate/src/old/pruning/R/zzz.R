.First.lib <- function(libname, package) {
  library.dynam("pruning", package,libname)
  methods:::bind_activation(TRUE)
}

.Last.lib <- function(libname, package) {
  methods:::bind_activation(FALSE)
}
