 
.onLoad = function(libname, pkgname) {
  library.dynam("pruning", pkgname)
  methods:::bind_activation(TRUE)
}

.onAttach <- function(libname, package) {
  methods:::bind_activation(FALSE)
}
