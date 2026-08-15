# Package load hook: populate the method registry with the built-in methods so
# they are available before any synthesis runs.
.onLoad <- function(libname, pkgname) {
  register_builtin_methods()
  invisible()
}
