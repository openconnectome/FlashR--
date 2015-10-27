.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to FlashR!")
}

.onLoad <- function(libname, pkgname) {
  op <- options()
  op.flashR <- list(
    flashR.path = "~/R",
    flashR.install.args = "",
    flashR.name = "FlashR",
    flashR.desc.author = 'person("Gregory", "Kiar", "gkiar@jhu.edu", role=c("aut", "cre"))',
    flashR.desc.license = "Apache 2.0",
    flashR.desc.suggests = NULL,
    flashR.desc = list()
  )
  toset <- !(names(op.flashR) %in% names(op))
  if(any(toset)) options(op.flashR[toset])

  invisible()
}
