.onLoad <- function(libname, pkgname) {
  op <- options()
  op.sclaser <- list(sclaser.silent_startup = TRUE)
  toset <- !(names(op.sclaser) %in% names(op))
  if (any(toset)) options(op.sclaser[toset])
  invisible()
}

.onAttach <- function(libname, pkgname) {
  if (isTRUE(getOption("sclaser.silent_startup"))) return(invisible())

  ascii_logo <- paste0(
    "         _      _   ___ ___ ___ \n",
    "  ___ __| |    /_\\ / __| __| _ \\\n",
    " (_-</ _| |__ / _ \\\\__ \\ _||   /\n",
    " /__/\\__|____/_/ \\_\\___/___|_|_\\\n",
    "                                 \n",
    "     scLASER  single-cell longitudinal pipeline\n"
  )


  authors_str <- tryCatch({
    ar <- utils::packageDescription(pkgname, fields = "Authors@R")
    if (is.na(ar) || !nzchar(ar)) stop("no Authors@R")
    ppl <- eval(str2expression(ar))
    paste(format(ppl, include = c("given", "family")), collapse = ", ")
  }, error = function(...) {
    utils::packageDescription(pkgname, fields = "Author")
  })

  ver <- as.character(utils::packageVersion(pkgname))
  lab <- "Fan Zhang Lab (CU Anschutz)"
  lab_url <- "https://fanzhanglab.org"
  cite_hint <- "Run citation('scLASER') for how to cite."

  msg <- paste(
    ascii_logo,
    sprintf("Version: %s", ver),
    sprintf("Authors: %s", authors_str %||% "N/A"),
    sprintf("Lab: %s - %s", lab, lab_url),
    cite_hint,
    "Tip: options(sclaser.silent_startup = TRUE) to mute this banner.",
    sep = "\n"
  )

  packageStartupMessage(msg)
}


`%||%` <- function(x, y) if (is.null(x) || is.na(x) || !nzchar(x)) y else x
