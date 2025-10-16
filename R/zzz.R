.onLoad <- function(libname, pkgname) {
  op <- options()
  op.sclaser <- list(sclaser.silent_startup = FALSE)
  toset <- !(names(op.sclaser) %in% names(op))
  if (any(toset)) options(op.sclaser[toset])
  invisible()
}

.onAttach <- function(libname, pkgname) {
  if (isTRUE(getOption("sclaser.silent_startup"))) return(invisible())

  ascii_logo <- paste0(
    "   ____   _____    _                _____  \n",
    "  / ___| / ____|  / \\    ___  ___  / ____| \n",
    "  \\___ \\| (___   / _ \\  / __|/ _ \\| |      \n",
    "   ___) |\\___ \\ / ___ \\ \\__ \\  __/| |___   \n",
    "  |____/ |____//_/   \\_\\|___/\\___| \\____|  \n",
    "      scLASER — single-cell longitudinal pipeline\n"
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
    sprintf("Lab: %s — %s", lab, lab_url),
    cite_hint,
    "Tip: options(sclaser.silent_startup = TRUE) to mute this banner.",
    sep = "\n"
  )

  packageStartupMessage(msg)
}


`%||%` <- function(x, y) if (is.null(x) || is.na(x) || !nzchar(x)) y else x
