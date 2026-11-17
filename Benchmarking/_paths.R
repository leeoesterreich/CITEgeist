# Path roots for the Benchmarking R code: cg_root("REPO"|"DATA"|"REF"|"ENVS"|"OUT").
# Reads CITEGEIST_<NAME>_ROOT; default is <repo>/<lowercase name>.
# cg_root("ENVS") defaults to <conda info --base>/envs when conda is on PATH, else NULL
# (meaning: use tools on PATH).
.cg_repo_root <- function() {
  here <- NULL
  for (i in rev(seq_len(sys.nframe()))) {
    of <- sys.frame(i)$ofile
    if (!is.null(of)) { here <- dirname(normalizePath(of)); break }
  }
  if (is.null(here)) {
    args <- commandArgs(trailingOnly = FALSE)
    f <- sub("^--file=", "", args[grepl("^--file=", args)])
    if (length(f)) here <- dirname(normalizePath(f[1]))
  }
  if (is.null(here)) {
    warning("_paths.R was not loaded via source() or Rscript; assuming the ",
            "repo root is getwd() = ", getwd(), ". Run from the repo root.",
            call. = FALSE)
    here <- file.path(getwd(), "Benchmarking")
  }
  normalizePath(file.path(here, ".."), mustWork = FALSE)
}
.cg_repo <- .cg_repo_root()

.cg_conda_envs <- function() {
  if (!nzchar(Sys.which("conda"))) return(NULL)
  base <- tryCatch(suppressWarnings(system2("conda", c("info", "--base"), stdout = TRUE, stderr = FALSE)),
                   error = function(e) character(0))
  if (length(base) == 0 || !nzchar(base[1])) return(NULL)
  file.path(base[1], "envs")
}

cg_root <- function(name) {
  if (!name %in% c("REPO", "DATA", "REF", "ENVS", "OUT")) {
    stop("unknown root '", name, "'; expected REPO, DATA, REF, ENVS or OUT")
  }
  value <- Sys.getenv(paste0("CITEGEIST_", name, "_ROOT"), unset = "")
  if (nzchar(value)) return(path.expand(value))
  if (name == "ENVS") return(.cg_conda_envs())
  if (name == "REPO") return(.cg_repo)
  file.path(.cg_repo, tolower(name))
}

# Python of conda env `name` under ENVS; bare "python" (on PATH) when ENVS is unset.
cg_env_python <- function(name) {
  envs <- cg_root("ENVS")
  if (is.null(envs)) return("python")
  file.path(envs, name, "bin", "python")
}

# ENVS/<...>; errors (naming the variable) when ENVS is unset.
cg_env_path <- function(...) {
  envs <- cg_root("ENVS")
  if (is.null(envs)) {
    stop("CITEGEIST_ENVS_ROOT is not set and conda was not found on PATH; set it to ",
         "the directory holding the conda environments.", call. = FALSE)
  }
  file.path(envs, ...)
}
