#' Register a decision function
#'
#' @param decisionFunction string or function
#'
#' @return decision function to be used
#'
#' @keywords internal
#'
registerDecFun <- function(decisionFunction) {
  if (is.character(decisionFunction)) {
    lowercase <- tolower(decisionFunction)
    split <- strsplit(lowercase, "")[[1]]
    if (split[1] == "c") {
      if (lowercase != "corminuspartcor") {
        message("Using the 'corMinusPartCor' decision function.")
      }
      return(corMinusPartCor)
    } else if (split[1] == "p") {
      if (lowercase != "prodcoef") {
        message("Using the 'prodCoef' decision function.")
      }
      return(prodCoef)
    } else {
      stop("Specified decisionFunction not found.")
    }
  } else if (!is.function(decisionFunction)) {
    stop("Please specify a decisionFunction.")
  } else {
    return(decisionFunction)
  }
}


#' Generate a (random) start
#'
#' @param n total number of mediators
#' @param m number of selected mediators
#'
#' @return binary vector of mediator selections
#'
#' @keywords internal

generateStart <- function(n, m) {
  out <- integer(n)
  out[sample(1:n, m, replace = FALSE)] <- 1L
  return(out)
}

#' Register a printing function for verbose output
#'
#' @param w half the width of the console
#' @param m the number of mediators
#'
#' @return callable printing function
#'
#' @keywords internal

registerCatFun <- function(w, m) {
  if (requireNamespace("crayon", quietly = TRUE)) {
    if (w > m) {
      catmsel <- function(msel, w) {
        prnt <- paste(msel)
        prnt[msel == 1] <- crayon::bgWhite(crayon::black("1"))
        prnt[msel == 0] <- crayon::silver("0")
        cat(prnt,"\n")
      }
    } else {
      catmsel <- function(msel, w) {
        prnt <- paste(msel[1:w])
        prnt[msel[1:w] == 1] <- crayon::bgWhite(crayon::black("1"))
        prnt[msel[1:w] == 0] <- crayon::silver("0")
        cat(prnt,"\n")
      }
    }
  } else {
    if (w > m) {
      catmsel <- function(msel, w) cat(msel,"\n")
    } else {
      catmsel <- function(msel, w) cat(msel[1:w],"\n")
    }
  }
  return(catmsel)
}


#' Register an lapply function for parallel or serial processing with or
#' without progress bar.
#'
#' @param parallel true/false
#' @param progressBar true/false
#'
#' @return drop-in replacement lapply function
#'
#' @keywords internal

registerLapplyFun <- function(parallel, progressBar, nCores) {
  if (parallel) {
    if (progressBar) {
      lapplyfun <- function(X, FUN) {
        clus <- parallel::makeCluster(nCores)
        on.exit(parallel::stopCluster(clus))

        # initialise the cluster with the exact environment of the cmf call
        parallel::clusterEvalQ(clus, {library(cmfilter)})
        vars <- ls(envir = parent.frame())
        parallel::clusterExport(clus, vars[vars != "clus"], parent.frame())

        # run the apply function
        out <- pbapply::pblapply(X = X, FUN = FUN, cl = clus)
        parallel::stopCluster(clus)
        return(out)
      }
    } else {
      lapplyfun <- function(X, FUN) {
        clus <- parallel::makeCluster(nCores)
        on.exit(parallel::stopCluster(clus))

        # initialise the cluster with the exact environment of the cmf call
        parallel::clusterEvalQ(clus, {library(cmfilter)})
        vars <- ls(envir = parent.frame())
        parallel::clusterExport(clus, vars[vars != "clus"], parent.frame())

        # run the apply function
        out <- parallel::parLapply(X = X, fun = FUN, cl = clus)
        parallel::stopCluster(clus)
        return(out)
      }
    }
  } else {
    if (progressBar) {
      lapplyfun <- pbapply::pblapply
    } else {
      lapplyfun <- lapply
    }
  }
  return(lapplyfun)
}
