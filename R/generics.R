#### Generic Functions ---------------------------------------------------------------####

#----------------------------------------------------------------------------------------#
# Author: Pedro Nascimento de Lima
# Purpose: This File implement common generic functions for the gerbil class.
# Creation Date: Sept 2020
#----------------------------------------------------------------------------------------#

#' Prints a \code{gerbil} object. Printed output includes a variable-by-variable summary of variable types and missingness rates. The implemented predictor matrix is also provided.
#'
#' @param x object of \code{gerbil} class
#' @param ... additional parameters to be passed down to inner functions.
#'
#' @return The functions \code{print.gerbil} and \code{summary.gerbil} display information about the \code{gerbil} object.
#'   Primarily, the variable type and missingness rate are displayed for each variable. The predictor matrix is also provided. 
#' 
#' @export
print.gerbil = function(x, ...) {
  m <- length(x$imputed)
  mcmciter <- dim(x$chainSeq[[1]])[2] - 1
  info <- x$summary
  info.tmp <- info[, "Miss.Rate"]
  info.tmp <- formatC(100 * info.tmp, format = "f", digits = 2)
  info.tmp <- paste0(info.tmp, "%")
  info1 <- info
  info1[, "Miss.Rate"] <- info.tmp
  if (is.element("Cens.Rate", colnames(info))) {
    info.tmp <- info[, "Cens.Rate"]
    info.tmp <- formatC(100 * info.tmp, format = "f", digits = 2)
    info.tmp <- paste0(info.tmp, "%")
    info1[, "Cens.Rate"] <- info.tmp
  }

  cat("Object Class: gerbil \n\n")
  if (m == 1) {
    cat("Includes ", m, " imputed dataset created using ", mcmciter, " iterations of MCMC. \n\n", sep = "")
  } else {
    cat("Includes ", m, " imputed datasets, each created using ", mcmciter, " iterations of MCMC. \n\n", sep = "")
  }
  cat("Predicted Variables, Types and Missing Rates: \n")
  print(info1)
  cat("\n")
  cat("Predictor Matrix: \n")
  print(x$predMat.initial)
}


#' Summarises a \code{gerbil} object. Printed output includes a variable-by-variable summary of variable types and missingness rates. The implemented predictor matrix is also provided.
#'
#' @param object An object of \code{gerbil} class
#' @param ... Additional parameters to be passed down to inner functions.
#'
#' @return The functions \code{print.gerbil} and \code{summary.gerbil} display information about the \code{gerbil} object.
#'   Primarily, the variable type and missingness rate are displayed for each variable. The predictor matrix is also provided. 
#' 
#' @export
#'
summary.gerbil = function(object, ...) {
  print(object, ...)
}


#' @title Extracting imputed datasets from gerbil objects
#'
#' @description
#' Using a \code{gerbil} object as an input, this function returns imputed datasets.
#'
#' @details
#' The function either return a single imputed dataset (if \code{imp} is a scalar) or a tall dataset if (if \code{imp} is a vector) with the individual datsets stacked on top of each other.
#'
#' @param gerb A \code{gerbil} object containing the imputed data.
#' @param imp The imputed datasets which are to be returned (defaults to \code{imp = 1}).  Letting \code{m} indicate the number of imputed datasets contained in \code{gerb}, \code{imp} should be a subset of \code{1:m}.
#'    If \code{imp} is a scalar, a single imputed dataset is returned.  If \code{imp} is a vector, then the individual datasets are stacked on top of each other and returned.
#'
#' @return \code{imputed()} returns a data frame or matrix. If \code{imp} has multiple elements, columns are added to indicate the imputation number and the case ID.
#'
#' @examples
#' \donttest{
#' #Load the India Human Development Survey-II dataset
#' data(ihd_mcar)
#' 
#' # Create a gerbil object
#' imps.gerbil <- gerbil(ihd_mcar, m = 5, ords = "education_level", semi = "farm_labour_days", 
#'        bincat = "job_field", n.cores = 1)
#'
#' # Return a single imputed datasets
#' imp.gerb <- imputed(imps.gerbil, imp = 2)
#'
#' # Return multiple (stacked) datasets
#' imp.gerb <- imputed(imps.gerbil, imp = 1:5)
#' }
#'
#' @export
#'
imputed <- function(gerb, imp = 1) {

  k <- length(imp)
  m <- length(gerb$imputed)
  n <- NROW(gerb$imputed[[1]])
  p <- NCOL(gerb$imputed[[1]])

  if(length(setdiff(imp, 1:m)) > 0) {
    warning("imp contains invalid elements.")
    imp <- intersect(imp, 1:m)
    k <- length(imp)
  }

  if (k > 1) {
    if (class(gerb$imputed[[1]]) == "matrix") {
      out <- matrix(NA, k * n, p + 2)
    } else {
      out <- data.frame(matrix(NA, k * n, p + 2))
    }
    colnames(out) <- c("imp", "ID", colnames(gerb$imputed[[1]]))

    nams2 <- c(t(matrix(imp, k, n)))
    if(length(rownames(gerb$imputed[[1]])) > 0) {
      nams1 <- rep(rownames(gerb$imputed[[1]]), k)
      rownames(out) <- paste0(nams1, ".", nams2, sep = "")
    } else {
      nams1 <- 1:n
    }

    out[, 1] <- nams2
    out[, 2] <- nams1

    for(i in 1:k) {
      out[((i - 1) * n + 1):(i * n), 3:(p + 2)] <- gerb$imputed[[imp[i]]]
    }
  } else {
    out <- gerb$imputed[[imp]]
  }

  return(out)

}


#' @title Write imputed datasets from gerbil objects to a file or files
#'
#' @description
#' Using a \code{gerbil} object as an input, this function writes imputed datasets to an output file.
#'
#' @details
#' The function writes imputed datasets to either an Excel (.xlsx) or a CSV (.csv) file, depending upon the extension of the parameter \code{file}.
#'   No other file types are supported.
#'   To write multiple imputed datasets simultaneously, specify \code{imp} as a vector with length greater than 1.
#'   Multiple imputed datasets are either written in a stacked format (if \code{tall = TRUE}) or written separately (if \code{tall = FALSE}).
#'
#' @param gerb A \code{gerbil} object containing the imputed data.
#' @param file The name of the file to which the imputed datasets are to be written.
#'   Which type of file (.xlsx or .csv) is created depends upon the extension of the parameter \code{file}.
#' @param imp The imputed datasets which are to be written.  Can be a scalar or, if multiple imputed datasets are to be written, a vector.
#'   All elements of \code{imp} should be integers greater than 0 but no greater than \code{m}, which is the number of imputed datasets in \code{gerb}.
#'   \code{imp} defaults to \code{1:m}.
#' @param tall A logical expression indicating whether the datasets are to be written in a tall (stacked) format
#'   or written separately.  When writing to an XLSX file with \code{tall = FALSE}, one tab is created for each imputed dataset.
#'   When writing to a CSV file with \code{tall = FALSE}, one file is created for each imputed dataset
#'   (in this case, several file names will be created from the base string given by \code{file}).
#' @param row.names A logical value indicating whether the row names of the datasets are to be written.
#' 
#' @return No returned value, but instead a data file is written to a specified directory.
#' 
#' @import openxlsx
#'
#' @examples
#' \donttest{
#' #Load the India Human Development Survey-II dataset
#' data(ihd_mcar)
#' 
#' # Create a gerbil object
#' imps.gerbil <- gerbil(ihd_mcar, m = 5, ords = "education_level", semi = "farm_labour_days", 
#'        bincat = "job_field", n.cores = 1)
#' 
#' # Write all imputed datasets to separate CSV files
#' write.gerbil(imps.gerbil, file.path(tempdir(), "gerbil_example.csv"), imp = 1:5, tall = FALSE)
#'
#' # Write all imputed datasets to a single CSV files
#' write.gerbil(imps.gerbil, file.path(tempdir(), "gerbil_example.csv"), imp = 1:5, tall = TRUE)
#'
#' # Write all imputed datasets to an XLSX file
#' write.gerbil(imps.gerbil, file.path(tempdir(), "gerbil_example.xlsx"), imp = 1:5, tall = FALSE)
#' }
#'
#' @export
#'
write.gerbil <- function(gerb, file = NULL, imp = NULL, tall = FALSE, row.names = FALSE) {

  m <- length(gerb$imputed)
  if (length(imp) == 0) {
    imp <- 1:m
  }
  k <- length(imp)
  n <- NROW(gerb$imputed[[1]])
  p <- NCOL(gerb$imputed[[1]])

  if (length(file) == 0) {
    file <- "gerbil_imputed.csv"
  }

  if (length(setdiff(imp, 1:m)) > 0) {
    warning("imp contains invalid elements.")
    imp <- intersect(imp, 1:m)
    k <- length(imp)
  }

  if (k == 0) {
    stop("Specify valid values for imp.")
  }

  if (substr(file, nchar(file) - 3, nchar(file)) == ".csv") {
    file <- substr(file, 1, nchar(file) - 4)
    use.xlsx <- FALSE
  } else if (substr(file, nchar(file) - 3, nchar(file)) == ".xls") {
    file <- substr(file, 1, nchar(file) - 4)
    use.xlsx <- TRUE
  } else if (substr(file, nchar(file) - 4, nchar(file)) == ".xlsx") {
    file <- substr(file, 1, nchar(file) - 5)
    use.xlsx <- TRUE
  } else {
    use.xlsx <- FALSE
  }

  if (use.xlsx) {
    suff <- ".xlsx"
  } else {
    suff <- ".csv"
  }

  if (k > 1 & !tall & !use.xlsx) {
    files <- paste0(file, imp, suff, sep = "")
  } else {
    files <- paste0(file, suff, sep = "")
  }

  if (use.xlsx) {
  # It is more straightforward to always require the package and use that as a
  # dependency.
  #   xlsx.loaded <- is.element("openxlsx", loadedNamespaces())
  #   if (!requireNamespace("openxlsx", quietly = TRUE)) {
  #     stop("The openxlsx package is needed when saving multiple datasets to a single file or when the file name specifies a .xlsx extension. Please install openxlsx, or if you'd like to write to .csv, only select a single post-intervention time and append the filename appropriately.",
  #          call. = FALSE)
  #   }
  #   if (!xlsx.loaded) {
  #     attachNamespace("openxlsx")
  #   }

    if (!tall) {
      wb <- openxlsx::createWorkbook()
      for(i in 1:k) {
        sheet <- openxlsx::addWorksheet(wb, sheetName = paste0("Imputed ", imp[i]))
        openxlsx::writeData(wb, sheet = i, gerb$imputed[[imp[i]]], rowNames = row.names, colNames = TRUE, startRow = 1)
      }
      openxlsx::saveWorkbook(wb, file = files, overwrite = TRUE)
      rm(wb)
    } else {
      out <- imputed(gerb, imp)
      wb <- openxlsx::createWorkbook()
      if (k > 1) {
        s.nam <- "Imputed (tall format)"
      } else {
        s.nam <- paste0("Imputed ", imp)
      }
      sheet <- openxlsx::addWorksheet(wb, sheetName = s.nam)
      openxlsx::writeData(wb, sheet = 1, out, rowNames = row.names, colNames = TRUE, startRow = 1)
      openxlsx::saveWorkbook(wb, file = files, overwrite = TRUE)
      rm(wb)
    }
  } else {
    if (!tall) {
      for(i in 1:k) {
        utils::write.csv(gerb$imputed[[imp[i]]], files[i], na = "", row.names = row.names)
      }
    } else {
      out <- imputed(gerb, imp)
      utils::write.csv(out, files, na = "", row.names = row.names)
    }
  }

}


