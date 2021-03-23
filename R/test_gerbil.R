#' @title Goodness-of-fit testing for \code{gerbil} objects
#'
#' @description
#' Using a \code{gerbil} object as an input, this function performs univariate and bivariate goodness-of-fit tests
#' to compare distributions of imputed and observed values.  
#'
#' @details
#' Goodness of fit is determined using contingency tables of counts across categories of the corresponding variable(s). 
#' For univariate testing (\code{type = 1}), a one-way table is calculated for observed cases and compared to an analogous table for imputed cases, 
#' whereas for bivariate testing (\code{type = 2}), two-way tables are calculated.  
#' Continuous variables are binned according to cut-points defined using the parameter \code{breaks}. 
#' Tests are performed using one of three methods (determined from the parameter \code{method}): 1) Chi-squared (the default); 2) Fisher's exact; and 3) A G-test. 
#' G-testing is implemented via the function \code{GTest()} from the \code{DescTools} package.
#' Note that for univariate testing of continuous variables, a Kolmogorov-Smirnov test may be performed instead by setting \code{ks = TRUE}.
#' 
#' The only required input is a parameter \code{x} which is a \code{gerbil} object.
#' 
#' Note that univariate differences between observed and imputed data may be explained by the missingness mechanism and are not necessarily indicative of poor imputations.
#' Note also that most imputation methods like gerbil (and mice and related methods) are not designed to capture complete bivariate distributions. As such, the bivariate tests may be likely to return small p-values.
#'
#' @param x A \code{gerbil} object containing the imputed data.
#' @param y A vector listing the column names of the imputed data for which tests should be run. See details. By default, \code{y} contains all columns of the data that required imputation.
#' @param type A scalar used to specify the type of tests that will be performed.  Options include univariate (marginal) tests (\code{type = 1}) and bivariate tests (\code{type = 2}). See details. Defaults to \code{type = 1}.
#' @param imp A scalar or vector indicating which of the multiply imputed datasets should be used for testing.  Defaults to \code{imp = 1}. 
#' @param breaks Used to determine the cut-points for binning of continuous variables into categories. Ideally, \code{breaks} is a named list, where the list names are the names of the continuous variables. 
#'               Each element of the list can be a vector giving the respective cutpoints or a scalar which is used to indicate the number of bins (in which case cutpoints are determined from percentiles in order to yield bins of approximately equal size). 
#'               If \code{breaks} is a scalar or a vector (and not a list), the binning strategy indicated by \code{breaks} is applied to each variable in accordance with the description above.  
#'               Defaults to \code{breaks = 4}. 
#' @param method The type of test that is used to compare contingency tables.  Options include \code{'chi-squared'} for chi-squared testing (the default), \code{'fisher'} for Fisher's exact test, and \code{'G'} for a G-test.  
#' @param ks If \code{TRUE}, a Kolmogorov-Smirnov test is used when for univariate comparisons with continuous variables.  This functionality is not enabled for bivariate testing.  Defaults to \code{FALSE}.
#' @param partial Indicates how partially imputed pairs are handled in bivariate testing. If \code{'imputed'}, cases with at least one missing variable in a pair are considered imputed. Otherwise (\code{partial = 'observed'}), only cases with both variables in the pair missing are considered imputed.  
#' @param ... Arguments to be passed to methods.
#' 
#' @return \code{gof_gerbil()} returns an object of the class \code{gof_gerbil} that has following slots:
#' 
#' \describe{
#'         \item{Stats}{A vector (when \code{type = 1}) or matrix (when \code{type = 2}) giving the value of the test statistic (or coefficient) for the corresponding variable (or variable pair).}
#'         \item{p.values}{A vector (when \code{type = 1}) or matrix (when \code{type = 2}) giving the value of the p-value for the test applied to the corresponding variable (or variable pair).}
#'         \item{Test}{A vector (when \code{type = 1}) or matrix (when \code{type = 2}) indicating the type of test applied to the corresponding variable (or variable pair).} 
#'         \item{Breaks}{A list giving the cutpoints used for binning each continuous or semi-continuous variable.}
#' }
#' 
#' @examples
#' \donttest{
#' #Load the India Human Development Survey-II dataset
#' data(ihd_mcar) 
#' 
#' imps.gerbil <- gerbil(ihd_mcar, m = 1, mcmciter = 200, ords = "education_level", 
#'        semi = "farm_labour_days", bincat = c("sex", "marital_status", "job_field", "own_livestock"))
#' 
#' #Run univariate tests
#' tests.gerbil.uni <- gof_gerbil(imps.gerbil, imp = 1, type = 1)
#' 
#' #Print a summary
#' tests.gerbil.uni
#' 
#' #Run bivariate tests
#' tests.gerbil.bi <- gof_gerbil(imps.gerbil, imp = 1, type = 2)
#' 
#' #Print a summary
#' tests.gerbil.bi
#' }
#'       
#' @export
#' 
#' @importFrom stats fisher.test chisq.test ks.test quantile
#' @importFrom DescTools GTest
#' 
gof_gerbil <- function (x, y = NULL, type = 1, imp = 1, breaks = NULL, method = c("chi-squared", "fisher", "G"), ks = FALSE, partial = "imputed", ...) 
{

  if (length(imp) > 1) {
    warning("'imp' must be a scalar.")
  }

  if (tolower(type) == "univariate" | type == 1) {
    out <- gof_gerbil.uni(gerb = x, vars = y, imp = imp, breaks = breaks, method = method, ks = ks, ...)
  } else if (tolower(type) == "bivariate" | tolower(type) == "multivariate" | type == 2) {
    out <- gof_gerbil.mult(gerb = x, vars = y, imp = imp, breaks = breaks, method = method, partial = partial, ...)  
  } else {
    stop("Input 'type' not recognized.")
  }

    gof_gerbil_object <- out
    class(gof_gerbil_object) = "gof_gerbil"
    return(gof_gerbil_object)
}


#' Prints a \code{gof_gerbil} object. Printed output pertains to the goodness-of-fit tests that are applied in order to compare the distribution between observed and imputed cases for relevant variables or variable pairs. 
#'
#' @param x object of \code{gof_gerbil} class
#' @param ... additional parameters to be passed down to inner functions.
#'
#' @return The functions \code{print.gof_gerbil} and \code{summary.gof_gerbil} display information 
#'   about the \code{cor_gerbil} object. The output displayed includes: 
#'   1) the average test statistic value across all variables or variable pairs contained in the object,
#'   2) the average p-value of all goodness-of-fit tests contained within the object, and
#'   3) the number of tests that yieled a p-value of less than 0.05.  
#' 
#' @export
print.gof_gerbil = function(x, ...) {

  stats <- x[[1]]
  pvals <- x[[2]]
  tests <- x[[3]]

  if (class(tests)[1] == "matrix") {
    type <- "bivariate"
    vars <- colnames(x[[3]])
    tests[upper.tri(tests)] <- NA
    stats[upper.tri(stats)] <- NA
    pvals[upper.tri(pvals)] <- NA
  } else {
    type <- "univariate"
    vars <- names(x[[3]])
  }
  
  cat("\nSummary analysis of ", type, " tests calculated using ", length(vars), " variables.\n\n", sep = "")

  meths <- c("Chi-Squared", "Kolmogorov-Smirnov", "Fisher's Exact", "G-Test")

  for (i in 1:length(meths)) {
    n.tmp <- sum(tests == meths[i], na.rm = TRUE)
    if (class(tests)[1] == "matrix") {
      n.na <- sum(tests == meths[i] & lower.tri(tests) & is.na(stats), na.rm = TRUE)
    } else {
      n.na <- sum(tests == meths[i] & is.na(stats), na.rm = TRUE)
    }
    if (n.tmp > 0) {
      stat.avg <- mean(stats[tests == meths[i]], na.rm = TRUE)
      pval.avg <- mean(pvals[tests == meths[i]], na.rm = TRUE)
      sig.tot <-  sum((pvals <= .05)[tests == meths[[i]]], na.rm = TRUE)

      cat("A total of ", n.tmp - n.na, " tests were performed using the ", meths[i]," method.\n", sep = "")
      if (n.na > 0) {
        cat("(A ", meths[i]," statistic could not be calculated for ", n.na, " variable pairs.)\n", sep = "")
      }
      cat("    Average statistic value: ", round(stat.avg, 4), "\n", sep = "")
      cat("    Average p-value: ", round(pval.avg, 4), "\n", sep = "")
      cat("    Number of tests with p-value less than 0.05: ", sig.tot, "\n", sep = "")
      cat("\n")
    } 
  }
}


#' Summarises a \code{gerbil} object. Printed output pertains to the goodness-of-fit tests that are applied in order to compare the distribution between observed and imputed cases for relevant variables or variable pairs. 
#'
#' @param object An object of \code{gof_gerbil} class
#' @param ... Additional parameters to be passed down to inner functions.
#'
#' @return The functions \code{print.gof_gerbil} and \code{summary.gof_gerbil} display information 
#'   about the \code{cor_gerbil} object. The output displayed includes: 
#'   1) the average test statistic value across all variables or variable pairs contained in the object,
#'   2) the average p-value of all goodness-of-fit tests contained within the object, and
#'   3) the number of tests that yieled a p-value of less than 0.05.  
#' 
#' @export
#'
summary.gof_gerbil = function(object, ...) {
  print(object, ...)
}


gof_gerbil.uni <- function (gerb, vars = NULL, imp = 1, breaks = NULL, method = c("chi-squared", "fisher", "G"), ks = TRUE, ...) {

  method <- method[1]

  #if(method == "G") {
  #  library(DescTools)
  #}

  if (is.logical(imp)) {
    if (imp) {
      imp <- as.numeric(dimnames(gerb)[[3]])
    } else {
      imp <- 1
    }
  }
  if(length(imp) == 0) {
    imp <- 1
  }

  imp <- imp[1]

  if (length(vars) == 0) {
    #vars <- rownames(gerb[[3]])
    #vars <- rownames(gerb$summary)[gerb$summary[, 4] > 0]
    vars <- colnames(gerb$missing)[colSums(gerb$missing == 1 | gerb$missing == 2) > 0]
  } else if (is.numeric(vars)) {
    vars <- rownames(gerb[[3]])[vars]
  }

  types <- as.character(gerb[[3]][vars, "Variable.Type"])

  cont.vars <- vars[grepl("continuous", types, fixed = TRUE) | types == "semicont"]

  breaks1 <- list()

  if (length(cont.vars) > 0) {
    breaks1[1:length(cont.vars)] <- 4
    names(breaks1) <- cont.vars
    if (!is.list(breaks) & length(names(breaks)) == 0 & length(breaks) > 0) {
      breaks1[1:length(cont.vars)] <- list(breaks)
    } else if (!is.list(breaks) & length(names(breaks)) > 0 & length(breaks) > 0) {
      for (i in 1:length(breaks)) {
        breaks1[names(breaks)[i]] <- breaks[i]
      }
    } else if (length(breaks) > 0) {
      if (length(names(breaks)) == 0 & length(breaks) == length(cont.vars)) {
        breaks1[cont.vars] <- breaks
      } else if (length(names(breaks)) == 0) {
        stop("Input 'breaks' must be named.")
      } else {
        for (i in 1:length(breaks)) {
          breaks1[names(breaks)[i]] <- breaks[[i]]
        }
      }
    }

    names(breaks1) <- cont.vars

    for (i in 1:length(breaks1)) {
      if (length(breaks1[[i]]) == 1) {
        d1 <- gerb[[1]][[imp]][, names(breaks1)[i]]

        if(is.element(names(breaks1)[i], names(gerb$mass.final))) {
          mass <- gerb$mass.final[names(breaks1)[i]]
          d1[d1 == mass] = NA
        }

        d1 <- d1[!is.na(d1)]
 
        breaks1[[i]] <- quantile(d1, probs = (1:(breaks1[[i]] - 1))/breaks1[[i]])
        breaks1[[i]] <- c(-Inf, breaks1[[i]], Inf)
      }
      if(sum(duplicated(breaks1[[i]])) > 0) {
        warning(paste0("Variable ", names(breaks1)[i], " has duplicated cutpoints."))
        breaks1[[i]] <- breaks1[[i]][!duplicated(breaks1[[i]])]
      }
    }
    breaks <- breaks1
  }

  test <- stats <- p.vals <- list()

  for(i in 1:length(vars)) {

    type <- types[i]
    var <- vars[i]
    
    if (length(imp) > 1) {
      dat <- matrix(NA, NROW(gerb[[1]][[imp[1]]]), length(imp))
      for(j in 1:length(imp)) {
        dat[, j] <- gerb[[1]][[imp[j]]][, var]
        if (length(gerb$ineligibles) > 0) {
          dat[gerb$ineligibles[, var], j] <- NA
        }
      }
    } else {
      dat <- gerb[[1]][[imp]][, var]
    }
    obs <- gerb[[2]][, var]
    obs <- as.numeric(obs == 1 | obs == 2)

    cens <- gerb[[2]][, var] == 2
    if(sum(cens, na.rm = TRUE) > 0) {
      if(length(imp) > 1) {
        dat[cens, ] <- NA
      } else {
        dat[cens] <- NA
      }
      #warning(paste0("Censored values detected when plotting ", var, ": Censored cases will be removed."))
      message(paste0("Censored values detected when plotting ", var, ": Censored cases will be removed."), "\n", appendLF = FALSE)
    }
    
    if (sum(obs == 1, na.rm = TRUE) == 0) {
      cont <- FALSE
      message(paste0("Will not run univariate analyses for ", var, ": No imputed values for the variable."), "\n", appendLF = FALSE)
    } else {
      cont <- TRUE
    }
     
    if ((type == "binary" | type == "categorical" | type == "ordinal") & cont) {
      
      tab <- table(as.matrix(dat)[, 1], obs)
      colnames(tab)[colnames(tab) == "1"] <- as.character(imp[1])
      if (NCOL(dat) > 1) {
        tab.tmp <- tab
        tab <- matrix(NA, NROW(tab), 1 + NCOL(dat))
        rownames(tab) <- rownames(tab.tmp)
        colnames(tab) <- c("0", as.character(imp))
        tab[, colnames(tab.tmp)] <- tab.tmp
        for (j in 2:length(imp)) {
          tab[, as.character(imp[j])] <- table(dat[, j], obs)[, "1"]
        }
      }
    
      if (tolower(method) == "fisher") {
        out <- fisher.test(tab)
        test[[i]] <- "Fisher's Exact"
        stats[[i]] <- out$estimate
        p.vals[[i]] <- out$p.value
      } else if (tolower(method) == "g") {
        out <- DescTools::GTest(tab)
        test[[i]] <- "G-Test"
        stats[[i]] <- out$statistic
        p.vals[[i]] <- out$p.value
      } else {
        out <- chisq.test(tab)
        test[[i]] <- "Chi-Squared"
        stats[[i]] <- out$statistic
        p.vals[[i]] <- out$p.value
      }

      names(test)[i] <- names(stats)[i] <- names(p.vals)[i] <- var
      
    } else if (type == "semicont" & cont) {
      
      mass <- gerb$mass.final[var]

      pos <- as.numeric(as.matrix(dat)[, 1] != mass)
      tab <- table(pos, obs)

      if (ks) {

        dat[dat == mass] <- NA

        d1 <- dat[obs == 0]
        d1 <- d1[!is.na(d1)]

        d2 <- dat[obs == 1]
        d2 <- d2[!is.na(d2)]

        if (tolower(method) == "fisher") {
          t1 <- fisher.test(tab)
          s1 <- t1$estimate
          p1 <- t1$p.value
          t1 <- "Fisher's Exact"
        } else if (tolower(method) == "g") {
          t1 <- DescTools::GTest(tab)
          s1 <- t1$statistic
          p1 <- t1$p.value
          t1 <- "G-Test"
        } else {
          t1 <- chisq.test(tab)
          s1 <- t1$statistic
          p1 <- t1$p.value
          t1 <- "Chi-Squared"
        }

        t2 <- ks.test(d1, d2)
        test[[i]] <- matrix(c(t1, "Kolmogorov-Smirnov"), 1, 2)
        stats[[i]] <- matrix(c(s1, t2$statistic), 1, 2)
        p.vals[[i]] <- matrix(c(p1, t2$p.value), 1, 2)

        colnames(test[[i]]) <- colnames(stats[[i]]) <- colnames(p.vals[[i]]) <- c(paste0(var, ".B"), var)
        #colnames(test[[i]]) <- colnames(stats[[i]]) <- colnames(p.vals[[i]]) <- c(paste0("", ".B"),"")

        #names(test)[i] <- names(stats)[i] <- names(p.vals)[i] <- var

      } else {

        dat1 <- cut(dat, breaks = breaks[var][[1]])
        if (sum(is.na(dat1)) > sum(is.na(dat))) {
          warning(paste0("NAs generated when applying 'cut()' to ", var))
        }

        tab1 <- table(dat1, obs)
        tab <- rbind(tab["0",],tab1)
        rownames(tab)[1] <- as.character(mass)

        if (tolower(method) == "fisher") {
          t1 <- fisher.test(tab)
          s1 <- t1$estimate
          p1 <- t1$p.value
          t1 <- "Fisher's Exact"
        } else if (tolower(method) == "g") {
          t1 <- DescTools::GTest(tab)
          s1 <- t1$statistic
          p1 <- t1$p.value
          t1 <- "G-Test"
        } else {
          t1 <- chisq.test(tab)
          s1 <- t1$statistic
          p1 <- t1$p.value
          t1 <- "Chi-Squared"
        }

        test[[i]] <- c(t1)
        stats[[i]] <- c(s1)
        p.vals[[i]] <- c(p1)

        names(test)[i] <- names(stats)[i] <- names(p.vals)[i] <- var
        #c(paste0(var, ".B"), var)

      }
      
    } else if (cont) {
     
      if (ks) {

        d1 <- dat[obs == 0]
        d1 <- d1[!is.na(d1)]

        d2 <- dat[obs == 1]
        d2 <- d2[!is.na(d2)]

        out <- ks.test(d1, d2)
        test[[i]] <- "Kolmogorov-Smirnov"
        stats[[i]] <- out$statistic
        p.vals[[i]] <- out$p.value

      } else {
        num.na <- sum(is.na(dat))
        dat1 <- cut(dat, breaks = breaks[var][[1]])
        if (sum(is.na(dat)) > num.na) {
          warning(paste0("NAs generated when applying 'cut()' to ", var))
        }

        tab <- table(dat1, obs)

        if (tolower(method) == "fisher") {
          t1 <- fisher.test(tab)
          s1 <- t1$estimate
          p1 <- t1$p.value
          t1 <- "Fisher's Exact"
        } else if (tolower(method) == "g") { 
          t1 <- DescTools::GTest(tab)
          s1 <- t1$statistic
          p1 <- t1$p.value
          t1 <- "G-Test"
        } else {
          t1 <- chisq.test(tab)
          s1 <- t1$statistic
          p1 <- t1$p.value
          t1 <- "Chi-Squared"
        }

        test[[i]] <- c(t1)
        stats[[i]] <- c(s1)
        p.vals[[i]] <- c(p1)

      }

      names(test)[i] <- names(stats)[i] <- names(p.vals)[i] <- var

    }
  }

stats <- as.matrix(data.frame(stats))[1, ]
p.vals <- as.matrix(data.frame(p.vals))[1, ]
test <- as.matrix(data.frame(test))[1, ]

return(list(Stats = stats, p.values = p.vals, Test = test, Breaks = breaks))

}


gof_gerbil.mult <- function (gerb, vars = NULL, imp = 1, breaks = NULL, method = c("chi-squared", "fisher", "G"), partial = "imputed", ...) {

  method <- method[1]

  #if(method == "G") {
  #  library(DescTools)
  #}

  if (is.logical(imp)) {
    if (imp) {
      imp <- as.numeric(dimnames(gerb)[[3]])
    } else {
      imp <- 1
    }
  }
  if(length(imp) == 0) {
    imp <- 1
  }

  imp <- imp[1]

  if (length(vars) == 0) {
    #vars <- rownames(gerb[[3]])
    #vars <- rownames(gerb$summary)[gerb$summary[, 4] > 0]
    #vars <- colnames(gerb$missing)[colSums(gerb$missing == 1 | gerb$missing == 2) > 0]
    vars <- colnames(gerb$missing)
  } else if (is.numeric(vars)) {
    vars <- rownames(gerb[[3]])[vars]
  }

  if (length(vars) <= 1) {
      stop("Need at least 2 variables for bivariate plotting.")
  }

  types <- as.character(gerb[[3]][vars, "Variable.Type"])

  cont.vars <- vars[grepl("continuous", types, fixed = TRUE) | types == "semicont"]

  breaks1 <- list()

  if (length(cont.vars) > 0) {
    breaks1[1:length(cont.vars)] <- 4
    names(breaks1) <- cont.vars
    if (!is.list(breaks) & length(names(breaks)) == 0 & length(breaks) > 0) {
      breaks1[1:length(cont.vars)] <- list(breaks)
    } else if (!is.list(breaks) & length(names(breaks)) > 0 & length(breaks) > 0) {
      for (i in 1:length(breaks)) {
        breaks1[names(breaks)[i]] <- breaks[i]
      }
    } else if (length(breaks) > 0){
      if (length(names(breaks)) == 0 & length(breaks) == length(cont.vars)) {
        breaks1[cont.vars] <- breaks
      } else if (length(names(breaks)) == 0) {
        stop("Input 'breaks' must be named.")
      } else {
        for (i in 1:length(breaks)) {
          breaks1[names(breaks)[i]] <- breaks[[i]]
        }
      }
    }

    names(breaks1) <- cont.vars

    for (i in 1:length(breaks1)) {
      if (length(breaks1[[i]]) == 1) {
        d1 <- gerb[[1]][[imp]][, names(breaks1)[i]]

        if(is.element(names(breaks1)[i], names(gerb$mass.final))) {
          mass <- gerb$mass.final[names(breaks1)[i]]
          d1[d1 == mass] = NA
        }

        d1 <- d1[!is.na(d1)]
 
        breaks1[[i]] <- quantile(d1, probs = (1:(breaks1[[i]] - 1))/breaks1[[i]])
        breaks1[[i]] <- c(-Inf, breaks1[[i]], Inf)
      }
      if(sum(duplicated(breaks1[[i]])) > 0) {
        warning(paste0("Variable ", names(breaks1)[i], " has duplicated cutpoints."))
        breaks1[[i]] <- breaks1[[i]][!duplicated(breaks1[[i]])]
      }
    }
    breaks <- breaks1
  }

  dat <- gerb[[1]][[imp]][, vars, drop = FALSE]

  if (length(breaks) > 0) {
    for(i in 1:length(breaks)) {
      var <- names(breaks)[i]
      dat1 <- dat[, var]

      if (types[vars == var] == "semicont") {

        mass <- gerb$mass.final[var]
        pos <- as.numeric(dat1 != mass)
        dat1[pos == 0] <- NA
        num.na <- sum(is.na(dat1))
        dat1 <- cut(dat1, breaks = breaks[var][[1]])
        if (sum(is.na(dat1)) > num.na) {
          warning(paste0("NAs generated when applying 'cut()' to ", var))
        }
        dat1 <- as.character(dat1)
        dat1[pos == 0] <- "mass"
        dat1 <- factor(dat1)

      } else {
        num.na <- sum(is.na(dat1))
        dat1 <- cut(dat1, breaks = breaks[var][[1]])
        if (sum(is.na(dat1)) > num.na) {
          warning(paste0("NAs generated when applying 'cut()' to ", var))
        }
      }

      dat[, var] <- dat1
    }
  }

  test <- stats <- p.vals <- matrix(NA, length(vars), length(vars))
  dimnames(test) <- dimnames(stats) <- dimnames(p.vals) <- list(vars, vars)

  for(i in 2:length(vars)) {
    for(j in 1:(i - 1)) {

      var1 <- vars[i]
      var2 <- vars[j]

      dat1 <- dat[, var1]
      dat2 <- dat[, var2]
    
      if (partial == "imputed") {
        miss <- gerb[[2]][, var1] == 1 | gerb[[2]][, var2] == 1
      } else {
        miss <- gerb[[2]][, var1] == 1 & gerb[[2]][, var2] == 1
      }

      cens1 <- gerb[[2]][, var1] == 2
      if(sum(cens1, na.rm = TRUE) > 0) {
        if(length(imp) > 1) {
          dat1[cens1, ] <- NA
        } else {
          dat1[cens1] <- NA
        }
        message(paste0("Censored values detected when plotting ", var1, ": Censored cases will be removed."), "\n", appendLF = FALSE)
      }

      cens2 <- gerb[[2]][, var2] == 2
      if(sum(cens2, na.rm = TRUE) > 0) {
        if(length(imp) > 1) {
          dat1[cens2, ] <- NA
        } else {
          dat1[cens2] <- NA
        }
        message(paste0("Censored values detected when plotting ", var2, ": Censored cases will be removed."), "\n", appendLF = FALSE)
      }
    
      if (sum(miss, na.rm = TRUE) == 0) {
        cont <- FALSE
        message(paste0("Will not perform bivariate test between ", 
           var1, " and ", var2, ": No imputed values for either variable."), 
          "\n", appendLF = FALSE)
      } else {
        cont <- TRUE
      }

      dat3 <- interaction(factor(dat1), factor(dat2))

      tab <- table(dat3, miss)

      if(NROW(tab) == 1 | NCOL(tab) == 1 | sum(tab) == 0) {
        test[i, j] <- NA
        stats[i, j] <- NA
        p.vals[i, j] <- NA
      } else if (tolower(method) == "fisher") {
        out <- fisher.test(tab)
        test[i, j] <- test[j, i] <- "Fisher's Exact"
        stats[i, j] <- stats[j, i] <- out$estimate
        p.vals[i, j] <- p.vals[j, i] <- out$p.value
      } else if (tolower(method) == "g") {
        out <- DescTools::GTest(tab)
        stats[i, j] <- stats[j, i] <- out$statistic
        p.vals[i, j] <- p.vals[j, i] <- out$p.value
        test[i, j] <- test[j, i] <- "G-Test"
      } else {
        out <- chisq.test(tab)
        test[i, j] <- test[j, i] <- "Chi-Squared"
        stats[i, j] <- stats[j, i] <- out$statistic
        p.vals[i, j] <- p.vals[j, i] <- out$p.value
      }
    }
  }

return(list(Stats = stats, p.values = p.vals, Test = test, Breaks = breaks))

}




