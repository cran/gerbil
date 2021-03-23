#' @title Correlation Analysis for \code{gerbil} Objects
#' 
#' @description
#' 
#' This function assesses the bivariate properties of imputed data using a correlation analysis.  
#' Specifically, it calculates pairwise correlations for observed cases and for imputed cases. 
#' The function also calculates the Fisher z-transformation for 
#' each correlation and performs a hypothesis test using the transformed correlations in order to 
#' compare correlations calculated using imputed cases to those calculated using observed cases.  
#' 
#' @details
#'
#' Cases are assigned a status of being observed or imputed in a pairwise fashion. That is, a specific 
#' data unit may be considered observed when calculating a correlation for one pair of variables and be
#' imputed when calculating a correlation for another pair. For a given pair of variables, cases that 
#' have both variables observed are always treated as observed, and cases that have both variables missing 
#' are always treated as imputed. Cases that have only one variable in the pair observed (i.e., those that are
#' partially imputed) are treated as imputed when the input \code{partial = 'imputed'} (the default) and are 
#' otherwise treated as observed.  
#' 
#' Correlations are calculated across an expanded dataset that creates binary indicators for categorical variables and for semicontinuous variables.
#' Unlike the algorithm used to calculate the imputations, missingness is not artificially imposed in any binary indicator.  
#' Missingness is imposed, however, in the variable corresponding to the continuous portion of a semicontinuous variable. 
#' 
#' Note that the hypothesis test based upon the Fisher z-transformation is based off of bivariate normal assumptions. 
#' As such, p-values may be misleading in data where this assumption does not hold.  
#' 
#' @param x A \code{gerbil} object containing the imputed data. 
#' @param y A vector listing the column names of the imputed data that will be included in the correlation analysis. By default, \code{y} contains all columns of the data that required imputation. If \code{TRUE}, all variables with missing values eligible for imputation are used.
#' @param imp A scalar indicating which of the multiply imputed datasets should be used for the analysis.  Defaults to \code{imp = 1}. 
#' @param log A character vector that includes names variable of which a log transformation is to be taken prior to calculating correlations.
#' @param partial Indicates how partially imputed pairs are handled when calculating correlations. If \code{partial = 'imputed'}, cases with at least one missing variable in a pair are considered imputed. Otherwise (\code{partial = 'observed'}), only cases with both variables in the pair missing are considered imputed.
#' 
#' @return \code{cor_gerbil()} retuns an object of the class \code{cor_gerbil} that has following slots:
#' 
#' \describe{
#'         \item{Correlations}{A list containing two elements -- these are named \code{Observed}, \code{Imputed}, and \code{All}. The first is a matrix giving the sample correlations when calculated across cases labeled as observed. The second and third are analogous correlation matrices calculated across only cases labeled as imputed and across all cases, respectively.}
#'         \item{n}{A list containing two elements -- these are named \code{Observed}, \code{Imputed}, and \code{All}. The first is a matrix giving number of cases in the respective pair of variables that have been labeled as observed. The second and third are analogous matrices indicating the number of cases labeled as imputed for each pair and indicating the total number of cases for each pair, respectively.}
#'         \item{Fisher.Z}{A list containing two elements -- these are named \code{Observed}, \code{Imputed}, and \code{All}. These matrices give the Fisher z-transformation of the correlations in the matrices provided in the slot \code{Correlations}.} 
#'         \item{Statistic}{A matrix that gives the value of the test statistic based on the Fisher z-transformation for each pair of variables. This statistic may be used to assess whether the correlations calculated across cases labeled as observed are statistically different from the correlations calculated across cases labeled as imputed.}
#'         \item{p.value}{A matrix that list the p-value for each test statistic provided in the matrix in the slot labeled \code{Statistic}.}
#' }
#' 
#' @export
#' 
#' @importFrom stats pnorm
#' 
#' @examples
#' \donttest{
#' #Load the India Human Development Survey-II dataset
#' data(ihd_mcar) 
#' 
#' imps.gerbil <- gerbil(ihd_mcar, m = 1, mcmciter = 100, ords = "education_level", 
#'        semi = "farm_labour_days", bincat = c("sex", "marital_status", "job_field", "own_livestock"))
#' 
#' #Run the correlation analysis
#' cors.gerbil <- cor_gerbil(imps.gerbil, imp = 1)
#' 
#' #Print a summary
#' cors.gerbil
#' }
cor_gerbil <- function (x, y = NULL, imp = 1, log = NULL, partial = "imputed") 
{
    bin.meth <- "binary"
    cat.meth <- "categorical"
    ord.meth <- "ordinal"
    cens.meth <- "censored"
    semi.meth <- "semicont"
    cont.meth <- "EMP"

    vars <- y

    if (is.logical(vars)) {
        if (vars) {
            vars <- colnames(x$missing.latent)[colSums(x$missing.latent == 
                1 | x$missing.latent == 2, na.rm = TRUE) > 0]
        }
        else {
            vars <- colnames(x$missing.latent)
        }
    }
    else if (length(vars) == 0) {
        vars <- colnames(x$missing.latent)
    }

    latent.vars <- colnames(x$missing.latent)
    orig.vars <- dimnames(x[[1]][[1]])[[2]]
    nams.out <- x$nams.out
    vars.old <- vars.new <- NULL

    for (i in 1:length(vars)) {
        if (is.element(vars[i], latent.vars)) {
            tmp <- vars[i]
            tmp1 <- nams.out[names(nams.out) == vars[i]]
        }
        else if (is.element(vars[i], orig.vars)) {
            tmp <- names(nams.out)[nams.out == vars[i]]
            tmp1 <- vars[i]
        }
        else {
            message(paste0("Variable name ", vars[i], " is invalid and will not be used for correlation analysis."), 
                "\n", appendLF = FALSE)
            tmp <- tmp1 <- NULL
        }
        if (length(vars.new) == 0 & length(tmp) > 0) {
            vars.new <- tmp
            vars.old <- tmp1
        }
        else if (length(tmp) > 0) {
            vars.new[(length(vars.new) + 1):(length(vars.new) + 
                length(tmp))] <- tmp
            vars.old[(length(vars.old) + 1):(length(vars.old) + 
                length(tmp1))] <- tmp1
        }
    }
    if (length(vars.new) == 0) {
        stop("No valid variable names were given.")
    }
    vars.old <- vars.old[!duplicated(vars.old)]
    vars <- vars.new
    dat1 <- x$imputed[[imp]]
    dat1 <- dat1[, vars.old, drop = FALSE]
    if (length(x$ineligibles) > 0) {
        dat1[x$ineligibles[, vars.old]] <- NA
    }
    type <- x$summary[, "Variable.Type"]
    names(type) <- rownames(x$summary)
    type <- type[vars.old]
    type[grepl(cont.meth, type)] <- cont.meth
    predictorMatrix <- x$predMat.initial
    predictorMatrix <- predictorMatrix[vars.old, vars.old, drop = FALSE]
    mass <- x$mass.final
    if (length(mass) > 0) {
        if (length(intersect(names(mass), vars.old)) > 0) {
            mass <- mass[intersect(names(mass), vars.old)]
        }
        else {
            mass <- NULL
        }
    }
    predictorMatrix[upper.tri(predictorMatrix)] <- 0
    non.cont <- which(type != cont.meth & type != cens.meth & 
        type != semi.meth)
    categories <- list()
    if (length(non.cont) > 0) {
        for (j in 1:length(non.cont)) {
            i <- non.cont[j]
            dat1[, i] <- factor(dat1[, i])
            categories[[j]] <- names(table(dat1[, i]))
            dat1[, i] <- as.numeric(dat1[, i]) - 1
        }
        names(categories) <- colnames(dat1)[non.cont]
    }
    catcols <- which(type == cat.meth)
    semicols <- which(type == semi.meth)
    if (length(catcols) > 0 | length(semicols) > 0) {
        tmpdat <- make.cat(dat1, predictorMatrix, catcols = catcols, 
            semicols = semicols, type = type, bin.meth = bin.meth, 
            cont.meth = cont.meth, mass = mass, categories = categories, 
            reord = TRUE, impose.miss = FALSE)
        predictorMatrix <- tmpdat[[3]]
        cat.nams <- tmpdat[[4]]
        cat.cols <- tmpdat[[5]]
        is.cat <- tmpdat[[6]]
        type <- tmpdat[[7]]
        semicols.new <- tmpdat[[8]]
        dat1 <- tmpdat[[1]]
    }
    else {
        cat.nams <- cat.cols <- is.cat <- semicols.new <- NULL
    }
    dat1 <- dat1[, vars.new, drop = FALSE]
    mis <- x$missing.latent
    mis <- mis[, colnames(dat1)]
    dat1 <- as.matrix(dat1)
    if (length(log) > 0) {
        for (i in 1:length(log)) {
            if (is.element(log[i], colnames(dat1))) {
                if (sum(dat1[, log[i]] <= 0, na.rm = TRUE) > 
                  0) {
                  dat1[dat1[, log[i]] <= 0, log[i]] <- NA
                  message(paste0("Non-positive values detected for ", 
                    log[i], ": Cases will be removed prior to log transformation."), 
                    "\n", appendLF = FALSE)
                }
                dat1[, log[i]] <- log(dat1[, log[i]])
                colnames(dat1)[colnames(dat1) == log[i]] <- paste0("log(", 
                  log[i], ")")
            }
            else {
                warning(paste0("Variable ", log[i], " does not appear in the dataset. No log will be taken."))
            }
        }
    }
    is.obs <- !is.na(dat1) & (mis != 1 & mis != 2)
    is.mis <- !is.na(dat1) & (mis == 1 | mis == 2)
    is.all <- is.obs | is.mis
    dat1[is.na(dat1)] <- 0
    dat12 <- dat1^2
    SSxy.all <- crossprod(dat1, dat1)
    SSxx.all <- crossprod(is.all, dat12)
    SSyy.all <- crossprod(dat12, is.all)
    SSx.all <- crossprod(is.all, dat1)
    SSy.all <- crossprod(dat1, is.all)
    n.all <- crossprod(is.all, is.all)
    if (partial == "imputed") {
        SSxy.obs <- crossprod(dat1 * is.obs, dat1 * is.obs)
        SSxx.obs <- crossprod(is.obs, dat12 * is.obs)
        SSyy.obs <- crossprod(dat12 * is.obs, is.obs)
        SSx.obs <- crossprod(is.obs, dat1 * is.obs)
        SSy.obs <- crossprod(dat1 * is.obs, is.obs)
        n.obs <- crossprod(is.obs, is.obs)
        SSxy.mis <- SSxy.all - SSxy.obs
        SSxx.mis <- SSxx.all - SSxx.obs
        SSyy.mis <- SSyy.all - SSyy.obs
        SSx.mis <- SSx.all - SSx.obs
        SSy.mis <- SSy.all - SSy.obs
        n.mis <- n.all - n.obs
    }
    else {
        SSxy.mis <- crossprod(dat1 * is.mis, dat1 * is.mis)
        SSxx.mis <- crossprod(is.mis, dat12 * is.mis)
        SSyy.mis <- crossprod(dat12 * is.mis, is.mis)
        SSx.mis <- crossprod(is.mis, dat1 * is.mis)
        SSy.mis <- crossprod(dat1 * is.mis, is.mis)
        n.mis <- crossprod(is.mis, is.mis)
        SSxy.obs <- SSxy.all - SSxy.mis
        SSxx.obs <- SSxx.all - SSxx.mis
        SSyy.obs <- SSyy.all - SSyy.mis
        SSx.obs <- SSx.all - SSx.mis
        SSy.obs <- SSy.all - SSy.mis
        n.obs <- n.all - n.mis
    }
    Sigma.xy.all <- (n.all * SSxy.all - SSx.all * SSy.all)/n.all
    Sigma.xx.all <- (n.all * SSxx.all - SSx.all * SSx.all)/n.all
    Sigma.yy.all <- (n.all * SSyy.all - SSy.all * SSy.all)/n.all
    R.all <- Sigma.xy.all/sqrt(Sigma.xx.all * Sigma.yy.all)
    R.all[n.all <= 3] <- NA
    R.all[abs(R.all) > 1] <- NA
    Sigma.xy.obs <- (n.obs * SSxy.obs - SSx.obs * SSy.obs)/n.obs
    Sigma.xx.obs <- (n.obs * SSxx.obs - SSx.obs * SSx.obs)/n.obs
    Sigma.yy.obs <- (n.obs * SSyy.obs - SSy.obs * SSy.obs)/n.obs
    R.obs <- Sigma.xy.obs/sqrt(Sigma.xx.obs * Sigma.yy.obs)
    R.obs[n.obs <= 3] <- NA
    R.obs[abs(R.obs) > 1] <- NA
    Sigma.xy.mis <- (n.mis * SSxy.mis - SSx.mis * SSy.mis)/n.mis
    Sigma.xx.mis <- (n.mis * SSxx.mis - SSx.mis * SSx.mis)/n.mis
    Sigma.yy.mis <- (n.mis * SSyy.mis - SSy.mis * SSy.mis)/n.mis
    R.mis <- Sigma.xy.mis/sqrt(Sigma.xx.mis * Sigma.yy.mis)
    R.mis[n.mis <= 3] <- NA
    R.mis[abs(R.mis) > 1] <- NA
    R.all1 <- R.all
    R.obs1 <- R.obs
    R.mis1 <- R.mis
    R.all1[abs(R.all) == 1] <- NA
    R.obs1[abs(R.obs) == 1] <- NA
    R.mis1[abs(R.mis) == 1] <- NA
    suppressWarnings(z.all <- 0.5 * log((1 + R.all1)/(1 - R.all1)))
    suppressWarnings(z.obs <- 0.5 * log((1 + R.obs1)/(1 - R.obs1)))
    suppressWarnings(z.mis <- 0.5 * log((1 + R.mis1)/(1 - R.mis1)))
    suppressWarnings(se.diff.r <- sqrt(1/(n.obs - 3) + 1/(n.mis - 
        3)))
    diff <- z.obs - z.mis
    z.stat <- abs(diff/se.diff.r)
    p.val <- (1 - pnorm(z.stat))
    twotailed <- TRUE
    if (twotailed) {
        p.val <- 2 * p.val
    }
    cor_gerbil_object <- list(Correlations = list(Observed = R.obs, 
        Imputed = R.mis, All = R.all), n = list(Observed = n.obs, 
        Imputed = n.mis, All = n.all), Fisher.Z = list(Observed = z.obs, 
        Imputed = z.mis, All = z.all), Statistic = z.stat, p.value = p.val)
    class(cor_gerbil_object) = "cor_gerbil"
    return(cor_gerbil_object)
}


#' Prints a \code{cor_gerbil} object. Printed output includes the average difference of correlations, as well as summaries of the test statistics based on Fisher's z and their p-values. 
#'
#' @param x object of \code{cor_gerbil} class
#' @param ... additional parameters to be passed down to inner functions.
#'
#' @return The functions \code{print.cor_gerbil} and \code{summary.cor_gerbil} display information 
#'   about the \code{cor_gerbil} object. The output displayed includes: 
#'   1) the average absolute difference in correlation between observed and imputed cases across all relevant variable pairs,
#'   2) the average value of the test statistic based on Fisher's z across all variable pairs,
#'   3) the largest test statistic observed across any variable pair, and
#'   4) the portion of p-values for the test based on Fisher's z that are less than 0.05.
#' 
#' @export
print.cor_gerbil = function(x, ...) {

  cors.obs <- x$Correlations$Observed
  cors.mis <- x$Correlations$Imputed
  cors.stats <- x$Statistic
  cors.p.vals <- x$p.value

  cor.diffs <- abs(cors.obs - cors.mis)
  diag(cor.diffs) <- NA
  
  cat("\nSummary analysis comparing correlations calculated between observed cases with \n   corresponding correlations calcluated between imputed cases:\n\n")

  cat("Average absolute difference in correlation between observed and imputed cases:\n\n")
  cat(round(mean(cor.diffs[upper.tri(cor.diffs)], na.rm = TRUE), 4))
  cat("\n\n")

  cors.stats[lower.tri(cors.stats)] <- 0
  cors.stats[diag(cors.stats)] <- 0

  colmn <- which(cors.stats == max(cors.stats, na.rm = TRUE)) %/% nrow(cors.stats) + 1
  row <- which(cors.stats == max(cors.stats, na.rm = TRUE)) %% nrow(cors.stats)

  cat("Average value of the test statistic based on Fisher's z across all variable pairs:\n\n")
  cat(round(mean(cors.stats[upper.tri(cors.stats)], na.rm = TRUE), 4))
  cat("\n\n")

  cat("The largest statistic (", round(cors.stats[row, colmn], 2), ") corresponds to variable pair ", rownames(cors.stats)[row], " and ", colnames(cors.stats)[colmn], ".\n\n", sep = "")

  cat("Portion of p-values for the test based on Fisher's z that are less than 0.05:\n\n", sep = "")
  cat(round(mean(cors.p.vals[upper.tri(cors.p.vals)] <= 0.05, na.rm = TRUE), 4))
  cat("\n\n")

}


#' Summarises a \code{gerbil} object. Printed output includes a variable-by-variable summary of variable types and missingness rates. The implemented predictor matrix is also provided.
#'
#' @param object An object of \code{cor_gerbil} class
#' @param ... additional parameters to be passed down to inner functions.
#'
#' @return The functions \code{print.cor_gerbil} and \code{summary.cor_gerbil} display information 
#'   about the \code{cor_gerbil} object. The output displayed includes: 
#'   1) the average absolute difference in correlation between observed and imputed cases across all relevant variable pairs,
#'   2) the average value of the test statistic based on Fisher's z across all variable pairs,
#'   3) the largest test statistic observed across any variable pair, and
#'   4) the portion of p-values for the test based on Fisher's z that are less than 0.05.
#' 
#' @export
#'
summary.cor_gerbil = function(object, ...) {
  print(object, ...)
}
