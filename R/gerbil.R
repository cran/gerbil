#### Gerbil Class Functions ----------------------------------------------------------####

#----------------------------------------------------------------------------------------#
# Author: Michael Robbins
# Purpose: This File contains the gerbil function itself and all functions used by the
# gerbil function.
# The gerbil function is the only function exported to package users and creates the
# gerbil S3 object, wich includes all results from the multiple imputation process.
# Creation Date: Sept 2020
#----------------------------------------------------------------------------------------#

#' @title General Efficient Regression-Based Imputation with Latent processes
#'
#' @description
#'
#' Coherent multiple imputation of general multivariate data as implemented through the GERBIL algorithm described by Robbins (2020).  The algorithm is
#'
#' \itemize{
#' \item \strong{coherent} in that imputations are sampled from a valid joint distribution, ensuring MCMC convergence;
#' \item \strong{general} in that data of general structure (binary, categorical, continuous, etc.) may be allowed;
#' \item \strong{efficient} in that computational performance is optimized using the SWEEP operator for both modeling and sampling;
#' \item \strong{regression-based} in that the joint distribution is built through a sequence of conditional regression models;
#' \item \strong{latent} in that a latent multivariate normal process underpins all variables; and
#' \item \strong{flexible} in that the user may specify which dependencies are enabled within the conditional models.
#' }
#'
#'
#' @details
#'
#' \code{gerbil} is designed to handle the following classes of variables:
#'
#' \itemize{
#'
#' \item \code{'continuous'}: Variables are transformed to be (nearly) standard normal prior to imputation. The default transformation method is based on empirical distributions (see Robbins, 2014) and ensures that imputed values of a variable are sampled from the observed values of that variable.
#' \item \code{'binary'}: Dichotomous variables are handled through probit-type models in that they are underpinned by a unit-variance normally distributed random variable.
#' \item \code{'categorical'}: Unordered categorical variables are handled by creating nested binary variables that underpin the categorical data.  Missingness is artificially imposed in the nested variables in order to ensure conditional independence between them.  See Robbins (2020) for details.
#' \item \code{'ordinal'}: Ordered categorical variables (ordinal) are handled through a probit-type model in that a latent normal distribution is assumed to underpin the ordinal observations.  See Robbins (2020) for details.
#' \item \code{'semicont'}: Mixed discrete/continuous (semi-continuous) variables are assumed to observe a mass at a specific value (most often zero) and are continuous otherwise.  A binary variable is created that indicates whether the semi-continuous variable takes on the point-mass value; the continuous portion is set as missing when the observed semi-continuous variable takes on the value at the point-mass.  See Robbins et al. (2013) for details.
#' }
#'
#' The parameter \code{type} allows the user to specify the class for each variable.  Routines are in place to establish the class by default for variables not stated in \code{type}. Note that it is not currently possible for a variable to be assigned a class of semi-continuous by default.
#'
#' \code{gerbil} uses a joint modeling approach to imputation that builds a joint model using a sequence of conditional models, as outlined in Robbins et al. (2013).
#' This approach differs from fully conditional specification in that the regression model for any given variable is only allowed to depend upon variables that preceed it in an index ordering.
#' The order is established by the parameter \code{visitSeq}. \code{gerbil} contains the flexibility to allow its user to establish which of the permissible dependencies are enabled within the conditional models.
#' Enabled dependencies are stated within the parameter \code{predMat}.  Note that the data matrix used for imputation is an expanded version of the data that are fed into the algorithm (variables are created that underpin unordered categorical and semi-continuous variables).
#' Note also that conditional dependencies between the nested binary variables of a single undordered categorical variables or the discrete and continuous portions of a semi-continuous variable are not permitted.
#'
#' The output of \code{gerbil} is an object of class \code{gerbil} which is a list that contains the imputed datasets (\code{imputed}), missingness indicators (\code{missing} and \code{missing.latent}), summary information (\code{summary}), output used for MCMC convergence diagostics (\code{chainSeq} and \code{R.hat}),
#' and modeling summaries (\code{visitSeq.initial}, \code{visitSeq.final}, \code{predMat.initial}, \code{predMat.final}, \code{drops}, and \code{forms}).
#' Some output regarding convergence diagnostics and modeling regards the expanded dataset used for imputation (the expanded dataset includes binary indicators for unordered categorical and semi-continuous variables).
#' Note that the nested binary variables corresponding to an unordered categorical variable \code{X} with categories labeled \code{a}, \code{b}, \code{c}, etc., are named \code{X.a}, \code{X.b}, \code{X.c}, and so forth in the expanded dataset.
#' Likewise, the binary variable indicating the point mass of a semi-continuous variable \code{Y} is named \code{Y.B} in the expanded dataset, and the positive portion (with missingness imposed) is left as being named \code{Y}.
#'
#' \code{gerbil} automatically checks each regression model for perfect collinearities and reduces the model as needed. 
#' Variables that have been dropped from a given model are listed in the element named \code{'drops'} in a \code{gerbil} object.  
#' 
#' @param dat The dataset that is to be imputed.  Missing values must be coded with \code{NA}. 
#' @param m The number of multiply imputed datasets to be created.  By default, \code{m = 1}.
#' @param mcmciter The number of iterations of Markov chain Monte Carlo that will be used to create each imputed dataset. By default, \code{m = 25}.
#' @param predMat A numeric matrix of \code{ncol(dat)} columns and no more than \code{nrow(dat)} rows, containing 0/1 data specifying the set of predictors to be used for each target row. Each row corresponds to a variable. A value of 1 means that the column variable is used as a predictor for the variable in the target row. By default, \code{predMat} is a square matrix of \code{ncol(dat)} rows and columns with 1's below the diagonal and 0's on and above the diagonal.  Any non-zero value on or above the diagonal will be set to zero.
#' @param type A named vector that gives the type of each variable contained in \code{dat}.  Possible types include \code{'binary'}, \code{'categorical'}, \code{'ordinal'}, \code{'semicont'} (semi-continuous), and \code{'continuous'}.  The vector type should be named where the names indicate the corresponding column of \code{dat}.  Types for variables not listed in type will be determined by default, in which case a variable with no more than \code{num.cat} possible values will be set as binary/categorical and is set as continuous otherwise.
#' @param visitSeq A vector of variable names that has (at least) contains all names of each column of \code{dat} that has missing values.  Within the I-Step and P-Step of gerbil, the variables will be modeled and imputed in the sequence given by \code{visitSeq}.   If \code{visitSeq = TRUE}, \code{visitSeq} is reset as being equal to the columns of \code{dat} ordered from least to most missingness.  If \code{visitSeq = NULL} (the default) or \code{visitSeq = FALSE} variables are ordered in accordance with the order of the rows of \code{predMat} or (if unavailable) the order in which they appear in the \code{dat}.
#' @param ords A character string giving a set of the column names of \code{dat} that indicate which variables are to be treated as ordinal. Elements of \code{ords} are overridden by any conflicting information in \code{type}. By default, \code{ords = NULL}.
#' @param semi A character string giving a set of the column names of \code{dat} that indicate which variables are to be treated as semi-continuous. Elements of \code{semi} are overridden by any conflicting information in \code{type}. By default, \code{semi = NULL}.
#' @param bincat  A character string giving a set of the column names of \code{dat} that indicate which variables are to be treated as binary or unordered categorical. Elements of \code{bincat} are overridden by any conflicting information in \code{type}. By default, \code{bincat = NULL}.
#' @param cont.meth The type of marginal transformation used for continuous variables.  Set to \code{"EMP"} by default for the empirical distribution transformation of Robbins (2014). The current version also includes an option for no transformation (\code{cont.meth = "none"}). Other transformation types will be available in future versions of \code{gerbil}.  .
#' @param num.cat Any variable that does not have a type specified by any of the other parameters will be treated as categorical if it takes on no more than \code{num.cat} possible values and as continuous if it takes on more than \code{num.cat} possible values. By default, \code{num.cat = 12}.
#' @param r The number of pairwise completely observed cases that must be available for any pair of variables to have dependencies enabled within the conditional models for imputation. By default, \code{r = 5}.
#' @param verbose If \code{TRUE} (the default), history is printed on console. Use \code{verbose = FALSE} for silent computation.
#' @param n.cores The number of CPU cores to use for parallelization. If \code{n.cores} is not specified by the user, it is guessed using the \code{detectCores} function in the parallel package.  If \code{TRUE}  (the default), it is set as \code{detectCores()}.  If \code{NULL}, it is set as \code{floor((detectCores()+1)/2)}.  If \code{FALSE}, it is set as \code{1}, in which case parallelization is not invoked.  Note that the documentation for \code{detectCores} makes clear that it is not failsafe and could return a spurious number of available cores. By default, \code{n.cores} is set as \code{floor((n + 1)/2)}, where \code{n} is the number of available clusters.
#' @param cl.type The cluster type that is passed into the \code{makeCluster()} function in the \code{parallel} package.  Defaults to \code{'PSOCK'}.
#' @param mass A named vector of the same length as the number of semi-continuous variables in \code{dat} that gives the location (value) of the point mass for each such variable. The point of mass for each semicontinuous variable is set to zero by default. 
#' @param ineligible Either a scalar or a matrix that is used to determined which values are to be considered missing but ineligible for imputation. Such values will be imputed internally within \code{gerbil} to ensure a coherent imputation model but will be reset as missing after imputations have been created. If \code{ineligible} is a scalar, all data points that take on the respective value will be considered missing but ineligible for imputation. If \code{ineligible} is a matrix (with the same number of rows as \code{dat} and column names that overlap with \code{dat}), entries of \code{TRUE} or \code{1} in \code{ineligible} indicate values that are missing but ineligible for imputation. If \code{ineligible = NULL} (the default), all missing values will be considered eligible for imputation. 
#' @param trace A logical that, if \code{TRUE}, implies that means and variances of variables are tracked across iterations. Set to \code{FALSE} to save computation time. However, trace plots and R hat statistics are disabled for \code{gerbil} objects created with \code{trace = FALSE}. Defaults to \code{TRUE}.
#' @param seed An integer that, when specified, is used to set the random number generator via \code{set.seed()}. 
#' 
#' @return \code{gerbil()} returns an object the class \code{gerbil} that contains the following slots:
#'
#'\describe{
#'        \item{imputed}{A list of length \code{m} that contains the imputed datasets.}
#'        \item{missing}{A matrix \code{0}s, \code{1}s, \code{2}s, and \code{4}s of the same dimension as \code{dat} that indicates which values were observed or missing.  A \code{0} indicates a fully observed value, a \code{1} indicates a missing value that was imputed, and a \code{4} indicates a missing value that was ineligible for imputation.}
#'        \item{summary}{A matrix with \code{ncol(dat)} number of rows that contains summary information, including the type of each variable and missingness rates. Note that for continuous variables, the type listed indicates the method of transformation used.}
#'        \item{chainSeq}{A list of six elements. Each element is a matrix with \code{mcmciter} columns and up to \code{ncol(dat)} rows. Objects \code{means.all} and \code{means.mis} give the variables means of data process across iterations of MCMC when all observations are incorporated and when only imputed values are incorporated, respectively. (Means of continuous variables are given on the transformed scale.) Similar objects are provided to track variances of variables. Variables are listed in the order provided by the \code{gerbil} object \code{visitSeq.latent}. Variables reported in this output are those contained in the dataset that has been expanded to include binary indicators for categorical and semi-continuous variables.}
#'        \item{R.hat}{The value of the R hat statistics of Gelman and Rubin (1992) for the means and variances of each variable. The R hat statistic is also provided for mean of binary variables.  Variables include those contained in the expanded dataset and are listed in the order provided by object \code{visitSeq.latent}. Only calculated if \code{m > 2} and \code{mcmciter >= 4}.}
#'        \item{missing.latent}{A matrix of the same dimensions as the expanded dataset, but used to indicate missingness in the expanded dataset. In this matrix, \code{0}s indicate fully observed values, \code{1}s indicate fully missing values, \code{3}s indicate values that have imposed missingness (for binary indicators corresponding to categorical or semi-continuous variables), and \code{4} indicates a missing value that is ineligible for imputation (as determined by the input \code{'ineligible'})..}
#'        \item{visitSeq.initial}{A vector of variable names giving the sequential ordering of variables that is used for imputation prior to expanding the dataset include nested binary and point-mass indicators.  Variables without missing values are excluded.}
#'        \item{visitSeq.final}{A vector of variable names giving the sequential ordering of variables in the expanded dataset that is used for imputation.  Variables without missing values are excluded.}
#'        \item{predMat.initial}{A matrix of ones and zeros indicating the dependencies enabled in the conditional models used for imputation.  This matrix is determined from the input 'predMat'. Rows corresponding to variables with no missing values are removed.}
#'        \item{predMat.final}{A matrix of ones and zeros indicating the dependencies enabled in the conditional models used for imputation.  This is of a similar format to the input 'predMat' but pertains to the expanded dataset.  Rows corresponding to variables with no missing values are removed.}
#'        \item{drops}{A list of length equal to the number of variables in the expanded dataset that have missing values.  Elements of the list indicate which variables were dropped from the conditional model for the corresponding variable due to either insufficient pairwise complete observations (see the input 'r') or perfect collinearities.}
#'        \item{forms}{A list of length equal to the number of variables in the expanded dataset that have missing values.  Elements of the list indicate the regression formula used for imputation of the respective variable.}
#'        \item{mass.final}{The final version of the input parameter \code{mass}.}
#'        \item{ineligibles}{A logical matrix with the same number of rows and columns as \code{dat} that indicates which elements are considered missing but ineligible for imputation.}
#'        \item{nams.out}{A vector used to link column names in the expanded data to corresponding names in the original data.}
#'}
#'
#' @references Gelman, A., & Rubin, D. B. (1992). Inference from iterative simulation using multiple sequences. \emph{Statistical Science}, 7(4), 457-472.
#'
#'   Robbins, M. W. (2014). The Utility of Nonparametric Transformations for Imputation of Survey Data. \emph{Journal of Official Statistics}, 30(4), 675-700.
#' 
#'   Robbins, M. W. (2020). A flexible and efficient algorithm for joint imputation of general data. arXiv preprint arXiv:2008.02243.
#'
#'   Robbins, M. W., Ghosh, S. K., & Habiger, J. D. (2013). Imputation in high-dimensional economic data as applied to the Agricultural Resource Management Survey. \emph{Journal of the American Statistical Association}, 108(501), 81-95.
#'
#' @export
#'
#' @importFrom stats cov lm model.matrix pnorm predict qnorm rchisq rnorm sd var
#' @importFrom utils data
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom pbapply pblapply
#' @importFrom truncnorm rtruncnorm
#' @importFrom mvtnorm rmvnorm
#' @import MASS
#'
#' @examples
#' #Load the India Human Development Survey-II dataset
#' data(ihd_mcar) 
#' 
#' # Gerbil without types specified
#' imps.gerbil <- gerbil(ihd_mcar, m = 1, mcmciter = 10)
#' \donttest{
#' # Gerbil with types specified (method #1)
#' types.gerbil <- c(
#'        sex = "binary", age = "continuous", 
#'        marital_status = "binary", job_field = "categorical", 
#'        farm_labour_days = "semicont", own_livestock = "binary", 
#'        education_level = "ordinal", income = "continuous")
#' imps.gerbil <- gerbil(ihd_mcar, m = 1, type = types.gerbil)
#' 
#' # Gerbil with types specified (method #2)
#' imps.gerbil <- gerbil(ihd_mcar, m = 1, ords = "education_level", semi = "farm_labour_days", 
#'        bincat = c("sex", "marital_status", "job_field", "own_livestock"))
#' 
#' # Gerbil with types specified (method #3)
#' types.gerbil <- c("binary", "continuous", "binary", "categorical", "semicont", 
#'        "binary", "ordinal", "continuous")
#' imps.gerbil <- gerbil(ihd_mcar, m = 1, type = types.gerbil)
#' 
#' # Variables of class factor are treated as binary/categorical by default
#' ihd.fac <- ihd_mcar
#' ihd.fac$sex <- factor(ihd_mcar$sex)
#' ihd.fac$marital_status <- factor(ihd_mcar$marital_status)
#' ihd.fac$job_field <- factor(ihd_mcar$job_field)
#' ihd.fac$own_livestock <- factor(ihd_mcar$own_livestock)
#' ihd.fac$education_level <- ordered(ihd_mcar$education_level)
#' imps.gerbil <- gerbil(ihd.fac, m = 1)
#' 
#' # Univariate plotting of one variable
#' plot(imps.gerbil, type = 1, y = "job_field")
#' 
#' # gerbil with predMat specified (method #1)
#' predMat <- matrix(c(1, 0, 0, 1), 2, 2)
#' dimnames(predMat) <- list(c("education_level", "income"), c("sex", "job_field"))
#' imps.gerbil <- gerbil(ihd_mcar, m = 1, type = types.gerbil, predMat = predMat)
#' 
#' # gerbil with predMat specified (method #2)
#' predMat <- rbind(
#'        c(0, 0, 0, 0, 0, 0, 0, 0), 
#'        c(1, 0, 0, 0, 0, 0, 0, 0), 
#'        c(1, 1, 0, 0, 0, 0, 0, 0), 
#'        c(1, 1, 1, 0, 0, 0, 0, 0), 
#'        c(1, 1, 1, 1, 0, 0, 0, 0), 
#'        c(1, 1, 1, 1, 1, 0, 0, 0), 
#'        c(1, 1, 1, 0, 1, 1, 0, 0), 
#'        c(0, 1, 1, 1, 1, 1, 1, 0) 
#'        )
#' imps.gerbil <- gerbil(ihd_mcar, type = types.gerbil, predMat = predMat)
#' 
#' # Multiple imputation with more iterations
#' imps.gerbil.5 <- gerbil(ihd_mcar, m = 5, mcmciter = 100, ords = "education_level", 
#'        semi = "farm_labour_days", bincat = "job_field", n.cores = 1)
#' 
#' plot(imps.gerbil.5, type = 1, y = "job_field", imp = 1:5) 
#' 
#' # Extract the first imputed dataset
#' imputed.gerb <- imputed(imps.gerbil.5, imp = 1)
#' 
#' # Write all imputed datasets to an Excel file
#' write.gerbil(imps.gerbil.5, file = file.path(tempdir(), "gerbil_example.xlsx"), imp = 1:5)
#' 
#' if(requireNamespace('mice')){
#' # Impute using mice for comparison
#' 
#' types.mice <- c("logreg", "pmm", "logreg", "polyreg", "pmm", "logreg", "pmm", "pmm")
#' imps.mice <- mice(ihd.fac, m = 1, method = types.mice, maxit = 100)
#' 
#' imps.mice1 <- mice(ihd.fac, m = 1, method = "pmm", maxit = 100)
#' 
#' imps.gerbil <- gerbil(ihd_mcar, m = 1, mcmciter = 100, ords = "education_level", 
#'     semi = "farm_labour_days", bincat = "job_field")
#' 
#' # Compare the performance of mice and gerbil
#' 
#' # Replace some gerbil datasets with mice datasets
#' imps.gerbil.m <- imps.gerbil.5
#' imps.gerbil.m$imputed[[2]] <- complete(imps.mice, action = 1)
#' imps.gerbil.m$imputed[[3]] <- complete(imps.mice1, action = 1)
#' 
#' # Perform comparative correaltion analysis
#' cor_gerbil(imps.gerbil.m, imp = 1, log = "income")
#' cor_gerbil(imps.gerbil.m, imp = 2, log = "income")
#' cor_gerbil(imps.gerbil.m, imp = 3, log = "income")
#' 
#' # Perform comparative univariate goodness-of-fit testing
#' gof_gerbil(imps.gerbil.m, type = 1, imp = 1)
#' gof_gerbil(imps.gerbil.m, type = 1, imp = 2)
#' gof_gerbil(imps.gerbil.m, type = 1, imp = 3)
#' 
#' # Perform comparative bivariate goodness-of-fit testing
#' gof_gerbil(imps.gerbil.m, type = 2, imp = 1)
#' gof_gerbil(imps.gerbil.m, type = 2, imp = 2)
#' gof_gerbil(imps.gerbil.m, type = 2, imp = 3)
#' 
#' # Produce univariate plots for comparisons 
#' plot(imps.gerbil.m, type = 1, file = file.path(tempdir(), "gerbil_vs_mice_univariate.pdf"), 
#'      imp = c(1, 2, 3), log = "income", lty = c(1, 2, 4, 5), col = c("blue4", "brown2", 
#'      "green3", "orange2"), legend = c("Observed", "gerbil", "mice: logistic", "mice: pmm"))
#' 
#' ### Produce bivariate plots for comparisons 
#' plot(imps.gerbil.m, type = 2, file = file.path(tempdir(), "gerbil_vs_mice_bivariate.pdf"), 
#'     imp = c(1, 2, 3), log = "income", lty = c(1, 2, 4, 5), col = c("blue4", "brown2", 
#'     "green3", "orange2"), pch = c(1, 3, 4, 5), legend = c("Observed", "gerbil", 
#'     "mice: logistic", "mice: pmm"))
#' 
#' }
#' }
gerbil <- function (dat, m = 1, mcmciter = 25, predMat = NULL, type = NULL,
    visitSeq = NULL, ords = NULL, semi = NULL, bincat = NULL,
    cont.meth = "EMP", num.cat = 12, r = 5, verbose = TRUE,
    n.cores = NULL, cl.type = NULL, mass = rep(0, length(semi)), ineligible = NULL, 
    trace = TRUE, seed = NULL)
{
    if (length(colnames(dat)) == 0) {
        colnames(dat) <- paste0("X", 1:NCOL(dat))
    }
    if (is.numeric(ords)) {
        ords <- colnames(dat)[ords]
    }
    cens <- NULL
    if (is.numeric(cens)) {
        cens <- colnames(dat)[cens]
    }
    if (is.numeric(semi)) {
        semi <- colnames(dat)[semi]
    }
    if (is.numeric(bincat)) {
        bincat <- colnames(dat)[bincat]
    }
    orig.nams <- colnames(dat)
    n.cores <- gerbilCluster(n.cores)
    if (!is.character(cens)) {
        cens <- colnames(dat)[cens]
    }
    is.df <- FALSE
    if (class(dat)[1] == "data.frame") {
        is.df <- TRUE
    }
    if (is.df) {
        classes <- as.matrix(data.frame(lapply(dat, base::class)))[1, ]
        facs <- which(classes == "factor")
        chars <- which(classes == "character")
        ints <- which(classes == "integer")
        ordereds <- which(classes == "ordered")
        facs <- setdiff(facs, which(is.element(colnames(dat), cens)))
        chars <- setdiff(chars, which(is.element(colnames(dat), cens)))
        ints <- setdiff(ints, which(is.element(colnames(dat), cens)))
        ordereds <- setdiff(ordereds, which(is.element(colnames(dat), cens)))
    }

    calc.var <- TRUE
    calc.latent <- FALSE

    if (length(ineligible) > 0) {
        ineligible.tmp <- matrix(FALSE, NROW(dat), NCOL(dat))
        colnames(ineligible.tmp) <- colnames(dat)
        if (length(ineligible) > 1) {
          if (length(colnames(ineligible)) == 0) {
            if (NCOL(ineligible) == NCOL(dat)) {
              colnames(ineligible) <- colnames(dat)
            } else {
              stop("Input 'ineligible' requires column names.")
            }            
          }
          if (sum(!is.element(colnames(ineligible), colnames(dat)))) {
            warning("Input 'ineligible' contains column names that don't appear in 'dat'.")
          }
          ineligible <- ineligible[, intersect(colnames(ineligible), colnames(dat))]
          ineligible <- ineligible == 1
          ineligible.tmp[, colnames(ineligible)] <- ineligible
        } else {
          ineligible.tmp[dat == ineligible | dat == as.character(ineligible)] <- TRUE
          ineligible.tmp[is.na(ineligible.tmp)] <- FALSE
        }
        if (sum(!is.na(dat[ineligible.tmp])) > 0) {
          warning("Some values marked as missing but ineligible for imputation (via the input 'ineligible') contain non-missing entries. \n   These will be set as NA.")
        }
        dat[ineligible.tmp] <- NA
    } else {
      ineligible.tmp <- NULL
    }

    if (length(seed) != 0) {
        set.seed(seed)
    }

    cens.dat <- matrix(0, NROW(dat), NCOL(dat))
    dimnames(cens.dat) <- dimnames(dat)
    for (i in 1:NCOL(dat)) {
        nchar.tmp <- nchar(as.character(dat[, i]))
        cens.dat[, i] <- substr(as.character(dat[, i]) ,1 ,1) == ">" | substr(as.character(dat[, i]), nchar.tmp , nchar.tmp) == "+"
    }
    cens.dat[is.na(cens.dat)] <- FALSE
    n.cens <- colSums(cens.dat, na.rm = TRUE)
    n.miss <- colSums(is.na(dat))
    n.obs <- NROW(dat) - n.miss - n.cens
    miss.rate <- n.miss/NROW(dat)
    cens.rate <- n.cens/NROW(dat)
    obs.rate <- n.obs/NROW(dat)
    if (sum(n.miss + n.cens) == 0) {
        stop("Input 'dat' contains no missing or censored values. (Missing values should be coded as NA.)")
    }
    miss.dat <- is.na(dat)

    bin.meth <- "binary"
    cat.meth <- "categorical"
    ord.meth <- "ordinal"
    cens.meth <- "censored"
    semi.meth <- "semicont"
    var.types <- c(bin.meth, cat.meth, ord.meth, cens.meth, semi.meth)
    var.types <- c(var.types, cont.meth)
    if (length(type) == NCOL(dat)) {
        if (length(names(type)) == 0) {
            names(type) <- colnames(dat)
            #warning("Input 'type' is unnamed.  Names will be inferred.")
        }
    }
    else if (length(type) == 0) {
        type <- rep(NA, NCOL(dat))
        names(type) <- colnames(dat)
    }
    else {
        if (length(names(type)) == 0) {
            stop("Input 'type' is unnamed.  Names are required given its length.")
        }
    }
    if (sum(!is.element(names(type), colnames(dat))) > 0) {
        warning("Input 'type' contains elements that do not match to the column names of input 'dat'.")
    }
    type1 <- rep(NA, NCOL(dat))
    names(type1) <- colnames(dat)
    type1[] <- type[match(names(type1), names(type))]
    type <- type1
    for (i in 1:NCOL(dat)) {
        if (is.na(type[i])) {
            num <- length(table(dat[, i]))
            if (is.element(colnames(dat)[i], bincat)) {
                if (num == 2) {
                  type[i] <- bin.meth
                }
                else {
                  type[i] <- cat.meth
                }
            }
            else if (is.element(colnames(dat)[i], ords)) {
                type[i] <- ord.meth
            }
            else if (is.element(colnames(dat)[i], cens)) {
                type[i] <- cens.meth
            }
            else if (is.element(colnames(dat)[i], semi)) {
                type[i] <- semi.meth
            }
            else if (n.cens[i] > 0) {
                type[i] <- cens.meth
            }
            else if (class(dat[, i])[1] == "ordered") {
                type[i] <- ord.meth
            }
            else if (num > 2 & (class(dat[, i])[1] == "factor" | class(dat[, i])[1] == "character")) {
                if (num == 2) {
                  type[i] <- bin.meth
                }
                else {
                  type[i] <- cat.meth
                }
                dat[, i] <- factor(dat[, i])
            }
            else {
                dat[is.nan(dat[, i]), i] <- NA
                if (num == 2) {
                  type[i] <- bin.meth
                }
                else if (num <= num.cat) {
                  type[i] <- cat.meth
                  dat[, i] <- factor(dat[, i])
                }
                else {
                  type[i] <- cont.meth
                  dat[, i] <- as.numeric(dat[, i])
                }
            }
        }
        else {
            if (type[i] == bin.meth) {
                num <- length(table(dat[, i]))
                if (num != 2) {
                  stop(paste0("Variable ", colnames(dat)[i],
                    " is marked as binary but does not take on exactly two values."))
                }
            } 
            else if (type[i] == cat.meth) {
                num <- length(table(dat[, i]))
                if (num == 2) {
                  warning(paste0("Variable ", colnames(dat)[i],
                    " is marked as categorical but takes on exactly two values. Resetting to binary."))
                  type[i] <- bin.meth
                }
            }
            else if (type[i] == "continuous") {
                type[i] <- cont.meth
            }
            else if (type[i] == cens.meth) {
                if (n.cens[i] == 0) {
                  stop(paste0("Variable ", colnames(dat)[i],
                    " is marked as censored but does not have any censored values."))
                }
            }
            else if (!is.element(type[i], var.types)) {
                stop(paste0("The variable type for ", names(type)[i],
                  " (", type[i], ") is not permissible."))
            }
        }
    }
    semi.vars <- which(type == "semicont")
    if (length(semi.vars) > 0) {
        if (length(semi.vars) == length(mass)) {
            if (length(names(mass)) == 0) {
                names(mass) <- names(type)[semi.vars]
            }
        }
        else if (length(mass) == 0) {
            mass <- rep(0, length(semi.vars))
            names(mass) <- names(type)[semi.vars]
        }
        else if (length(mass) == 1) {
            mass <- rep(mass, length(semi.vars))
            names(mass) <- names(type)[semi.vars]
        }
        else {
            if (sum(!is.element(names(mass), colnames(dat)[semi.vars])) >
                0) {
                warning("Input 'mass' contains elements that do not match to the assigned semi-continuous variables.")
            }
            mass1 <- rep(NA, length(semi.vars))
            names(mass1) <- colnames(dat)[semi.vars]
            mass1[] <- mass[match(names(mass1), names(mass))]
            mass1[is.na(mass1)] <- 0
            mass <- mass1
        }
        for (i in 1:length(semi.vars)) {
            num <- sum(dat[, semi.vars[i]] == mass[i], na.rm = TRUE)
            if (num == 0) {
                stop(paste0("Variable ", colnames(dat)[semi.vars[i]],
                  " is marked as semi-continuous but has no cases with value ",
                  mass[i], "."))
            }
        }
    }
    cens <- names(type)[type == cens.meth]
    semi <- names(type)[type == semi.meth]
    ords <- names(type)[type == ord.meth]
    dat1 <- dat
    continu <- which(type != "binary" & type != "categorical" & type != "ordinal" & type != "censored" & type != "semicont")
    type.info <- type
    type.info[continu] <- paste0("continuous (", type.info[continu], ")")

    if (length(cens) == 0) {
        info <- data.frame(Variable.Type = type.info, Num.Observed = n.obs,
            Num.Miss = n.miss, Miss.Rate = miss.rate)
    }
    else {
        info <- data.frame(Variable.Type = type, Num.Fully.Observed = n.obs,
            Num.Miss = n.miss, Num.Censored = n.cens, Miss.Rate = miss.rate,
            Cens.Rate = cens.rate)
    }
    rownames(info) <- colnames(dat1)
    if (verbose) {
        info.tmp <- info[, "Miss.Rate"]
        info.tmp <- formatC(100 * info.tmp, format = "f", digits = 2)
        info.tmp <- paste0(info.tmp, "%")
        info1 <- info
        info1[, "Miss.Rate"] <- info.tmp
        if (length(cens) > 0) {
            info.tmp <- info[, "Cens.Rate"]
            info.tmp <- formatC(100 * info.tmp, format = "f",
                digits = 2)
            info.tmp <- paste0(info.tmp, "%")
            info1[, "Cens.Rate"] <- info.tmp
        }
        message("Variable Summary: ", "\n", sep = "", appendLF = FALSE)
        print(info1)
        message("\n", appendLF = FALSE)
    }
    obs.vars <- rownames(info)[(1 - obs.rate) == 0]
    imp.vars <- rownames(info)[(1 - obs.rate) != 0]
    all.vars <- rownames(info)
    predictorMatrix1 <- (1 - diag(1, NCOL(dat)))
    predictorMatrix1[upper.tri(predictorMatrix1)] <- 0
    dimnames(predictorMatrix1) <- list(colnames(dat), colnames(dat))
    if (length(predMat) == 0) {
        predictorMatrix <- predictorMatrix1
    }
    else {
        predictorMatrix <- predMat
    }
    if (length(rownames(predictorMatrix)) == 0) {
        if (NROW(predictorMatrix) == length(all.vars)) {
            rownames(predictorMatrix) <- all.vars
        }
        else if (NROW(predictorMatrix) == length(imp.vars)) {
            rownames(predictorMatrix) <- imp.vars
        }
        else {
            stop("Rows of input 'predMat' cannot be matched to specific variables listed in input 'dat'.")
        }
    }
    if (length(colnames(predictorMatrix)) == 0) {
        if (NCOL(predictorMatrix) == length(all.vars)) {
            colnames(predictorMatrix) <- all.vars
        }
        else {
            stop("Columns of input 'predMat' cannot be matched to specific variables listed in input 'dat'.")
        }
    }
    over.row <- setdiff(rownames(predictorMatrix), all.vars)
    if (length(over.row) > 0) {
        warning("Rows in 'predMat' correspond to variables that do not exist in 'dat'. These rows will be removed.")
    }
    under.row <- setdiff(imp.vars, rownames(predictorMatrix))
    if (length(under.row) > 0) {
        warning("Variables in input 'dat' have missing values but do not correspond to a row of 'predMat'. Dependence structure will be assigned by default.")
    }
    predictorMatrix <- predictorMatrix[intersect(all.vars, rownames(predictorMatrix)),
        , drop = FALSE]
    over.col <- setdiff(colnames(predictorMatrix), all.vars)
    if (length(over.col) > 0) {
        warning("Columns in 'predMat' correspond to variables that do not exist in 'dat'. These columns will be removed.")
    }
    under.col <- setdiff(all.vars, colnames(predictorMatrix))
    if (length(under.col) > 0) {
        warning("Variables in input 'dat' do not correspond to a column of 'predMat'. Dependence structure will be assigned by default.")
    }
    predictorMatrix <- predictorMatrix[, intersect(all.vars,
        colnames(predictorMatrix)), drop = FALSE]
    predictorMatrix1[rownames(predictorMatrix), colnames(predictorMatrix)] <- predictorMatrix
    predictorMatrix <- predictorMatrix1
    if (is.numeric(visitSeq)) {
        visitSeq <- colnames(dat)[visitSeq]
    } else if (is.logical(visitSeq)) {
        if (visitSeq) {
            visitSeq <- colnames(dat1)[order(1 - obs.rate, decreasing = FALSE)]
            #visitSeq <- colnames(dat1)
        }
        else {
            visitSeq <- colnames(dat1)
        }
    } else if (length(visitSeq) == 0 & length(rownames(predMat)) == 0) {
        #visitSeq <- colnames(dat1)[order(1 - obs.rate, decreasing = FALSE)]
        visitSeq <- colnames(dat1)
    } else if (length(visitSeq) == 0 & length(rownames(predMat)) > 0) {
        visitSeq <- rownames(predictorMatrix)
    } else {
        if (is.numeric(visitSeq)) {
            if (max(visitSeq) > NCOL(dat1) | min(visitSeq) <
                1) {
                visitSeq <- visitSeq[visitSeq >= 1 & visitSeq <=
                  NCOL(dat1)]
                warning("Input 'visitSeq' indicates variables that are not contained in input 'dat'.")
            }
            visitSeq <- colnames(dat1)[visitSeq]
        }
        else {
            visitSeq <- intersect(visitSeq, colnames(dat1))
            if (sum(!is.element(visitSeq, colnames(dat1)))) {
                warning("Input 'visitSeq' indicates variables that are not contained in input 'dat'.")
            }
        }
        if (sum(!is.element(rownames(info)[(1 - obs.rate) > 0],
            visitSeq)) > 0) {
            stop("There are variables with missing data that do not appear in input 'visitSeq'.")
        }
        visitSeq <- setdiff(visitSeq, obs.vars)
        visitSeq <- c(obs.vars, visitSeq)
    }
    if (sum(diag(predictorMatrix) != 0) > 0) {
        diag(predictorMatrix) <- 0
        warning("Setting diagonal elements of the (resorted) input 'predMat' to zero.")
    }
    if (sum(predictorMatrix[upper.tri(predictorMatrix)] != 0) >
        0) {
        predictorMatrix[upper.tri(predictorMatrix)] <- 0
        warning("Some elements in the upper triangle of input 'predMat' are non-zero.  These will be set to zero.")
    }
    predictorMatrix[upper.tri(predictorMatrix)] = t(predictorMatrix)[upper.tri(t(predictorMatrix))]
    if (sum(rownames(predictorMatrix) != visitSeq) > 0 & length(predMat) >
        1) {
        warning("Input 'visitSeq' implies different variable ordering than input 'predictorMatrix'.  Permissible variable dependencies will be inferred.")
    }
    type <- type[visitSeq]
    dat1 <- dat1[, visitSeq]
    predictorMatrix <- predictorMatrix[visitSeq, visitSeq]
    visitSeq <- intersect(visitSeq, imp.vars)
    predictorMatrix[upper.tri(predictorMatrix)] <- 0
    predictorMatrix.old <- predictorMatrix[visitSeq, ]
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
            cont.meth = cont.meth, mass = mass, categories = categories)
        obs.latent <- tmpdat[[2]]
        predictorMatrix <- tmpdat[[3]]
        cat.nams <- tmpdat[[4]]
        cat.cols <- tmpdat[[5]]
        is.cat <- tmpdat[[6]]
        type <- tmpdat[[7]]
        semicols.new <- tmpdat[[8]]
        nams.out <- tmpdat[[9]]
        dat1 <- tmpdat[[1]]
    } else {
        obs.latent <- matrix(as.numeric(miss.dat), NROW(dat1), NCOL(dat1))
        dimnames(obs.latent) <- dimnames(miss.dat)
        cat.nams <- cat.cols <- is.cat <- semicols.new <- NULL
        nams.out <- colnames(dat1)
        names(nams.out) <- colnames(dat1)
    }
    if (length(ineligible) > 0) {
      obs.latent[ineligible.tmp[, nams.out]] <- 4
    }
    bincols <- which(type == bin.meth)
    ordcols <- which(type == ord.meth)
    censcols <- which(type == cens.meth)
    transcols <- which(type == cont.meth)
    cens <- colnames(cens.dat)[colSums(cens.dat) > 0]
    if (length(cens) > 0) {
      obs.latent[, cens] <- 2 * cens.dat[, cens] + obs.latent[, cens]
    }
    if (length(censcols) > 0) {
        cens <- matrix(as.character(dat1[, censcols, drop = FALSE]), NROW(dat1), length(censcols))
        cens <- substr(cens, 1, 1) == ">" | substr(cens, nchar(cens), nchar(cens)) == "+"
        dat1[, censcols] <- gsub(">", "", as.matrix(dat1[, censcols]), fixed = TRUE)
        dat1[, censcols] <- gsub("+", "", as.matrix(dat1[, censcols]), fixed = TRUE)
        dat1[, censcols] <- as.numeric(dat1[, censcols])
        colnames(cens) <- colnames(dat1)[censcols]
        cens[is.na(cens)] <- FALSE
    }
    dat1 <- as.matrix(dat1)
    parms <- list()
    varTrans <- NULL
    censTrans <- NULL
    if (length(transcols) > 0) {
        varTrans <- rep(cont.meth, length(transcols))
        tmpdat0 <- gerbil.trans(dat1, transcols = transcols,
            transtype = varTrans, nams = colnames(dat1)[transcols],
            verbose = verbose)
        parms <- tmpdat0[[2]]
        tmpdat0 <- tmpdat0[[1]]
    }
    else {
        tmpdat0 <- dat1
    }
    if (length(censcols) > 0) {
        censTrans <- rep(cens.meth, length(censcols))
        tmpdat0 <- gerbil.trans.cens(tmpdat0, transcols = censcols,
            transtype = censTrans, nams = colnames(dat1)[censcols],
            cens = cens, verbose = verbose)
        parms[(length(transcols) + 1):(length(transcols) + length(censcols))] <- tmpdat0[[2]]
        tmpdat0 <- tmpdat0[[1]]
    }
    verbose1 <- verbose
    if (n.cores > 1 & m > 1) {
        verbose1 <- FALSE
        if (verbose) {
            message("Parallelizing with n.cores = ", n.cores,
                "...\n", sep = "", appendLF = FALSE)
        }
        requireNamespace("parallel", quietly = TRUE)
        if (length(cl.type) == 0) {
          cl <- parallel::makeCluster(n.cores)
        } else {
          cl <- parallel::makeCluster(n.cores, type = cl.type)
        }
        first <- proc.time()
        list.out <- parallel::parLapply(cl = cl, X = 1:m, gerbil.inner,
            data = tmpdat0, predictorMatrix = predictorMatrix,
            obs.latent = obs.latent, cat.cols = cat.cols, is.cat = is.cat,
            mcmciter = mcmciter, bincols = bincols, ordcols = ordcols,
            censcols = censcols, cens = cens, r = r, verbose = FALSE,
            semicols.new = semicols.new, calc.var = calc.var, calc.latent = calc.latent, 
            trace = trace)
        last <- proc.time()
        parallel::stopCluster(cl)
        if (verbose) {
            message("Created ", m, " imputed datasets.  Total time = ",
                formatC((last[3] - first[3]), format = "f", digits = 2),
                ".\n", sep = "", appendLF = FALSE)
        }
    }
    tmp.final <- proc.time()
    out <- list()
    for (i in 1:m) {
        if (n.cores == 1 | m == 1) {
            tmpdat <- gerbil.inner(X = i, tmpdat0, predictorMatrix,
                obs.latent = obs.latent, cat.cols = cat.cols, is.cat = is.cat,
                mcmciter = mcmciter, bincols = bincols, ordcols = ordcols,
                censcols = censcols, cens = cens, r = r, verbose = verbose,
                semicols.new = semicols.new, calc.var = calc.var, calc.latent = calc.latent, 
                trace = trace)
        }
        else {
            tmpdat <- list.out[[i]]
        }
        if (trace) {
          if (i == 1) {
              means.all <- means.mis <- means.all.latent <- means.mis.latent <- array(NA, c(dim(tmpdat$means.all), m))
              dimnames(means.all) <- dimnames(means.mis) <- dimnames(means.all.latent) <- dimnames(means.mis.latent) <- list(dimnames(tmpdat$means.all)[[1]],
                  dimnames(tmpdat$means.all)[[2]], 1:m)
              if (calc.var) {
                vars.all <- vars.mis <- vars.all.latent <- vars.mis.latent <- array(NA, c(dim(tmpdat$vars.all), m))
                dimnames(vars.all) <- dimnames(vars.mis) <- dimnames(vars.all.latent) <- dimnames(vars.mis.latent) <- list(dimnames(tmpdat$vars.all)[[1]],
                    dimnames(tmpdat$vars.all)[[2]], 1:m)
              } else {
                vars.all <- vars.mis <- vars.all.latent <- vars.mis.latent <- NULL
              }
          }
          if (i == 1) {
            means.obs <- tmpdat$means.obs
            vars.obs <- tmpdat$vars.obs
          }
          means.all[, , i] <- tmpdat$means.all
          means.mis[, , i] <- tmpdat$means.mis
          if (length(tmpdat$means.mis.latent) > 0) {
            means.all.latent[, , i] <- tmpdat$means.all.latent
            means.mis.latent[, , i] <- tmpdat$means.mis.latent
          }
          if (calc.var) {
            vars.all[, , i] <- tmpdat$vars.all
            vars.mis[, , i] <- tmpdat$vars.mis
            if (length(tmpdat$vars.mis.latent) > 0) {
              vars.all.latent[, , i] <- tmpdat$vars.all.latent
              vars.mis.latent[, , i] <- tmpdat$vars.mis.latent
            }
          }
        }
        if (i == 1) {
            drops <- tmpdat$drops
            predMat <- tmpdat$predictorMatrix
        }
        tmpdat <- tmpdat$imputes
        if ((length(transcols) + length(censcols)) > 0) {
            out1 <- gerbil.untrans(dat1, tmpdat, parms = parms,
                transcols = c(transcols, censcols), transtype = c(varTrans,
                  censTrans), nams = colnames(dat1)[c(transcols,
                  censcols)], cens = cens, verbose = verbose1,
                X = i)
        }
        else {
            out1 <- tmpdat
        }
        if (length(catcols) > 0 | length(semicols) > 0) {
            out1 <- make.uncat(out1, cat.cols, is.cat, cat.nams,
                mass)
        }
        if (is.df) {
            out1 <- data.frame(out1)
        }
        if (length(categories) > 0) {
            for (k in 1:length(categories)) {
                cat.tmp <- categories[[k]]
                suppressWarnings(is.num <- sum(is.na(as.numeric(cat.tmp))) == 0)
                tmp <- out1[, names(categories)[k]]
                tmp1 <- rep(NA, length(tmp))
                for (j in 1:length(cat.tmp)) {
                  tmp1[tmp == j - 1] <- cat.tmp[j]
                }
                if (is.num) {
                  tmp1 <- as.numeric(as.character(tmp1))
                }
                else {
                  tmp1 <- factor(tmp1, levels = cat.tmp)
                }
                out1[, names(categories)[k]] <- tmp1
            }
        }

        if (length(ineligible) == 1) {
            out1[ineligible.tmp[, colnames(out1)]] <- ineligible
        } else if (length(ineligible) > 1) {
            out1[ineligible.tmp[, colnames(out1)]] <- NA
        }

        if (is.df) {
            if (length(facs) > 0) {
                for (j in facs) {
                  out1[, j] <- factor(out1[, j])
                }
            }
            if (length(ordereds) > 0) {
                for (j in ordereds) {
                  out1[, j] <- factor(out1[, j], ordered = T)
                }
            }
            if (length(chars) > 0) {
                for (j in chars) {
                  out1[, j] <- as.character(out1[, j])
                }
            }
            if (length(ints) > 0) {
                for (j in ints) {
                  out1[, j] <- as.integer(out1[, j])
                }
            }
        }
        out[[i]] <- out1[, orig.nams]

    }
    names(out) <- 1:m
    tmp.final <- proc.time() - tmp.final
    if (verbose & !verbose1) {
        message("Completed untransformations and post-processing, Time = ",
            formatC(tmp.final[3], format = "f", digits = 2),
            "\n", sep = "", appendLF = FALSE)
    }
    if (m > 2 & mcmciter >= 4 & trace) {
      R.hat <- matrix(NA, 4, dim(means.all)[1])
      dimnames(R.hat) <- list(c("Means", "Latent.Process.Means", "Variances", "Latent.Process.Variances"), dimnames(means.all)[[1]])
      use <- floor((mcmciter + 1)/2):mcmciter
      R.hat[1, ] <- apply(means.all[, use, , drop = FALSE], 1, get.R.hat)
      if (calc.latent) {
        R.hat[2, ] <- apply(means.all.latent[, use, , drop = FALSE], 1, get.R.hat)
      }
      if (calc.var) {
        R.hat[3, ] <- apply(vars.all[, use, , drop = FALSE], 1, get.R.hat)
        if (calc.latent) {
          R.hat[4, ] <- apply(vars.all.latent[, use, , drop = FALSE], 1, get.R.hat)
        }
      } else {
        R.hat <- R.hat[1:2, ]
      }
      R.hat <- t(R.hat)
      if (!calc.latent) {
        R.hat <- R.hat[, intersect(colnames(R.hat), c("Means", "Variances")), drop = FALSE]
      }
    }
    else {
      R.hat <- NULL
    }
    forms <- rep(NA, NROW(predMat))
    names(forms) <- rownames(predMat)
    for (i in 1:NROW(predMat)) {
        suff <- paste(colnames(predMat)[predMat[i, ] != 0], collapse = " + ")
        if (length(suff) == 0 | suff == "") {
            suff <- as.character(1)
        }
        forms[i] <- paste0(rownames(predMat)[i], " ~ ", suff)
    }
    if(!calc.var & trace) {
      if (calc.latent) {
        seq <- list(means.all = means.all, means.mis = means.mis, means.all.latent = means.all.latent, means.mis.latent = means.mis.latent, means.obs = means.obs)
      } else {
        seq <- list(means.all = means.all, means.mis = means.mis, means.obs = means.obs)
      }
    } else if (trace) {
      if (calc.latent) {
        seq <- list(means.all = means.all, means.mis = means.mis, means.all.latent = means.all.latent, means.mis.latent = means.mis.latent, means.obs = means.obs,
                    vars.all = vars.all, vars.mis = vars.mis, vars.all.latent = vars.all.latent, vars.mis.latent = vars.mis.latent, vars.obs = vars.obs)
      } else {
        seq <- list(means.all = means.all, means.mis = means.mis, means.obs = means.obs,
                    vars.all = vars.all, vars.mis = vars.mis, vars.obs = vars.obs)
      }
    } else {
      seq <- NULL
    }

    miss.dat.out <- miss.dat + 2 * cens.dat
    if (length(ineligible.tmp) > 0) {
      miss.dat.out <- miss.dat.out + 3 * ineligible.tmp
    }

    gerbil_object = list(imputed = out, missing = miss.dat.out, summary = info, chainSeq = seq,
        R.hat = R.hat, missing.latent = obs.latent, visitSeq.initial = visitSeq, visitSeq.final = rownames(predMat),
        predMat.initial = predictorMatrix1, predMat.final = predMat, drops = drops, forms = forms,
        mass.final = mass, ineligibles = ineligible.tmp, nams.out = nams.out)
    class(gerbil_object) = "gerbil"
    return(gerbil_object)
}


gerbil.trans <- function(dat, transcols = 1:NCOL(dat), transtype = rep("NONE", NCOL(dat)),
                         nams = (1:NCOL(dat))[transcols], verbose = TRUE) {

  tmp.time.all <- proc.time()
  dat <- as.data.frame(dat)
  out <- dat

  parms <- list()
  for (i in 1:length(transcols)) {
    if (transtype[i] == "EMP") {
      tmptime <- proc.time()
      tmp <- dat[!is.na(dat[, transcols[i]]), transcols[i]]
      tmp1 <- pemp(tmp)
      tab <- tmp1[[2]]
      tmp1 <- tmp1[[1]]
      tmp1 <- qnorm(tmp1)
      parms[[i]] <- list(tab,c(mean(tmp1), sd(tmp1)))
      out[!is.na(dat[, transcols[i]]), transcols[i]] <- tmp1
      tmptime <- proc.time() - tmptime
      #cat("Transformed ",nams[i]," w/ Emp. Dist., Time = ",formatC(tmptime[3], format = "f", digits = 2),", n = ",length(tmp),"\n",sep="")
    }
  }

  tmp.time.all <- proc.time() - tmp.time.all
  if (verbose) {
    message("Completed transformations, Time = ",formatC(tmp.time.all[3], format = "f", digits = 2), "\n", sep = "", appendLF = FALSE)
  }

  return(list(out, parms))

}


make.cat <- function(dat, predictorMatrix = NULL, catcols = NULL, semicols = NULL,
                     type = rep("EMP", 1:NCOL(data)), bin.meth = "trunc",
                     cont.meth = "EMP", mass = rep(0, length(semicols)), categories = NULL, 
                     reord = TRUE, impose.miss = TRUE)
{

  nams.out <- colnames(dat)

  allcols <- c(catcols, semicols)
  allcols <- allcols[order(allcols, decreasing = FALSE)]

  dat <- data.frame(dat)
  obs <- matrix(as.numeric(is.na(dat)), NROW(dat), NCOL(dat))
  dimnames(obs) <- dimnames(dat)
  obs <- as.data.frame(obs)
  var.nams <- colnames(dat)[allcols]

  cat.nams <- cols <- list()
  cat.nams[length(allcols)] <- NA

  old.opt <- options()
  on.exit(options(old.opt))
  options(na.action = 'na.pass')

  for (j in 1:length(allcols)) {
    col <- var.nams[j]
    col1 <- which(colnames(dat)==col)
    #col1 <- catcols[j]
    #col <- colnames(dat)[col1]

    predictorMatrix <- data.frame(predictorMatrix)
    colnames(predictorMatrix) <- rownames(predictorMatrix) <- colnames(dat)

    num.na <- sum(is.na(dat[, col]))

    if (is.element(allcols[j], catcols)) {
      #options(na.action = 'na.pass')
      tmp <- model.matrix(~x-1, data = data.frame(x = factor(dat[, col])))
      #options(na.action = 'na.omit')

      out.nams <- levels(factor(dat[, col]))
      if(length(categories) > 0) {
        colnames(tmp) <- categories[[col]]
      } else {
        colnames(tmp) <- out.nams
      }
      if (reord) {
        ords <- order(colSums(tmp, na.rm = TRUE), decreasing = FALSE)
        tmp <- tmp[,ords]
        cat.nams[[j]] <- out.nams[ords]
      } else {
        cat.nams[[j]] <- out.nams
      }
      tmp1 <- tmp
      obs.tmp <- matrix(obs[, col], NROW(obs), NCOL(tmp))
      dimnames(obs.tmp) <- dimnames(tmp)
      #obs.tmp <- as.data.frame(obs.tmp)
      if (num.na > 0 & impose.miss) {
        for (i in 2:NCOL(tmp)) {
          here.tmp <- rowSums(tmp[, 1:(i - 1), drop = FALSE]) == 1
          tmp1[here.tmp, i] <- NA
          obs.tmp[here.tmp, i] <- 3
        }
      }
      tmp1 <- tmp1[, -NCOL(tmp1), drop = FALSE]
      obs.tmp <- obs.tmp[, -NCOL(obs.tmp), drop = FALSE]
    } else {
      mass.tmp <- mass[which(semicols == allcols[j])]
      obs.tmp <- matrix(obs[, col], NROW(dat), 2)
      here.tmp <- dat[, col] == mass.tmp
      bin <- as.integer(!here.tmp)
      cont <- dat[, col]
      cont[here.tmp] <- NA
      obs.tmp[here.tmp, 2] <- 3
      tmp1 <- as.matrix(cbind(bin, cont))
      colnames(tmp1) <- colnames(obs.tmp) <- paste0(col, c(".B", ""))
      #colnames(tmp1) <- paste0("ZZZZXCVWFEDHZXCV",c(".B",""))
      #colnames(tmp1) <- c(".B","")
    }

    cols[[j]] <- col1:(col1 + NCOL(tmp1) - 1)

    dat <- as.list(dat)
    dat[[col]] <- tmp1
    dat <- data.frame(dat)

    obs <- as.list(obs)
    obs[[col]] <- obs.tmp
    obs <- data.frame(obs)

    tmp2 <- matrix(predictorMatrix[, col], nrow = NROW(predictorMatrix), ncol = NCOL(tmp1))
    colnames(tmp2) <- colnames(tmp1)

    nams <- rownames(predictorMatrix)
    predictorMatrix <- as.list(predictorMatrix)
    predictorMatrix[[col]] <- tmp2
    predictorMatrix <- data.frame(predictorMatrix)
    rownames(predictorMatrix) <- nams

    predictorMatrix <- t(predictorMatrix)
    predictorMatrix <- as.matrix(predictorMatrix)
    predictorMatrix <- data.frame(predictorMatrix)

    tmp2 <- matrix(predictorMatrix[, col], nrow = NROW(predictorMatrix), ncol = NCOL(tmp1))
    colnames(tmp2) <- colnames(tmp1)

    nams <- rownames(predictorMatrix)
    predictorMatrix <- as.list(predictorMatrix)
    predictorMatrix[[col]] <- tmp2
    predictorMatrix <- data.frame(predictorMatrix)
    rownames(predictorMatrix) <- nams

    predictorMatrix <- t(predictorMatrix)
    predictorMatrix <- as.matrix(predictorMatrix)

    if (is.element(allcols[j], semicols)) {
      colnames(dat)[cols[[j]]] <- colnames(obs)[cols[[j]]] <- colnames(predictorMatrix)[cols[[j]]] <- rownames(predictorMatrix)[cols[[j]]] <- colnames(tmp1)
    }

    lower <- type[(1:length(type)) < col1]
    upper <- type[(1:length(type)) > col1]
    lower1 <- nams.out[(1:length(nams.out)) < col1]
    upper1 <- nams.out[(1:length(nams.out)) > col1]
    if (is.element(allcols[j], catcols)) {
      middle <- rep(bin.meth, NCOL(tmp1))
      middle1 <- rep(col, NCOL(tmp1))
    } else {
      middle <- c(bin.meth, cont.meth)
      middle1 <- c(col, col)
    }
    type <- c(lower, middle, upper)
    nams.out <- c(lower1, middle1, upper1)
  }

  names(cat.nams) <- names(cols) <- var.nams

  names(nams.out) <- colnames(dat)

  is.cat <- which(is.element(allcols,catcols))
  #is.semi <- which(is.element(allcols,semicols))
  #cat.cols <- cols[is.cat]
  #cat.nams <- cat.nams[is.cat]
  #semi.cols <- cols[is.semi]

  if(length(is.cat) < length(cols)) {
    semicols.new <- as.matrix(data.frame(cols[-is.cat]))
  } else {
    semicols.new <- NULL
  }

  return(list(dat, obs, predictorMatrix, cat.nams, cols, is.cat, type, semicols.new = semicols.new, nams.out = nams.out))

}


make.uncat  <-  function(dat, cols, is.cat, nams, mass)
{

  #dat <- data.frame(dat)
  is.semi <- setdiff(1:length(cols), is.cat)

  rm <- NULL
  for (i in 1:length(cols)) {
    cols.tmp <- cols[[i]]
    rm <- c(rm, cols.tmp[-1])

    if (is.element(i, is.cat)) {
      nams.tmp <- nams[[i]]
      tmp <- rep(NA,NROW(dat))
      here <- dat[, cols.tmp[1]] == 1
      tmp[here] <- nams.tmp[1]
      for (j in 2:length(cols.tmp)) {
        here <- is.na(tmp) & dat[, cols.tmp[j]] == 1
        tmp[here] <- nams.tmp[j]
      }
      tmp[is.na(tmp)] <- nams.tmp[length(nams.tmp)]
      dat[, cols.tmp[1]] <- as.numeric(as.character(tmp))
      colnames(dat)[cols.tmp[1]] <- names(cols)[i]
    } else {
      mass.tmp <- mass[which(is.semi == i)]
      tmp <- dat[, cols.tmp[2]]
      tmp[dat[, cols.tmp[1]] == 0] <- mass.tmp
      dat[, cols.tmp[1]] <- tmp
      nam <- colnames(dat)[cols.tmp[2]]
      colnames(dat)[cols.tmp] <- paste0(nam, c("", ".Y"))
    }
  }

  dat <- dat[, -rm]

  return(dat)

}


gerbil.trans.cens <- function (dat, transcols = 1:NCOL(dat), transtype = rep("NONE", NCOL(dat)),
                               nams = (1:NCOL(dat))[transcols], cens = NULL, verbose = TRUE)
{

  tmp.time.all <- proc.time()
  dat <- as.data.frame(dat)
  out <- dat

  parms <- list()

  for (i in 1:length(transcols)) {
    if (transtype[i] == "censored") {
      tmptime <- proc.time()
      tmp <- dat[!is.na(dat[,transcols[i]]),transcols[i]]
      cens.tmp <- cens[!is.na(dat[, transcols[i]]), colnames(dat)[transcols[i]]]
      tmp1 <- pemp.cens(tmp, cens.tmp)
      tab <- tmp1[[2]]
      tmp1 <- tmp1[[1]]
      tmp1 <- qnorm(tmp1)
      parms[[i]] <- list(tab,c(mean(tmp1),sd(tmp1)))
      out[!is.na(dat[,transcols[i]]),transcols[i]] <- tmp1
      tmptime <- proc.time()-tmptime
      #cat("Transformed ", nams[i], " w/ Emp. Cens. Dist., Time = ", formatC(tmptime[3], format = "f", digits = 2), ", n = ", length(tmp), "\n", sep = "")
    }
  }

  tmp.time.all <- proc.time() - tmp.time.all
  if (verbose) {
    message("Completed transformations, Time = ", formatC(tmp.time.all[3], format = "f", digits = 2), "\n", sep = "", appendLF = FALSE)
  }

  return(list(out, parms, cens))

}


gerbil.inner <- function (X = 1, data, predictorMatrix, mcmciter = 100, obs = NULL, obs.latent = NULL, cat.cols = NULL, is.cat = NULL,
                          mis = NULL, bincols = NULL, ordcols = NULL, censcols = NULL, r = 5,
                          cens = NULL, verbose = TRUE, semicols.new = NULL, calc.var = FALSE, calc.latent = FALSE, trace = TRUE)
{

  ### Define gerbil.init within gerbil.inner
  gerbil.init <- function (dat, obs, predictorMatrix ,mis = NULL, p = length(mis), r = 5,
                           bincols = NULL, ordcols = NULL, censcols = NULL, cens = NULL, imp = 1, verbose = TRUE) {

    ### Define get.bad within gerbil.init
    get.bad <- function(dat, obs, cens = NULL, predictorMatrix, mis, p, r = 5, verbose = TRUE, imp = 1) {

      ### Define colspace within get.bad
      colspace <- function(A, r = nrow(A), k = ncol(A), verbose = TRUE) {

        ### Define RREF within colspace
        RREF <- function(X, ...) {

          GaussianElimination <- function(A, B, tol=sqrt(.Machine$double.eps),
                                        verbose=FALSE, fractions=FALSE){
            # A: coefficient matrix
            # B: right-hand side vector or matrix
            # tol: tolerance for checking for 0 pivot
            # verbose: if TRUE, print intermediate steps
            # fractions: try to express nonintegers as rational numbers
            # If B is absent returns the reduced row-echelon form of A.
            # If B is present, reduces A to RREF carrying B along.

            # CRAN Check Problem - 1

            # PNL - This is a bad idea - package dependencies should not be resolved in run time.

            # I'm removing this and importing mass by default, but we should disambiguate calls to the MASS package.

            # if (fractions) {
            #   mass <- require(MASS)
            #   if (!mass) stop("fractions=TRUE needs MASS package")
            # }

            if ((!is.matrix(A)) || (!is.numeric(A)))
              stop("argument must be a numeric matrix")
            n <- nrow(A)
            m <- ncol(A)
            if (!missing(B)){
              B <- as.matrix(B)
              if (!(nrow(B) == nrow(A)) || !is.numeric(B))
                stop("argument must be numeric and must match the number of row of A")
              A <- cbind(A, B)
            }
            i <- j <- 1
            while (i <= n && j <= m){
              while (j <= m){
                currentColumn <- A[,j]
                currentColumn[1:n < i] <- 0
                # find maximum pivot in current column at or below current row
                which <- which.max(abs(currentColumn))
                pivot <- currentColumn[which]
                if (abs(pivot) <= tol) { # check for 0 pivot
                  j <- j + 1
                  next
                }
                if (which > i) A[c(i, which),] <- A[c(which, i),] # exchange rows
                A[i,] <- A[i,]/pivot # pivot
                row <- A[i,]
                A <- A - outer(A[,j], row) # sweep
                A[i,] <- row # restore current row
                if (verbose) if (fractions) print(fractions(A))
                else print(round(A, round(abs(log(tol,10)))))
                j <- j + 1
                break
              }
              i <- i + 1
            }
            # 0 rows to bottom
            zeros <- which(apply(A[, 1:m], 1, function(x) max(abs(x)) <= tol))
            if (length(zeros) > 0){
              zeroRows <- A[zeros, ]
              A <- A[-zeros, ]
              A <- rbind(A, zeroRows)
            }
            rownames(A) <- NULL
            if (fractions) fractions (A) else round(A, round(abs(log(tol, 10))))
          }

          GaussianElimination(X, ...)
          # returns the reduced row-echelon form of X
        }

        ## Returns the columns that create singularity in the matrix A
        if (k > r) {
          if (verbose){print("WARNING: More columns than rows in colspace.")}
        }
        rA <- RREF(A)
        pivot <- 1
        out <- NULL
        for (i in 2:k) {
          if (rA[(pivot + 1), i] == 1) {
            pivot <- pivot + 1
          } else {
            out <- c(out, i)
          }
        }
        return(out)
      }

      drops <- list()

      if (length(cens) > 0) {
        obs[, colnames(cens)] <- obs[, colnames(cens)] & !cens
      }

      myobs <- crossprod(as.matrix(obs))

      for (i in 1:p) {

        nam <- colnames(dat)[mis[i]]
        here <- which(colnames(predictorMatrix) == nam)
        here1 <- which(rownames(predictorMatrix) == nam)

        drops[[i]] <- NA

        use <- obs[, nam]
        new <- !use
        use <- which(use)
        new <- which(new)

        if (here == 1) {
          keep <- NULL
        } else {
          rem <- predictorMatrix[here1, 1:(here-1)]
          rem <- which(rem == 0)
          if (length(rem) > 0) {
            keep <- setdiff(1:(here-1), rem)
          } else {
            keep <- 1:(here-1)
          }
        }

        if (length(keep) > 0) {

          X <- as.matrix(cbind(1, dat[use, keep, drop = FALSE]))
          out <- colspace(crossprod(X))
          if (length(out) > 0) {
            out <- out - 1
          }

          out <- c(keep)[out]
          out1 <- keep[myobs[nam, keep] < r]
          bad <- union(out, out1)

          if (length(bad) > 0) {
            predictorMatrix[here1, bad] <- 0
            if (imp == 1 & verbose) {
              message("Dropping from conditional model for ", rownames(predictorMatrix)[here1], ": ", paste(colnames(predictorMatrix)[bad], collapse = ", "), "\n", appendLF = FALSE)
            }
            drops[[i]] <- colnames(predictorMatrix)[bad]
          }

        }

        dat[new, nam] <- sample(dat[use, nam], length(new) ,replace = TRUE)

      }

      names(drops) <- colnames(dat)[mis]

      return(list(predictorMatrix = predictorMatrix, drops = drops))

    }

    colnames(obs) <- colnames(dat)

    if (length(mis) == 0) {
      mis <- which(colSums(obs == 0) > 0)
      p <- length(mis)
    }

    drops <- get.bad(dat = dat, obs = !is.na(dat), cens = cens, predictorMatrix = predictorMatrix, mis = mis, p = p, r = r, imp = imp, verbose = verbose)
    predictorMatrix <- drops$predictorMatrix
    drops <- drops$drops

    cutpts <- lowupp <- list()

    for (i in 1:p) {
      nam <- colnames(dat)[mis[i]]
      here <- which(colnames(predictorMatrix) == nam)
      here1 <- which(rownames(predictorMatrix) == nam)

      use <- !is.na(dat[, nam])
      if (is.element(nam, colnames(cens))) {
        use <- use & !is.na(cens[, nam])
      }
      new <- !use
      use <- which(use)
      new <- which(new)

      if (here == 1) {
        keep <- NULL
      } else {
        keep <- predictorMatrix[here1, 1:(here-1)]
        keep <- which(keep == 1)
      }

      if (length(keep) == 0) {
        if (is.element(mis[i], bincols)) {
          Y <- as.numeric(as.character(dat[use, nam]) == levels(factor(dat[, nam]))[2])
          alldat <- cbind(rep(1, NROW(dat)))
          coeff <- qnorm(mean(Y))
          means <- c(tcrossprod(alldat, t(coeff)))

          lims <- as.character(dat[, nam])
          upper <- rep(Inf, length(lims))
          lower <- rep(-Inf, length(lims))
          upper[as.character(lims) == levels(factor(lims))[1] & !is.na(lims)] <- 0
          lower[as.character(lims) == levels(factor(lims))[2] & !is.na(lims)] <- 0
          lowupp[[i]] <- cbind(lower = lower, upper = upper)

          dat[, nam] <- truncnorm::rtruncnorm(NROW(dat), a = lower, b = upper, mean = means, sd = 1)
        } else if (is.element(mis[i], censcols)) {
          use.c <- obs[, nam] == 0
          upper <- rep(Inf, sum(use.c))
          lower <- cens[use.c, nam]
          lower[is.na(lower)] <- -Inf
          lowupp[[i]] <- cbind(lower = lower, upper = upper)

          dat[use.c, nam] <- truncnorm::rtruncnorm(sum(use.c), a = lower, b = upper, mean = 0, sd = 1)
        } else if (is.element(mis[i], ordcols)) {
          Y <- ordered(dat[, nam])
          means <- rep(0, NROW(dat))
          zetas <- table(Y[use])
          levs <- names(zetas)
          zetas <- cumsum(zetas)
          zetas <- zetas[-length(zetas)]/zetas[length(zetas)]
          zetas <- c(-Inf, qnorm(zetas), Inf)
          cutpts[[which(ordcols == mis[i])[1]]] <- zetas
          zetas.l <- zetas[-length(zetas)]
          zetas.u <- zetas[-1]

          lower <- zetas.l[match(Y, levs)]
          upper <- zetas.u[match(Y, levs)]
          lower[is.na(lower)] <- -Inf
          upper[is.na(upper)] <- Inf
          lowupp[[i]] <- cbind(lower = lower, upper = upper)

          dat[, nam] <- truncnorm::rtruncnorm(NROW(dat), a = lower, b = upper, mean = means, sd = 1)
        } else {
          mod <- lm(dat[use, nam]~1)
          S2 <- (summary(mod)$sigma)^2
          beta.hat <- mod$coef
          dum <- is.na(beta.hat)
          df <- sum(!dum)
          if (!is.na(S2)) {
            S2 <- (length(use)-df)*S2/rchisq(n = 1, df = (length(use)-df))
          } else {
            S2 <- 0
          }
          newdata <- data.frame(rep(1, length(new)))

          lowupp[[i]] <- FALSE
          dat[new, nam] <- predict(mod, newdata) + rnorm(NROW(newdata), sd = sqrt(S2))
        }
      } else {

        if (is.element(mis[i], bincols)) {
          Y <- as.numeric(as.character(dat[use, nam]) == levels(factor(dat[, nam]))[2])
          alldat <- cbind(rep(1, NROW(dat)))
          coeff <- qnorm(mean(Y))
          means <- c(tcrossprod(alldat, t(coeff)))

          lims <- as.character(dat[, nam])
          upper <- rep(Inf, length(lims))
          lower <- rep(-Inf, length(lims))
          upper[as.character(lims) == levels(factor(lims))[1] & !is.na(lims)] <- 0
          lower[as.character(lims) == levels(factor(lims))[2] & !is.na(lims)] <- 0
          lowupp[[i]] <- cbind(lower = lower, upper = upper)

          dat[, nam] <- truncnorm::rtruncnorm(NROW(dat), a = lower, b = upper, mean = means, sd = 1)
        } else if (is.element(mis[i], censcols)) {
          use.c <- obs[, nam] == 0
          upper <- rep(Inf, sum(use.c))
          lower <- cens[use.c, nam]
          lower[is.na(lower)] <- -Inf
          lowupp[[i]] <- cbind(lower = lower, upper = upper)

          dat[use.c, nam] <- truncnorm::rtruncnorm(sum(use.c), a = lower, b = upper, mean = 0, sd = 1)
        } else if (is.element(mis[i], ordcols)) {
          Y <- ordered(dat[, nam])
          means <- rep(0, NROW(dat))
          zetas <- table(Y[use])
          levs <- names(zetas)
          zetas <- cumsum(zetas)
          zetas <- zetas[-length(zetas)]/zetas[length(zetas)]
          zetas <- c(-Inf, qnorm(zetas), Inf)
          cutpts[[which(ordcols == mis[i])[1]]] <- zetas
          zetas.l <- zetas[-length(zetas)]
          zetas.u <- zetas[-1]

          lower <- zetas.l[match(Y, levs)]
          upper <- zetas.u[match(Y, levs)]
          lower[is.na(lower)] <- -Inf
          upper[is.na(upper)] <- Inf
          lowupp[[i]] <- cbind(lower = lower, upper = upper)

          dat[, nam] <- truncnorm::rtruncnorm(NROW(dat), a = lower, b = upper, mean = means, sd = 1)
        } else {
          olddat <- data.frame(Y = dat[use, nam], dat[use, keep, drop = FALSE])
          mod <- lm(Y~., data = olddat)
          S2 <- (summary(mod)$sigma)^2
          beta.hat <- mod$coef
          dum <- is.na(beta.hat)
          df <- sum(!dum)
          if (!is.na(S2)) {
            S2 <- (NROW(olddat)-df)*S2/rchisq(n = 1, df = (NROW(olddat) - df))
          } else {
            S2 <- 0
          }
          newdata <- data.frame(dat[new, keep, drop = FALSE])

          lowupp[[i]] <- FALSE
          dat[new, nam] <- predict(mod, newdata) + rnorm(NROW(newdata), sd = sqrt(S2))
        }
      }
    }

    if (length(ordcols) > 0) {
      names(cutpts) == ordcols
    }

    return(list(dat = dat, predictorMatrix = predictorMatrix, cutpts = cutpts, drops = drops, lowupp = lowupp))
  }

  ### Define IStep within gerbil.inner
  IStep <- function (data, obs, params, n = nrow(data), p1 = ncol(data), mis = NULL, p = length(mis), olddat = NULL, lowupp = NULL) {

    ### Define getEX within IStep
    getEX <- function (data, betas, mis, n = nrow(data), p1 = ncol(data), p = length(mis)) {
      means <- as.matrix(cbind(rep(1, n), data))
      mis <- mis + 1
      for (i in 1:p) {
        #means[, mis[i]] <- means %*% as.matrix(betas[i, ])
        #means[, mis[i]] <- means %*% t(betas[i, , drop = FALSE])
        means[, mis[i]] <- tcrossprod(means, betas[i, , drop = FALSE])
      }
      return(means[, mis, drop = FALSE])
    }

    ### Define getCovs within IStep
    getCovs <- function (betas, sigmas, p = length(sigmas)) {
      ### betas is a p x p matrix with 0's in the diagonal and the upper triangle.
      covars <- matrix(0, p, p)

      for (j in 1:p) {
        if (j == 1) {
          covars[j, j] <- sigmas[j]
        }
        if (j > 1) {
          #coefs <- matrix(betas[j, 1:(j-1)], 1, (j-1))
          #covars[j, 1:(j-1)] <- coefs %*% matrix(covars[1:(j-1), 1:(j-1)], (j-1), (j-1))
          coefs <- betas[j, 1:(j-1), drop = FALSE]
          #covars[j, 1:(j-1)] <- coefs %*% covars[1:(j-1), 1:(j-1), drop = FALSE]
          covars[j, 1:(j-1)] <- tcrossprod(coefs, t(covars[1:(j-1), 1:(j-1), drop = FALSE]))
          covars[1:(j-1), j] <- covars[j, 1:(j-1)]
          #covars[j, j] <- sigmas[j] + covars[1:(j-1), 1:(j-1)]*(t(coefs) %*% coefs)
          #covars[j, j] <- sigmas[j] + coefs %*% covars[1:(j-1), 1:(j-1)] %*% t(coefs)
          #covars[j, j] <- sigmas[j] + coefs %*% tcrossprod(covars[1:(j-1), 1:(j-1)], coefs)
          covars[j, j] <- sigmas[j] + crossprod(t(coefs), tcrossprod(covars[1:(j-1), 1:(j-1)], coefs))
        }
      }
      return(covars)
    }

    ### Define sweepimpstep wtihin IStep
    sweepimpstep <- function (data, x, obs, means, invSigma, n = nrow(data), p = ncol(data), lims = NULL) {
      ## data is n x p.
      ## means is n x p.
      ## covar is p x p.
      ## obs is n x 1.
      ## x is in (1:p)

      myrswp <- function (V, b) {
        p <- ncol(V)
        u <- is.na(match(1:p, b))
        a <- (1:p)[u]
        out <- 0 * V
        dimnames(out) <- dimnames(V)
        if (length(a) == 0)
          return(-solve(V))
        else if (length(a) == p)
          return(V)
        else {
          Saa <- V[a, a, drop = FALSE]
          Sab <- V[a, b, drop = FALSE]
          Sbb <- V[b, b, drop = FALSE]
          Sbb <- solve(Sbb)
          #B <- Sab %*% (Sbb)
          B <- crossprod(t(Sab),Sbb)
          #out[a, a] <- Saa - B %*% t(Sab)
          out[a, a] <- Saa - tcrossprod(B, Sab)
          out[a, b] <- -B
          out[b, a] <- -t(B)
          out[b, b] <- -(Sbb)
          return(out)
        }
      }

      obs0 <- obs == 0

      imp <- sum(I(obs0))
      y <- setdiff((1:p), x)
      mu1 <- means[, x]
      mu1 <- mu1[obs0]
      invSigma <- (-1) * invSigma
      invSigma <- myrswp(invSigma, x)
      VarX1 <- invSigma[x, x]

      if (length(y) > 0) {
        #mu2 <- means[,y]
        #mu2 <- matrix(mu2[obs == 0], ncol = p - 1)
        #c2 <- as.matrix(data[, y])
        #c2 <- matrix(c2[obs == 0], ncol  =p - 1)
        mu2 <- means[obs0, y, drop = FALSE]
        c2 <- data[obs0, y, drop = FALSE]

        #EX1 <- mu1 + (c2 - mu2) %*% invSigma[y, x]
        EX1 <- mu1 + tcrossprod(as.matrix(c2 - mu2), t(invSigma[y, x]))
      } else {
        EX1 <- mu1
      }

      if(!is.logical(lims)) {
        output <- truncnorm::rtruncnorm(imp, a = lims[, 1], b = lims[, 2], mean = EX1, sd = sqrt(VarX1))
      } else {
        output <- rnorm(imp, mean = EX1, sd = sqrt(VarX1))
      }

      return(output)
    }

    if (length(mis) == 0) {
      mis <- which(colSums(obs == 0) > 0)
      p <- length(mis)
    }

    betas <- params[[1]]
    ### betas is a p x p1 + 1 matrix where row i gives the regression coefficients for the mis[i]^th variable
    sigmas <- params[[2]]
    means <- getEX(data ,betas, mis)
    covar <- getCovs(betas[, (mis + 1)], sigmas ,p)
    covar <- solve(covar)
    for (j in 1:p) {
      i <- mis[j]
      data[,i][obs[,i] == 0] <- sweepimpstep(data = data[, mis], x = j, obs = obs[,i], means = means, invSigma = covar, n = n, p = p, lims = lowupp[[j]])
    }
    return(data)
  }

  ### Define PStep within gerbil.inner
  PStep <- function(dat, predictorMatrix, p = NCOL(dat), n = NROW(dat), mis = 1:p, bincols = NULL, ordcols = NULL) {

    myswp <- function (V, b) {
      p <- ncol(V)
      u <- is.na(match(1:p, b))
      a <- (1:p)[u]
      out <- 0 * V
      dimnames(out) <- dimnames(V)
      if (length(a) == 0)
        return(-solve(V))
      else if (length(a) == p)
        return(V)
      else {
        Saa <- V[a, a, drop = FALSE]
        Sab <- V[a, b, drop = FALSE]
        Sbb <- V[b, b, drop = FALSE]
        Sbb=solve(Sbb)
        #B <- Sab %*% Sbb
        B <- crossprod(t(Sab),Sbb)
        #out[a, a] <- Saa - B %*% t(Sab)
        out[a, a] <- Saa - tcrossprod(B, Sab)
        out[a, b] <- B
        out[b, a] <- t(B)
        out[b, b] <- -Sbb
        return(out)
      }
    }

    myrswp <- function (V, b) {
      p <- ncol(V)
      u <- is.na(match(1:p, b))
      a <- (1:p)[u]
      out <- 0 * V
      dimnames(out) <- dimnames(V)
      if (length(a) == 0)
        return(-solve(V))
      else if (length(a) == p)
        return(V)
      else {
        Saa <- V[a, a, drop = FALSE]
        Sab <- V[a, b, drop = FALSE]
        Sbb <- V[b, b, drop = FALSE]
        Sbb <- solve(Sbb)
        #B <- Sab %*% (Sbb)
        B <- crossprod(t(Sab),Sbb)
        #out[a, a] <- Saa - B %*% t(Sab)
        out[a, a] <- Saa - tcrossprod(B, Sab)
        out[a, b] <- -B
        out[b, a] <- -t(B)
        out[b, b] <- -(Sbb)
        return(out)
      }
    }

    X <- as.matrix(cbind(1, dat))
    #X <- cbind(1,dat)
    X <- crossprod(X)

    betas <- matrix(0, p, p + 1)
    sigmas <- rep(0, p)
    colnames(betas) <- c("B_0", colnames(dat))
    rownames(betas) <- c(colnames(dat))
    betas <- betas[mis, , drop=FALSE]
    sigmas <- sigmas[mis]

    #never <- which(colSums(predictorMatrix) == 0) + 1
    never <- which(colSums(predictorMatrix[mis, , drop = FALSE]) == 0) + 1

    X <- myswp(X, 1)
    X.int <- X
    start <- setdiff(1:p, mis) + 1
    start <- setdiff(start, never)
    if (length(start) > 0) {
      X <- myswp(X, start)
    }

    for (i in 1:length(mis)) {
      if (mis[i] == 1) {
        rem <- NULL
      } else {
        rem <- predictorMatrix[mis[i], 1:(mis[i]-1)]
        rem <- which(rem == 0)
      }

      rem <- rem + 1
      keep <- setdiff(1:mis[i], rem)
      rem <- setdiff(rem, never)
      keep <- setdiff(keep, never)

      if (length(keep) == 1) {
        newX <- X.int
      } else if ((length(keep) - 1) < length(rem)) {
        newX <- myswp(X.int, keep[-1])
      } else {
        newX <- myrswp(X, rem)
      }

      Xinv <- (-1)*newX[keep, keep, drop = FALSE]
      df <- n - length(keep)
      beta.hat <- newX[mis[i] + 1,keep]
      #if (!is.element(mis[i], bincols) & !is.element(mis[i], ordcols)) {
      if (!is.element(mis[i], bincols)) {
        S2 <- newX[mis[i] + 1,mis[i] + 1]/df
        S2 <- df * S2/rchisq(n = 1, df = df)
      } else {
        S2 <- 1
      }
      beta.hat <- mvtnorm::rmvnorm(1, beta.hat, S2*Xinv)
      betas[i, keep] <- beta.hat
      sigmas[i] <- S2
      if (!is.element(mis[i] + 1, never)) {
        X <- myswp(X, mis[i] + 1)
      }
    }

    return(list(betas, sigmas))
  }

  ### Define get.means within gerbil.inner
  get.means <- function(dat, obs, mis, bincols, cols, is.cat, ordcols, cutpts, semicols = NULL, calc.var = FALSE, calc.obs = FALSE, calc.latent = FALSE) {

    colVars <- function(x, na.rm = FALSE) {
      N <- colSums(!is.na(x))
      sumsq <- colSums(x^2, na.rm = na.rm)
      sums <- colSums(x, na.rm = na.rm)
      return((sumsq - sums^2/N)/(N - 1))
    }

    dat[obs == 4] <- NA

    vars.mis.latent <- vars.all.latent <- vars.mis <- vars.all <- NULL
    means.mis.latent <- means.all.latent <- NULL

    if(calc.latent) {
      means.all.latent <- colMeans(dat[, mis, drop = FALSE], na.rm = TRUE)
      dat1 <- dat
      dat1[obs == 0 | obs >= 3] <- NA
      means.mis.latent <- colMeans(dat1[, mis, drop = FALSE], na.rm = TRUE)

      if(calc.var) {
        vars.all.latent <- colVars(dat[, mis, drop = FALSE], na.rm = TRUE)
        vars.mis.latent <- colVars(dat1[, mis, drop = FALSE], na.rm = TRUE)
      }

      non.cont <- c(bincols, ordcols, semicols)
      non.cont <- intersect(non.cont, mis)
      non.cont <- colnames(dat)[non.cont]
    }

    dat[, bincols] <- dat[, bincols] >= 0

    if (length(cols) > 0) {
      for (i in 1:length(cols)) {
        cols.tmp <- cols[[i]]
        if (is.element(i, is.cat) & sum(is.element(cols.tmp, mis) > 0)) {
          dat.tmp <- dat[, cols.tmp]
          dat.tmp[obs[, cols.tmp] >= 3] <- 0
          here <- obs[, cols.tmp[1]] == 1
          if (sum(here) > 0) {
            dat.tmp1 <- dat.tmp[here, ]
            dat.tmp2 <- matrix(0, NROW(dat.tmp1), NCOL(dat.tmp1))
            colnames(dat.tmp2) <- colnames(dat.tmp1)
            dat.tmp2[, 1] <- dat.tmp1[, 1]
            is.last <- dat.tmp2[, 1] == 0
            for (j in 2:NCOL(dat.tmp2)) {
              dat.tmp2[is.last, j] <- dat.tmp1[is.last, j]
              is.last[dat.tmp1[, j] == 1] <- FALSE
            }
            dat.tmp[here, ] <- dat.tmp2
            dat[, cols.tmp] <- dat.tmp
          }
        } else if (sum(is.element(cols.tmp, mis) > 0)) {
          dat[dat[, cols.tmp[1]] == 0, cols.tmp[2]] <- NA
        }
      }
    }

    new.ords <- intersect(mis, ordcols)
    if (length(new.ords) > 0) {
      for(i in new.ords) {
        cutpts.tmp <- cutpts[[which(ordcols == i)]]
        dat[, i] <- as.numeric(cut(dat[, i], breaks = cutpts.tmp))
      }
    }

    if (calc.latent) {
      means.all <- means.all.latent
      vars.all <- vars.all.latent
      means.mis <- means.mis.latent
      vars.mis <- vars.mis.latent
      if (length(non.cont) > 0) {

        dat1 <- dat
        dat1[obs == 0 | obs >= 3] <- NA

        means.all[non.cont] <- colMeans(dat[, non.cont, drop = FALSE], na.rm = TRUE)
        means.mis[non.cont] <- colMeans(dat1[, non.cont, drop = FALSE], na.rm = TRUE)

        if(calc.var) {
          vars.all[non.cont] <- colVars(dat[, non.cont, drop = FALSE], na.rm = TRUE)
          vars.mis[non.cont] <- colVars(dat1[, non.cont, drop = FALSE], na.rm = TRUE)
        }
      }
    } else {
      dat1 <- dat
      dat1[obs == 0 | obs >= 3] <- NA
      means.all <- colMeans(dat[, mis, drop = FALSE], na.rm = TRUE)
      means.mis <- colMeans(dat1[, mis, drop = FALSE], na.rm = TRUE)
      if(calc.var) {
        vars.all <- colVars(dat[, mis, drop = FALSE], na.rm = TRUE)
        vars.mis <- colVars(dat1[, mis, drop = FALSE], na.rm = TRUE)
      }
    }

    means.obs <- vars.obs <- NULL
    if (calc.obs) {
      dat2 <- dat
      dat2[obs == 1 | obs == 2] <- NA
      means.obs <- colMeans(dat2, na.rm = TRUE)
      if(calc.var) {
        vars.obs <- colVars(dat2, na.rm = TRUE)
      }
    }

    return(list(means.all = means.all, vars.all = vars.all, means.mis = means.mis, vars.mis = vars.mis,
                means.all.latent = means.all.latent, vars.all.latent = vars.all.latent,
                means.mis.latent = means.mis.latent, vars.mis.latent = vars.mis.latent,
                means.obs = means.obs, vars.obs = vars.obs))

  }

  ### Start code for gerbil.inner function
  imp <- X

  if (length(cens) > 0) {
    cens <- as.matrix(data.frame(cens))
  }
  if (length(obs) == 0) {
    obs <- !is.na(data)
  }
  if (length(mis) == 0) {
    mis <- sort(union(which(colSums(obs == 0) > 0), censcols))
  }
  bincols <- intersect(bincols, mis)
  ordcols <- intersect(ordcols, mis)
  censcols <- intersect(censcols, mis)
  nmis <- setdiff(1:NCOL(data), mis)

  if (length(mis) == 0) {
    stop("dat has no missing values")
  }

  colnames(obs) <- colnames(predictorMatrix) <- rownames(predictorMatrix) <- colnames(data)

  new.ord <- c(nmis,mis)

  predictorMatrix[upper.tri(predictorMatrix)] <- t(predictorMatrix)[upper.tri(predictorMatrix)]

  data <- data[, new.ord]
  obs <- obs[, new.ord]
  obs.latent <- obs.latent[, new.ord]
  if (length(cat.cols) > 0) {
    for (i in 1:length(cat.cols)) {
      cat.cols[[i]] <- match(cat.cols[[i]], new.ord)
    }
  }
  if (length(semicols.new) > 0) {
    semicols.new <- semicols.new[2, ]
    semicols.new <- match(semicols.new, new.ord)
  }
  predictorMatrix <- predictorMatrix[new.ord, new.ord]
  predictorMatrix[upper.tri(predictorMatrix)] <- 0

  new.mis <- (length(nmis) + 1):NCOL(data)

  if (length(bincols) > 0) {
    bincols <- which(is.element(1:NCOL(data), bincols)[new.ord])
    obs[,bincols] <- 0
  }
  if (length(ordcols) > 0) {
    ordcols <- which(is.element(1:NCOL(data), ordcols)[new.ord])
    obs[,ordcols] <- 0
  }
  if (length(censcols) > 0) {
    censcols <- which(is.element(1:NCOL(data), censcols)[new.ord])
    obs.tmp <- obs[, censcols]
    obs.tmp[cens] <- 0
    cens.new <- matrix(NA,NROW(cens),NCOL(cens))
    cens.new[cens] <- data[, censcols][cens]
    obs[, censcols] <- obs.tmp
    colnames(cens.new) <- colnames(cens)
    #data[, enscols] <- NA
  } else {
    cens.new <- NULL
  }

  first <- proc.time()
  imputes <- gerbil.init(data, obs, predictorMatrix, mis = new.mis, bincols = bincols, ordcols = ordcols, censcols = censcols, cens = cens.new, imp = X, verbose = verbose, r = r)
  lowupp <- imputes[[5]]
  drops <- imputes[[4]]
  cutpts <- imputes[[3]]
  predictorMatrix <- imputes[[2]]
  imputes <- imputes[[1]]
  params <- PStep(imputes, predictorMatrix, mis = new.mis, bincols = bincols, ordcols = ordcols)
  last <- proc.time()
  if (verbose) {
    message("Imp. ",imp,": gerbil initialized.  Time = ", formatC(last[3] - first[3], format = "f", digits = 2), "\n", sep = "", appendLF = FALSE)
  }

  if (trace) {
    means.tmp <- get.means(imputes, obs = obs.latent, mis = new.mis, bincols = bincols, cols = cat.cols, is.cat = is.cat, ordcols = ordcols, cutpts = cutpts, semicols = semicols.new, calc.var = calc.var, calc.obs = TRUE, calc.latent = calc.latent)

    means.all <- means.mis <- means.all.latent <- means.mis.latent <- matrix(NA, length(new.mis), mcmciter + 1)
    dimnames(means.all) <- dimnames(means.mis) <- dimnames(means.all.latent) <- dimnames(means.mis.latent) <- list(colnames(data)[new.mis], c(0, 1:mcmciter))

    means.obs <- means.tmp$means.obs
    vars.obs <- means.tmp$vars.obs

    means.all[, 1] <- means.tmp$means.all
    means.mis[, 1] <- means.tmp$means.mis
    if (length(means.tmp$means.mis.latent) > 0) {
      means.all.latent[, 1] <- means.tmp$means.all.latent
      means.mis.latent[, 1] <- means.tmp$means.mis.latent
    }

    if (calc.var) {
      vars.all <- vars.mis <- vars.all.latent <- vars.mis.latent <- matrix(NA, length(new.mis), mcmciter + 1)
      dimnames(vars.all) <- dimnames(vars.mis) <- dimnames(vars.all.latent) <- dimnames(vars.mis.latent) <- list(colnames(data)[new.mis], c(0, 1:mcmciter))
      vars.all[, 1] <- means.tmp$vars.all
      vars.mis[, 1] <- means.tmp$vars.mis
      if (length(means.tmp$vars.mis.latent) > 0) {
        vars.all.latent[, 1] <- means.tmp$vars.all.latent
        vars.mis.latent[, 1] <- means.tmp$vars.mis.latent
      }
    } else {
      vars.all <- vars.mis <- vars.all.latent <- vars.mis.latent <- NULL
    }
  } else {
    means.obs <- means.all <- means.mis <- means.all.latent <- means.mis.latent <- NULL
    vars.obs <- vars.all <- vars.mis <- vars.all.latent <- vars.mis.latent <- NULL
  }

  for (i in 1:mcmciter) {

    first <- proc.time()
    imputes <- IStep(imputes, obs, params, mis = new.mis, olddat = data, lowupp = lowupp)
    middle <- proc.time()
    params <- PStep(dat = imputes, predictorMatrix = predictorMatrix, mis = new.mis, bincols = bincols, ordcols = ordcols)
    last <- proc.time()
    if (verbose) {
      message("Imp. ", imp, ": MCMC iteration ", i, " completed. Total time = ", formatC(last[3] - first[3], format = "f", digits = 2),", I-Step: ", formatC(middle[3] - first[3], format = "f", digits = 2), ", P-Step: ", formatC(last[3] - middle[3], format = "f", digits = 2), "\n", sep = "", appendLF = FALSE)
    }

    if (trace) {
      means.tmp <- get.means(imputes, obs = obs.latent, mis = new.mis, bincols = bincols, cols = cat.cols, is.cat = is.cat, ordcols = ordcols, cutpts = cutpts, semicols = semicols.new, calc.var = calc.var, calc.obs = FALSE, calc.latent = calc.latent)

      means.all[, i + 1] <- means.tmp$means.all
      means.mis[, i + 1] <- means.tmp$means.mis
      if (length(means.tmp$means.mis.latent) > 0) {
        means.all.latent[, i + 1] <- means.tmp$means.all.latent
        means.mis.latent[, i + 1] <- means.tmp$means.mis.latent
      }

      if (calc.var) {
        vars.all[, i + 1] <- means.tmp$vars.all
        vars.mis[, i + 1] <- means.tmp$vars.mis
        if (length(means.tmp$vars.mis.latent) > 0) {
          vars.all.latent[, i + 1] <- means.tmp$vars.all.latent
          vars.mis.latent[, i + 1] <- means.tmp$vars.mis.latent
        }
      }
    }

  }

  if (length(bincols) > 0) {
    for (i in bincols) {
      levs <- levels(factor(data[, i]))
      is.num <- sum(is.na(levs)) == 0
      tmp <- rep(NA,NROW(imputes))
      tmp[imputes[, i] < 0] <- levs[1]
      tmp[imputes[, i] >= 0] <- levs[2]
      if (is.num) {
        tmp <- as.numeric(as.character(tmp))
      }
      imputes[, i] <- tmp
    }
  }

  if (length(ordcols) > 0) {
    for (i in 1:length(ordcols)) {
      levs <- levels(factor(data[, ordcols[i]]))
      is.num <- sum(is.na(levs)) == 0
      cutpts.tmp <- cutpts[[i]]
      tmp <- as.character(cut(imputes[, ordcols[i]], breaks = cutpts.tmp, labels = levs))
      #tmp <- rep(NA, NROW(imputes))
      #tmp1 <- imputes[, ordcols[i]]
      #for (j in 1:(length(cutpts.tmp) - 1)) {
      #  tmp[tmp1 <= cutpts.tmp[j + 1] & tmp1 > cutpts.tmp[j]] <- levs[j]
      #}
      if (is.num) {
        tmp <- as.numeric(tmp)
      }
      imputes[,ordcols[i]] <- tmp
    }
  }

  imputes <- imputes[, order(new.ord, decreasing = FALSE)]
  predictorMatrix <- predictorMatrix[colnames(data)[new.mis], , drop = FALSE]

  return(list(imputes = imputes, predictorMatrix = predictorMatrix, drops = drops,
              means.obs = means.obs, vars.obs = vars.obs,
              means.all = means.all, vars.all = vars.all,
              means.mis = means.mis, vars.mis = vars.mis,
              means.all.latent = means.all.latent, vars.all.latent = vars.all.latent,
              means.mis.latent = means.mis.latent, vars.mis.latent = vars.mis.latent))

}


gerbil.untrans <- function(olddat, newdat, parms, transcols = 1:NCOL(olddat), transtype = rep("NONE", NCOL(olddat)),
                           nams = (1:NCOL(olddat))[transcols], cens = NULL, verbose = TRUE, X = 1)
{

  tmp.time.all <- proc.time()
  olddat <- as.data.frame(olddat)
  newdat <- as.data.frame(newdat)
  out <- olddat

  for (i in 1:length(transcols)) {
    imped <- is.na(olddat[,transcols[i]])
    if (is.element(colnames(olddat)[transcols[i]],colnames(cens))) {
      imped <- imped|cens[,colnames(olddat)[transcols[i]]]
    }
    #|substr(olddat[,transcols[i]],1,1)==">"
    if (sum(imped)>0) {
      if (transtype[i] == "EMP" | transtype[i] == "censored") {
        tmptime <- proc.time()
        tmp <- newdat[imped,transcols[i]]
        if (transtype[i] == "EMP") {
          tmp <- pnorm(tmp, mean = parms[[i]][[2]][1], sd  =parms[[i]][[2]][2])
        } else {
          tmp <- pnorm(tmp, mean = 0, sd = 1)
        }
        tmp1 <- olddat[, transcols[i]]
        tmp <- qemp(dat = tmp, tab = parms[[i]][[1]])
        tmptime <- proc.time() - tmptime
        #cat("Untransformed ", nams[i], " w/ Emp. Dist., Time = ", formatC(tmptime[3], format = "f", digits = 2), ", m = ", length(tmp), "\n", sep = "")
        if (class(tmp) == "factor") {
          tmp <- as.character(tmp)
          tmp1[imped] <- tmp
          tmp1 <- factor(tmp1)
        } else {
          tmp1[imped] <- tmp
        }
        out[,transcols[i]] <- tmp1
      } else {
        out[,transcols[i]] <- newdat[,transcols[i]]
      }
    }
  }

  non.trans <- setdiff(1:NCOL(olddat), transcols)
  if (length(non.trans) > 0) {
    out[, non.trans] <- newdat[, non.trans]
  }

  tmp.time.all <- proc.time() - tmp.time.all
  if (verbose) {
    message("Completed untransformations for imputed dataset ", X, ", Time = ", formatC(tmp.time.all[3], format = "f", digits = 2), "\n", sep="", appendLF = FALSE)
  }
  return(out)

}


pemp <- function(dat, n = length(dat))
{

  tab <- table(dat)
  ctab <- cumsum(tab)
  ctab1 <- c(0, ctab[-length(tab)])
  ctab2 <- (.5 * tab + ctab1)/n
  out <- ctab2[match(as.character(dat), names(ctab2))]
  out <- unname(out)

  return(list(out, cumsum(tab)/n))

}


pemp.cens <- function(dat, cens = rep(FALSE, NROW(dat)))
{

  #cens <- substr(dat,1,1) == ">"
  #dat <- as.numeric(gsub(">", "", dat, fixed = TRUE))

  tab.o <- table(factor(dat[!cens]))
  ctab.o <- sum(!cens) - cumsum(tab.o)
  vals.o <- as.numeric(names(tab.o))
  ctab.o <- c('-Inf' = sum(tab.o), ctab.o)

  tab.c <- table(factor(dat[cens]))
  #ctab.c <- sum(cens) - cumsum(tab.c)
  ctab.c <- cumsum(tab.c[length(tab.c):1])[length(tab.c):1]
  vals.c <- as.numeric(names(tab.c))

  cens.val <- c(-Inf, vals.c, Inf)

  out.o1 <- out.o <- rep(NA, length(tab.o))
  names(out.o1) <- names(out.o) <- names(tab.o)

  out.c <- rep(NA,length(tab.c))
  names(out.c) <- names(tab.c)

  ctab.c <- c(ctab.c, 'Inf' = 0)

  start.p <- 1

  for (i in 2:length(cens.val)) {
    here.o <- as.character(vals.o[vals.o > cens.val[i-1] & vals.o <= cens.val[i]])
    #here.c <- as.character(vals.c[vals.c > cens.val[i-1] & vals.c <= cens.val[i]])
    here.c <- names(ctab.o)[as.numeric(names(ctab.o))<=cens.val[i]]
    here.c <- here.c[length(here.c)]
    here.c1 <- names(ctab.o)[as.numeric(names(ctab.o))<=cens.val[i-1]]
    here.c1 <- here.c1[length(here.c1)]

    n.tmp <- ctab.o[here.c1] + ctab.c[as.character(cens.val[i])]

    if (length(here.o) > 0) {
      out.o[here.o] <- start.p * (.5 * tab.o[here.o] + ctab.o[here.o] + ctab.c[as.character(cens.val[i])])/n.tmp
      out.o1[here.o] <- start.p * (ctab.o[here.o] + ctab.c[as.character(cens.val[i])])/n.tmp
    }
    if (ctab.c[as.character(cens.val[i])]>0) {
      out.c[as.character(cens.val[i])] <- start.p*(ctab.o[here.c] + ctab.c[as.character(cens.val[i])])/n.tmp
      start.p <- out.c[as.character(cens.val[i])]
    }
  }

  dat1 <- dat2 <- dat
  dat1[cens] <- NA
  dat2[!cens] <- NA

  out <- rowSums(cbind(out.o[match(as.character(dat1), names(out.o))], out.c[match(as.character(dat2), names(out.c))]), na.rm = TRUE)
  out <- 1-out
  out <- unname(out)

  out1 <- 1 - out.o1
  if (out1[length(out1)] != 1) {
    out1 <- c(out1, 1)
    names(out1)[length(out1)] <- paste0(">", max(vals.c))
  }

  out <- list(out, out1)
  if (sum(cens) > 0) {
    out[[3]] <- cens
  }

  return(out)

}


qemp <- function(dat,tab)
{

  out <- cut(dat, breaks = c(0, tab), labels = names(tab))
  out <- as.character(out)
  suppressWarnings(out1 <- as.numeric(out))
  if (sum(is.na(out1)) > 0) {
    return(factor(out))
  } else {
    return(out1)
  }

}



get.R.hat <- function(tmp, N = NROW(tmp), m = NCOL(tmp))
{

  colVars <- function(x, na.rm = FALSE) {
    N <- colSums(!is.na(x))
    sumsq <- colSums(x^2, na.rm = na.rm)
    sums <- colSums(x, na.rm = na.rm)
    return((sumsq - sums^2/N)/(N - 1))
  }

  theta.bar.m <- colMeans(tmp, na.rm = TRUE)
  theta.bar.. <- mean(theta.bar.m, na.rm = TRUE)
  B <- N*var(theta.bar.m)
  s.hat <- colVars(tmp, na.rm = TRUE)
  W <- mean(s.hat)

  sigma2 <- (N-1)/N * W+1/N * B

  V.hat <- sigma2 + B/(m * N)

  #sqrt(((N-1)/N*W+1/N*B)/W)

  var.si2 <- var(s.hat)
  cov.si2.xi2 <- cov(s.hat, theta.bar.m^2)
  cov.si2.xi <- cov(s.hat, theta.bar.m)

  var.V <- ((N - 1)/N)^2 * (1/m) * var.si2
  var.V <- var.V + ((m + 1)/(m * N))^2 * (2/(m - 1)) * B^2
  var.V <- var.V  +2 * (m+1) * (N-1)/(m * N^2) * (N/m) * (cov.si2.xi2 - 2 * theta.bar.. * cov.si2.xi)
  df <- 2 * V.hat^2/var.V

  return(sqrt(V.hat/W * df/(df - 2)))

}


gerbilCluster <- function(n) {
  #requireNamespace("parallel", quietly = TRUE)
  chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  chk <- as.numeric(as.character(chk))
  if (is.logical(n)) {
    if (n) {
      n1 <- parallel::detectCores()
    } else {
      n1 <- 1
    }
  } else if (length(n) == 0) {
    #n1 <- parallel::detectCores() - 1
    n1 <- parallel::detectCores()
    n1 <- floor((n1+1)/2)
  } else {
    n1 <- n
  }
  return(min(n1, chk, na.rm = TRUE))
}
