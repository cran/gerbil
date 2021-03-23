
#### Data Documentation --------------------------------------------------------------####

#----------------------------------------------------------------------------------------#
# Author: Pedro Nascimento de Lima
# Purpose: This File documents the datasets shipped with the gerbil package
# Creation Date: Sept 2020
#----------------------------------------------------------------------------------------#


#' Example data from the India Human Development Survey
#'
#' This dataset is a subset from the India Human Development survey. This dataset is included in the package only for demonstration purposes and should not be used for other purposes.
#'
#' @format A data frame with 42155 rows and 8 variables:
#' \describe{
#'		\item{sex}{0 = individual is male; 1 = individual is female}
#'		\item{age}{Age of the individual, between 0 & 99}
#'		\item{marital_status}{Individual’s marital status. 0 = Unmarried; 1 = Married}
#'		\item{job_field}{Refer’s to the field of the individual’s profession or job status (i.e. agricultural worker; small business owner; student; unemployed; etc.}
#'		\item{farm_labour_days}{Number of days a year the individual worked on a farm. Can take on the value zero}
#'		\item{own_livestock}{0 = Individual does not own livestock. 1 = Individual does own livestock}
#'		\item{education_level}{Years of schooling attained by the individual. Censored for values above 16.}
#'		\item{income}{Household’s income, in rupees. Value can be negative}
#' }
#' @source
#' Desai, Sonalde, and Vanneman, Reeve. India Human Development Survey-II (IHDS-II), 2011-12.
#' Inter-university Consortium for Political and Social Research [distributor],
#' 2018-08-08.
#' \doi{10.3886/ICPSR36151.v6}
'ihd'


#' Missing at Random example data from the India Human Development Survey
#'
#' This dataset is a subset from the India Human Development survey. This dataset is included in the package only for demonstration purposes and should not be used for other purposes.
#'
#' @format A data frame with 42155 rows and 8 variables:
#' \describe{
#'		\item{sex}{1 = individual is male; 2 = individual is female}
#'		\item{age}{Age of the individual, between 0 & 99}
#'		\item{marital_status}{Individual’s marital status. 0 = married, absent spouse, 1 = Married, 2 = Unmarried, 3 = Widowed, 4 = Divorced/Separated, 5 = married, no gauna}
#'		\item{job_field}{Refer’s to the field of the individual’s profession or job status (i.e. agricultural worker; small business owner; student; unemployed; etc.}
#'		\item{farm_labour_days}{Number of days a year the individual worked on a farm. Can take on the value zero}
#'		\item{own_livestock}{0 = Individual does not own livestock. 1 = Individual does own livestock}
#'		\item{education_level}{Years of schooling attained by the individual. Censored for values above 16.}
#'		\item{income}{Household’s income, in rupees. Value can be negative}
#' }
#' @source
#' Desai, Sonalde, and Vanneman, Reeve. India Human Development Survey-II (IHDS-II), 2011-12.
#' Inter-university Consortium for Political and Social Research [distributor],
#' 2018-08-08.
#' \doi{10.3886/ICPSR36151.v6}
'ihd_mar'


#' Missing Completely at Random example data from the India Human Development Survey
#'
#' This dataset is a subset from the India Human Development survey. This dataset is included in the package only for demonstration purposes and should not be used for other purposes.
#'
#' @format A data frame with 42155 rows and 8 variables:
#' \describe{
#'		\item{sex}{1 = individual is male; 2 = individual is female}
#'		\item{age}{Age of the individual, between 0 & 99}
#'		\item{marital_status}{Individual’s marital status. 0 = married, absent spouse, 1 = Married, 2 = Unmarried, 3 = Widowed, 4 = Divorced/Separated, 5 = married, no gauna}
#'		\item{job_field}{Refer’s to the field of the individual’s profession or job status (i.e. agricultural worker; small business owner; student; unemployed; etc.}
#'		\item{farm_labour_days}{Number of days a year the individual worked on a farm. Can take on the value zero}
#'		\item{own_livestock}{0 = Individual does not own livestock. 1 = Individual does own livestock}
#'		\item{education_level}{Years of schooling attained by the individual. Censored for values above 16.}
#'		\item{income}{Household’s income, in rupees. Value can be negative}
#' }
#' @source
#' Desai, Sonalde, and Vanneman, Reeve. India Human Development Survey-II (IHDS-II), 2011-12.
#' Inter-university Consortium for Political and Social Research [distributor],
#' 2018-08-08.
#' \doi{10.3886/ICPSR36151.v6}
'ihd_mcar'
