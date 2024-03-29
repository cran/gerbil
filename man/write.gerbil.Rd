% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R
\name{write.gerbil}
\alias{write.gerbil}
\title{Write imputed datasets from gerbil objects to a file or files}
\usage{
write.gerbil(gerb, file = NULL, imp = NULL, tall = FALSE, row.names = FALSE)
}
\arguments{
\item{gerb}{A \code{gerbil} object containing the imputed data.}

\item{file}{The name of the file to which the imputed datasets are to be written.
Which type of file (.xlsx or .csv) is created depends upon the extension of the parameter \code{file}.}

\item{imp}{The imputed datasets which are to be written.  Can be a scalar or, if multiple imputed datasets are to be written, a vector.
All elements of \code{imp} should be integers greater than 0 but no greater than \code{m}, which is the number of imputed datasets in \code{gerb}.
\code{imp} defaults to \code{1:m}.}

\item{tall}{A logical expression indicating whether the datasets are to be written in a tall (stacked) format
or written separately.  When writing to an XLSX file with \code{tall = FALSE}, one tab is created for each imputed dataset.
When writing to a CSV file with \code{tall = FALSE}, one file is created for each imputed dataset
(in this case, several file names will be created from the base string given by \code{file}).}

\item{row.names}{A logical value indicating whether the row names of the datasets are to be written.}
}
\value{
No returned value, but instead a data file is written to a specified directory.
}
\description{
Using a \code{gerbil} object as an input, this function writes imputed datasets to an output file.
}
\details{
The function writes imputed datasets to either an Excel (.xlsx) or a CSV (.csv) file, depending upon the extension of the parameter \code{file}.
  No other file types are supported.
  To write multiple imputed datasets simultaneously, specify \code{imp} as a vector with length greater than 1.
  Multiple imputed datasets are either written in a stacked format (if \code{tall = TRUE}) or written separately (if \code{tall = FALSE}).
}
\examples{
\donttest{
#Load the India Human Development Survey-II dataset
data(ihd_mcar)

# Create a gerbil object
imps.gerbil <- gerbil(ihd_mcar, m = 5, ords = "education_level", semi = "farm_labour_days", 
       bincat = "job_field", n.cores = 1)

# Write all imputed datasets to separate CSV files
write.gerbil(imps.gerbil, file.path(tempdir(), "gerbil_example.csv"), imp = 1:5, tall = FALSE)

# Write all imputed datasets to a single CSV files
write.gerbil(imps.gerbil, file.path(tempdir(), "gerbil_example.csv"), imp = 1:5, tall = TRUE)

# Write all imputed datasets to an XLSX file
write.gerbil(imps.gerbil, file.path(tempdir(), "gerbil_example.xlsx"), imp = 1:5, tall = FALSE)
}

}
