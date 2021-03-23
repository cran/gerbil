
#### Automated Tests ---------------------------------------------------------------####

#----------------------------------------------------------------------------------------#
# Author: Pedro Nascimento de Lima
# Purpose: This File contains automated tests applied to the gerbil class.
# Creation Date: Sept 2020
#----------------------------------------------------------------------------------------#

# Test 1 - Low-Dimensional Manually Created Dataset:

# Step 1 - Generate Dummy Test Datasets:
library(dplyr)
library(gerbil)

# Creating a missing data frame:
missing_data = data.frame(
  X1 = rnorm(n = 1000, mean = 10, sd = 3),
  X2 = runif(n = 1000, min = 10, max = 40)
) %>%
  mutate(X3 = X1 * (2 * rnorm(n = 1000)) + X2)

# Create Missing at Random Data:
missing_data$X1[sample(1000, 100)] <- NA
missing_data$X2[sample(1000, 100)] <- NA
missing_data$X3[sample(1000, 100)] <- NA


# Step 2 - Test Gerbil and alternative package:

# test_that("mice can handle this missing data file", {
#
#   library(mice, quietly = T)
#
#   # Testing with MICE:
#   mice_test = mice(data = missing_data, printFlag = F)
#
#   expect_equal(class(mice_test), "mids")
# })


test_that("gerbil can handle this missing data file", {

  gerbil_object = gerbil::gerbil(dat = missing_data)

  expect_equal(class(gerbil_object), "gerbil")

})



### Using the IHD internal dataset

data("ihd_mcar")

#This will load the ihd_ii_miss dataset (the other dataset will be used for diagnostic purposes)

test_that("Gerbil works without types specified", {

  gerbil_object = gerbil(ihd_mcar, m = 1)

  expect_equal(class(gerbil_object), "gerbil")

})

