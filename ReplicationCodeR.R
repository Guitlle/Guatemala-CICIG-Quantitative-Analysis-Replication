library(data.table)
library(ebal)
library(matlib)

# After many errors and warnings, fixed the ebal method
# The issue was that this data only have one treatment case
# and when ebalance tries to sum over the columns, it actually
# sums over one only row, generating one number, instead of a vector
source("ebal.fixed.R")

trends = read.csv("trends.csv")

tx = as.integer(trends$Country.Code == "GTM")
X = trends[,c("gdp_per_capita_ppp_2011", 
              "homicide",
              "poverty_headcount_320", 
              "trend_homicide")]
#X$gdp_per_capita_ppp_2011 = X$gdp_per_capita_ppp_2011/1000
X = scale(X)

synthCtrl <- ebalance2(tx, X, print.level = 3, constraint.tolerance = 0.1, max.iterations = 200)
names(synthCtrl$w) = as.character(trends$Country.Code[tx==0])
synthCtrl$w
# Stata allows to specify targets explicitly, while R method does not. The crisis group
# analysis used the option "target(1)" in their Stata code. 
# R optimizer is always different from Stata's so the weights are not very close.
# Furthermore, the large numbers in one of the columns kept the method from converging.
# After rescaling, and adjusting tolerance, it converged with satisfactory results.
# Despite this, the same pattern arises: El Salvador, Rep Dominicana and 
# El Salvador are the top similar countries).
# COL         CRI         DOM         HND         MEX         NIC         PAN         SLV         VEN 
# 0.004456396 0.054981864 0.219866243 0.045199024 0.002299287 0.207508392 0.015020446 0.515068975 0.002768139 