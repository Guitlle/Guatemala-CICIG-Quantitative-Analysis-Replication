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
# These are the covariates selected in the Crisisgroup analysis
X1 = trends[,c("gdp_per_capita_ppp_2011", 
              "homicide",
              "poverty_headcount_320", 
              "trend_homicide")]
# I have put all covariates and remove those that had collinearity issues
# to have an alternative analysis.
X2 =trends[,c("trend_gdp_per_capita_ppp_2011", 
    "trend_homicide", 
    "trend_household_consumption", 
    "trend_poverty_headcount_320",
    #"adult_literacy_rate", 
    #"gdp_per_capita_ppp_2011", 
    "homicide" 
    #"household_consumption", 
    #"poverty_headcount_320", 
    #"under5_mortality_rate",
    #"youth_literacy_rate"
)
]
# Looking for collinearity
cor(X2)
#X$gdp_per_capita_ppp_2011 = X$gdp_per_capita_ppp_2011/1000

# Without scaling, the optimization process in the entropy balance method
# cannot converge well. Why Stata optimizer is ok without resclaing?
X1 = scale(X1)
X2 = scale(X2)

synthCtrl1 <- ebalance2(tx, X1, print.level = 3, constraint.tolerance = 0.05, norm.constant = 1)
names(synthCtrl1$w) = as.character(trends$Country.Code[tx==0])
round(synthCtrl1$w/sum(synthCtrl1$w), 4)
sum(synthCtrl1$w)

write.csv(synthCtrl1$w/sum(synthCtrl1$w), "ebalanced.csv")
# Stata allows to specify targets explicitly, while R method does not. The crisis group
# analysis used the option "target(1)" in their Stata code. 
# R optimizer is always different from Stata's so the weights are not very close.
# Furthermore, the large numbers in one of the columns kept the method from converging.
# After rescaling, and adjusting tolerance, it converged with satisfactory results.
# Despite this, the same pattern arises: El Salvador, Rep Dominicana and 
# El Salvador are the top similar countries).
# COL    CRI    DOM    HND    MEX    NIC    PAN    SLV    VEN 
# 0.0036 0.0329 0.2417 0.0287 0.0001 0.2087 0.0021 0.4823 0.0000 
# Also, the ebalance algorithm is not respecting the normalization constraint.
# Probably because of the lack of the bad convergence and scaling arbitrarity.
# The alternative model below has a bigger weights sum, (1.9) which is probably because 
# it doesn't even converge.
# This could also be due to the very small sample size for this dataset.
# Only 9 controls, 1 intervened case, and a few covariates with means and trends
# Changing the scaling (dividing the covariates by a constant) has certain effect over the output weights.
# However, the same top 3 countries are found: SLV, NIC, DOM. 

# The alternative covariates selections has many convergence problems. 
# This definitely influenced the decision to select the original covariates above. 
synthCtrl2 <- ebalance2(tx, X2, print.level = 3, constraint.tolerance = 1, max.iterations = 20, norm.constant = 1)
names(synthCtrl2$w) = as.character(trends$Country.Code[tx==0])
round(synthCtrl2$w/sum(synthCtrl2$w), 4)
sum(synthCtrl2$w)
# This selection of covariates has the following results:
#    COL    CRI    DOM    HND    MEX    NIC    PAN    SLV    VEN 
# 0.0000 0.0000 0.1081 0.5590 0.0162 0.0407 0.1779 0.0121 0.0861 
# These results are quite different from the original selection. 
# Small sample size could also influence this bad results.
# The original covariates selection shows stable results across two
# different software packages (R and Stata) and show much better results than
# alternative covariates selections, which suggests that they are and optimal selection
# of covariates for this analysis.
# In stata, ebalance complains about the last 2 covariates in X2 and after removing them, it does not converge.

timeseries = read.csv("data_time_series.csv")

timeseries
