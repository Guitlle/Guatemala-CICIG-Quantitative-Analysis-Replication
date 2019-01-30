# 30/01/2019 
# 
# Microsynth analysis for synthetic counterfactual estimation.
library(microsynth)
library(haven)
library(readxl)

gtmpop = read_excel("gtm_pop.xlsx")
gt1 = read_stata("./guatemala1.dta")
colnames(gt1)
gt1$interv = as.integer(gt1$CountryCode == "GTM" & gt1$year>=2007)
table(gt1$year)
unstack(gt1[, c("homicide", "CountryCode", "year")], form = homicide~ CountryCode)

gt1$log_homicide = log(gt1$homicide + 1)
gt1$gdp_per_capita_ppp_2011 = as.numeric(gt1$gdp_per_capita_ppp_2011)

temp = aggregate(gdp_per_capita_ppp_2011 ~ CountryCode, gt1, mean)
names(temp) = c("CountryCode", "mean_gdp")
gt1 = merge(gt1, temp, by = "CountryCode")

temp = aggregate(household_consumption ~ CountryCode, gt1, mean)
names(temp) = c("CountryCode", "mean_hc")
gt1 = merge(gt1, temp, by = "CountryCode")

temp = aggregate(poverty_headcount_320 ~ CountryCode, gt1, mean)
names(temp) = c("CountryCode", "mean_povhc")
gt1 = merge(gt1, temp, by = "CountryCode")

temp = aggregate(under5_mortality_rate ~ CountryCode, gt1, mean)
names(temp) = c("CountryCode", "mean_under5_mm")
gt1 = merge(gt1, temp, by = "CountryCode")

temp = aggregate(youth_literacy_rate ~ CountryCode, gt1, mean)
names(temp) = c("CountryCode", "mean_youthlr")
gt1 = merge(gt1, temp, by = "CountryCode")

# input
gt2 <- data.frame(gt1)
temp = aggregate(homicide ~ CountryCode, gt2, mean)
gt2$homicide[is.na(gt2$homicide)] = merge(gt2[, c("CountryCode", "CountryName")], temp, all = T, by = "CountryCode")$homicide[is.na(gt2$homicide)]

msresult = microsynth(data = gt2[gt2$year >= 2000 & gt2$year <= 2015, ], idvar = "CountryCode", 
                      timevar = "year", intvar = "interv", 
                      match.covar = c("mean_gdp", "mean_hc", "mean_povhc", "mean_under5_mm", "mean_youthlr"),
                      match.out = c("homicide"), result.var = c("homicide"),
                      test = "lower", use.backup = T ,
                      perm=250, jack = T)
summary(msresult)
# Trt    Con Pct.Chng Linear.pVal Linear.Lower Linear.Upper Perm.pVal
# homicide 343.84 386.33   -11.0%      0.1736       -26.7%         8.1%    0.4444
# Omnibus      --     --       --      0.2106           --           --    0.4444
round(msresult$w$Weights[,1], 4)
# Main weights to build the synthetic control were:
# COL    CRI    DOM    GTM    HND    MEX    NIC    PAN    SLV    VEN 
# 0.0000 0.0000 0.3443 1.0000 0.2989 0.0000 0.1969 0.0000 0.1598 0.0000 
diff = data.frame(diff = msresult$Plot.Stats$Difference[1,1,])
diff$year = row.names(diff)

gtmpop_d = merge(gtmpop, diff, by = "year")
sum(gtmpop_d$GTM_pop100k * gtmpop_d$diff)
# -6805.191
