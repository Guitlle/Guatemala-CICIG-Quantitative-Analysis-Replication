# 30/01/2019 
# 
# Microsynth analysis for synthetic counterfactual estimation with more Worldbank data.
library(microsynth)
library(readxl)
library(stringr)
library(data.table)

wbdata = read.csv("../../DATOS/Worldbank/WDI/WDIData.csv")

homicides = "VC.IHR.PSRC.P5"
data_sub = wbdata[wbdata$Indicator.Code == homicides & wbdata$Country.Code %in% c(
    "GTM", "SLV", "MEX", "HND", "NIC", "CRI", "PAN", "VEN", "COL", "ECU", "HTI", "BOL", "BLZ", "PER", "BRA", "ARG", "GUY", "CHL", "URY"
),]
colns = colnames(data_sub) 
yearsc = str_sub(colns, 1,1) == "X" & colns != "X"
names(data_sub) = ifelse(yearsc, str_sub(colns, 2,5), colns)
data_sub = data.table(data_sub)
data_wide = melt(data_sub, id.vars = c("Country.Code", "Indicator.Code"), measure.vars = names(data_sub)[yearsc])
data_wide[, year := as.numeric(as.character(variable))]
data_wide = data_wide[year>1996 & year<2017,]
data_wide[, mean_homicide := mean(value, na.rm = T), by=Country.Code]
data_wide[is.na(value), value := mean_homicide]
data_wide[, interv := as.numeric(Country.Code == "GTM" & year >= 2007)]
data_wide[, homicide_rate := value]

msresult = microsynth(data = data_wide, idvar = "Country.Code", 
                      timevar = "year", intvar = "interv", 
                      #match.covar.min = c("mean_gdp", "mean_hc", "mean_povhc", "mean_under5_mm", "mean_youthlr"),
                      match.out = c("homicide_rate"), result.var = c("homicide_rate"),
                      use.backup = T, test = "lower", 
                      perm=250, jack = T)

summary(msresult)
# Trt    Con Pct.Chng Linear.pVal Linear.Lower Linear.Upper Perm.pVal
# homicide_rate 366.73 391.71    -6.4%      0.2399       -19.4%         8.8%    0.2778
# Omnibus           --     --       --      0.2447           --           --    0.2778
round(msresult$w$Weights[,1], 4)
# Main weights to build the synthetic control were:
# ARG    BLZ    BOL    BRA    CHL    COL    CRI    ECU    SLV    GTM    GUY    HTI    HND    MEX    NIC    PAN    PER    URY 
# 0.0000 0.0251 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.1968 1.0000 0.1857 0.0000 0.0000 0.0000 0.1899 0.0000 0.0000 0.0000 
# VEN 
# 0.4026 
diff = data.frame(diff = msresult$Plot.Stats$Difference[1,1,])
diff$year = row.names(diff)

gtmpop_d = merge(gtmpop, diff, by = "year")
sum(gtmpop_d$GTM_pop100k * gtmpop_d$diff)
# -4484.148

