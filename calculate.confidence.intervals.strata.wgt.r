calculate.confidence.intervals.strata.wgt = function (species=species, data.yrs) {
#-----------------------------------------------------------------
#methods for calculating confidence intervals for grounfish survey indices
#based on the paper by Stephen Smith (1997) 
#"Bootstrap confidence limits for groundfish trawl survey estimates of mean abundance. Can. J. Fish. Aquat. Sci.
#This script uses the mirror match bootstrap method (BWR), described in the paper, which was shown to perform better on average with respect to accuracy of the confidence
#intervals being estimated. And the Bias Corrected (BC) method of calculating confidence intervals is used.  
#Three methods of bootstrapping and calculating confidence intervals were compared in the paper. 
print("correct")  
#Setup
p = setup.parameters()
loadfunctions()

#Extract data from groundfish survey
#Extract Cat data from groundfish survey
k = groundfish.db(DS="cat.base") #export from grounfish survey database

j = survey.process(k, species, data.yrs)
#Gulf Strata
j = subset(j, !(strat == '558' | strat == '559'))
#Georges Bank Strata
j = subset(j, !(strat == "5Z1"| strat == "5Z2"| strat == "5Z3"| strat == "5Z4"| strat == "5Z8"| strat == "5Z7" |
                  strat == "5Z6"| strat == "5Z5"| strat == "5Z9" ))
#Spring Strata
j = subset(j, !(strat == "406" | strat == "401" | strat == "402" | strat == "400" | strat == "397" | strat == "399" 
                | strat == "398" | strat ==  "404" | strat == "403"| strat ==  "407" | strat == "405" | strat == "408"
                | strat == "409" | strat ==  "410" | strat == "411" ))
#Deep Water Strata
j = subset(j, !(strat =="503" | strat == "501" | strat == "502" | strat == "505" | strat == "504" | strat == "558" | strat == "559"))
j = subset(j, !(strat =="498" |strat =="497" |strat == "496"))

#Test
#j = j[which(j$yr == 2006), ]

j$Strata = j$strat
species = j$name.common[1]

#choose variable name
#variable.n = "totno_sd"
variable.n = "totwgt_sd"

#Load file from groundfish.stratum where area is converted to trawlable units divide area by 0.011801 (41ft by 1.75 nm)
#Already converted to list NH = trawalable units, Strata = strata
#Extract Strata Area Information
load("C:\\ss\\data\\strat.rdata")

sc.t <- j
st = strat1

#Calculate stratified estimates from the groundfish survey
#---------------------------------------------------------------------------

ci.n = data.frame()
ci.w = data.frame()
yrs = data.yrs
alpha.t = 0.05 #confidence interval eg. 0.05 = 95%, 0.1 = 90%
alpha.b = 0.05
nresamp = 1000
prints = T
method = "BWR"
CI.method = 'BC'
m <- match.call()

for (yr in yrs){
  #subset by year
  sc = sc.t[which(sc.t$yr == yr), ]
  
  #Choose Variable
  #variable = sc$totno_sd
  variable = sc$totwgt_sd
  
  s.group <- is.element(st$Strata, sc$Strata)
  s.group.Strata <- st$Strata[s.group]
  s.group.NH <- st$NH[s.group]
  s.obj <- is.element(sc$Strata, s.group.Strata)
  var <- variable[s.obj]
  sc <- sc[s.obj, ]

  yhi <- split(var, sc$Strata) #split the variable by strata
  nh <- as.vector(sapply(yhi, length)) #numer of tows per strata
  nhws <- sapply(yhi, function(x) length(x [x > 0])) #calculate the number of samples > 1 in each strata
  
  Strata = s.group.Strata #list of strata
  Nh = s.group.NH #trawalable units in each strata
  Wh = s.group.NH/sum(s.group.NH) #strata percent of the total area in trawlable units
  Sh = sum(s.group.NH)
  call = m
  
  #Calculate Stratified Estimates and Confidence intervals 
  #-------------------------------------------------------------------------------------
  yh <- as.vector(sapply(yhi, mean)) #mean of variable for each strata
  yst <- sum(Wh * yh, na.rm = TRUE) #sum of the mean of the variable for each strata, multiplied by percent area of each strata

  sh <- as.vector(sapply(yhi, var)) #calculate variance of each variable, per strata
  se.yst <- sqrt(sum((((Nh * (Nh - nh))/sum(Nh)^2) * sh)/nh, na.rm = TRUE)) #calculate standard error
  ah <- (Nh * (Nh - nh))/nh #
  df.yst <- (sum(ah * sh, na.rm = TRUE)^2)/(sum(((ah * sh)^2)/(nh - 1), na.rm = TRUE)) #degrees of freedom

  #Calculate the confidence interval, based on the t distribution, rather than the normal distribution
  #qt fuction calculates the value for the 95% confidence interval by looking up the t distribution,
  #based on the degrees of freedom calculated above. This is multiplied by the standard error calculated above
  #formulas are: Lower limit = M - (tCL)(sM), Upper limit = M + (tCL)(sM)
  #------------------------------------------------------------------------------------
  ci.yst <- yst + (c(qt(alpha.t/2, df.yst), -qt(alpha.t/2, df.yst)) * se.yst) #confidence interval
  
  #Calculate Design Weighted Area Occupied
  dwao <- sum(Wh*(nhws / nh)) * sum(Nh) * 0.011801
  
  #Calculate Gini Index
  gi=NA
  
  gi <- gini(x=yh,y=Nh)
  
  #Calculate the total population
  Yst = yst * sum(Nh) 
 
  #Use mirror match method (BWR) to calculate confidence intervals in Bias Corrected (BC) method to calculate confidence interval
  #-------------------------------------------------------------------------
  call <- match.call(expand = FALSE)
  
  out <- matrix(0, nrow = nresamp + 1, ncol = 3, dimnames = list(c("Actual", 1:nresamp), c("Mean", "Variance",'gini')))
  out[1, ] <- c(yst, (se.yst)^2, gi)
  
  
  fh <- nh/Nh
  kh <- (nh - 1)/(1 - fh)
  ph <- ((1/kh) - (1/ceiling(kh)))/((1/floor(kh)) - (1/ceiling(kh)))
  
  for (i in 1:nresamp) {
      yhib <- bwr.boot(yhi, kh, ph, sample, replace = TRUE, simplify = FALSE)
      yhib[nh == 1] <- yhi[nh == 1]
      nhws = sapply(yhib, FUN = function(x) sum(x > 0))
      out[i + 1, ] <- c(sum(Wh * as.vector(sapply(yhib, mean)), na.rm = TRUE),
                      (sum((((Nh * (Nh - as.vector(sapply(yhib, length))))/sum(Nh)^2) *
                      (as.vector(sapply(yhib, var))))/as.vector(sapply(yhib, length)), na.rm = TRUE)),
                       gini(x = as.vector(sapply(yhib, mean)), y = Nh))
    }
  
  
  orig.mean = out[1, 1]
  orig.var = out[1, 2]
  boot.means = out[c(2:(nresamp + 1)), 1]
  boot.vars = out[c(2:(nresamp + 1)), 2]
  gi = out[c(2:(nresamp + 1)), 3]
  call = call
  method = method
  
  #Summary of Bootstrapped means
  #-------------------------------------------------------------------------------------------------------
  options(digits = 4)
  boot.est <- mean(boot.means)
  gini.mean <- mean(gi)
  ci.boot=list()
  
  loc.bc <- sum(boot.means < boot.est)
  lim.bc <- sort(boot.means)[c(loc.bc, loc.bc + 1)]
  z0 <- (loc.bc + ((boot.est - lim.bc[1])/(lim.bc[2] - lim.bc[1])))
  z0 <- qnorm(z0/length(boot.means))
  probs.z0 <- pnorm(qnorm(c(alpha.b/2, (1 - alpha.b/2), 0.5)) + 2 * z0)
  ci.boot[[1]] <- ci.boot.mean <- quantile(boot.means, probs = probs.z0)
  ci.boot.gini <- quantile(gi, probs = c(alpha.b/2, (1 - alpha.b/2), 0.5), na.rm=T)
  
    
  #Print out the yearly estimates and write them to a data frame  
  #--------------------------------------------------------------------------------------------------------
  options(digits = max(options()$digits - 5, 5))
  
  
  if(prints) {cat("\n", "Pop Total =", format(Yst), "\n", 
                    "Original Mean =", format(orig.mean), "\n",
                    "Year =", yr, "\n",
                    "Original Variance =", format(orig.var), "\n",
                    "Number of bootstraps = ", length(boot.means), "\n",
                    "Bootstrap Mean=", format(boot.est), "\n",
                    "Variance of Bootstrap Mean=", format(var(boot.means)), "\n",
                    "CI Method=", c(CI.method), "\n",
                    "CI's for alpha=", alpha.b, "are ", format(ci.boot.mean[1:2]), "\n",
                    "Length =", format(ci.boot.mean[2] - ci.boot.mean[1]), "\n",
                    "Shape=", format(log((ci.boot.mean[2] - ci.boot.mean[3])/(ci.boot.mean[3] - ci.boot.mean[1]))), "\n",
                    "Resample Method = ", method, "\n")
    }
  
  ci.boot[[2]] = ci.boot.gini
  
  n.row <- data.frame(species,
            year = as.numeric(yr),
            pop.total = Yst,
            variable = variable.n,
            orig.mean = as.numeric(format(orig.mean)),
            boot.mean = as.numeric(format(boot.est)),
            var.boot.mean = as.numeric(format(var(boot.means))),
            lower.ci = as.numeric(format(ci.boot.mean[1])),
            upper.ci = as.numeric(format(ci.boot.mean[2])),
            length = as.numeric(format(ci.boot.mean[2] - ci.boot.mean[1])),
            dwao = dwao,
            gini = gini.mean,
            lower.ci.gini = format(ci.boot.gini[1]),
            upper.ci.gini = format(ci.boot.gini[2]))
  

  ci.w = rbind(n.row, ci.w )

}

print(ci.w)

return(ci.w)

}









