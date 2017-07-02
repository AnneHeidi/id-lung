
### Comprehensive smoking index
# As proposed in http://onlinelibrary.wiley.com/doi/10.1002/sim.2680/full

# I have a data frame named covs that include all the covariates = questionnaire information. 
# In that data I have made variables that have a categorical smoking status (covs$smoking_status_cat), 
# estimated duration of smoking in years (covs$duration.adj), 
# time in years since quit smoking for former smokers (covs$tsq_smoking),
# and mean intensity across the smoking years in sigarettes per day (covs$cum.mean.int)
# 

covs$caco<-ifelse(covs$case=="ctrl", 0, 1)

results_CSI_AIC <- matrix(NA, (length(seq(0.1, 2.7, by=0.1))*length(11:55)),3,
                          dimnames=list(1:1215, c("Tau", "Delta","AIC")))
f <- "caco ~ CSI + age.sample"

tau_range <- rep(11:55, 27)
delta_range <- rep.int(seq(0.1, 2.7, by=0.1), c(rep(45, 27)))

by(covs$tsq_smoking, covs$smoking_status_cat, summary)
by(covs$cum.mean.int, covs$smoking_status_cat, summary)

dur <- covs$duration.adj
tsc <- covs$tsq_smoking
int <- covs$cum.mean.int 

for (i in 1:(length(seq(0.1, 2.7, by=0.1))*length(11:55))) {
  
  tau<-tau_range[i]
  delta<-delta_range[i]
  results_CSI_AIC[i,"Tau"] <- tau
  results_CSI_AIC[i,"Delta"] <- delta
  
  TSC <- pmax(tsc - delta, 0)
  DUR <- pmax(dur + tsc - delta,0) - TSC
  
  CSI <- ifelse(covs$smoking_status_cat=="Never", 0, (1-0.5^(DUR/tau))*(0.5^(TSC/tau))*(log(int+1)))
  
  
  model <- glm(as.formula(f), family="binomial", data=covs)   #  if by clogit: model <- clogit(covs$case=="case" ~ smoking_status_cat + strata(pairs), method="exact", data=covs)
  results_CSI_AIC[i,"AIC"] <- model$aic                       #results_CSI_AIC[i,"AIC"] <- extractAIC(mod)[2]
  print(i)
}

results_CSI_AIC<-as.data.frame(results_CSI_AIC)
summary(results_CSI_AIC$AIC)

min<-min(results_CSI_AIC$AIC, na.rm=TRUE)
tau <- results_CSI_AIC[results_CSI_AIC$AIC==min,]$Tau # 25
delta <- results_CSI_AIC[results_CSI_AIC$AIC==min,]$Delta #2.7 
TSC <- pmax(tsc - delta, 0)
DUR <- pmax(dur + tsc - delta,0) - TSC
covs$CSI <- ifelse(covs$smoking_status_cat=="Never", 0, (1-0.5^(DUR/tau))*(0.5^(TSC/tau))*(log(int+1)))

summary(covs$CSI)
table(round(covs$CSI, 1), covs$smoking_status_cat)
boxplot(covs$CSI, covs$smoking_status_cat)

##

