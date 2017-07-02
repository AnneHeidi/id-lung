##-------------------------------------------------------------------------------------------------------------------------------## 
#          DNA methylation in C-C series for Lung cancer        Part 2: Logistic regression More efficient code
##-------------------------------------------------------------------------------------------------------------------------------## 
# Preprocessed by Gianluca, residuals made by Florence using chip and 6cat chip position
rm(list=ls())

covs <- readRDS("/home/tno010/Desktop/DNAm GE Lung cancer case-control 2016/NOWAC DNAm Therese/Analyser/output/Lung_covs_080316.rds")  # 263 obs, 261 vars

##-------------------------------------------------------------------------------------------------------------------------------## 
### Differential methylation analysis
# residual data set
#DNAmres_set <- readRDS("/home/tno010/Desktop/NOWAC DNAm Therese/Analyser/Lung cancer data set/res_NOWAC_lung_MM_Methyl_techcovar_residuals.rds")
#DNAmres <- DNAmres_set[[2]]

#DNAmres_set <- readRDS("/home/tno010/Desktop/NOWAC DNAm Therese/Analyser/Lung cancer data set/res_NOWAC_lung_MM_Methyl_techcovar_residuals_storagetime.rds")
#DNAmres <- DNAmres_set[[2]]

DNAmres_set <- readRDS("/home/tno010/Desktop/DNAm GE Lung cancer case-control 2016/NOWAC DNAm Therese/Analyser/Lung cancer data set/res_NOWAC_lung_MM_Methyl_techcovar_residuals_storagetime_age.rds")
DNAmres <- DNAmres_set[[3]]

##-------------------------------------------------------------------------------------------------------------------------------## 
# Excluding Subjects

# Remove subjects with more than 5% NA
InclSubj <- readRDS("/home/tno010/Desktop/DNAm GE Lung cancer case-control 2016/NOWAC DNAm Therese/Analyser/output/Subjects_to_include_in_analyses.rds")  
# betas
DNAmres <- DNAmres[rownames(DNAmres) %in% InclSubj,]  
dim(DNAmres) # 260 485512 
# Sample sheet
covs <- covs[covs$sampleID %in% InclSubj ,]  
# Align Methylation data to covariates
dim(DNAmres) # 260 485512
DNAmres <- DNAmres[covs$sampleID,]
all(rownames(DNAmres)==covs$sampleID)

##-------------------------------------------------------------------------------------------------------------------------------## 
# Standardizing betas
DNAmres =t(DNAmres) # 485000 rows
dim(DNAmres)
DNAmres_std <- array(NA, c(nrow(DNAmres), ncol(DNAmres)),
                     dimnames=list(rownames(DNAmres),
                                   colnames(DNAmres)))
for (i in 1:nrow(DNAmres)) {
  beta <- DNAmres[i,]
  std_beta <- (beta/sd(beta, na.rm=TRUE))
  DNAmres_std[i,] <- std_beta
  print(i)
}

#saveRDS(DNAmres_std, file="/home/tno010/Desktop/NOWAC DNAm Therese/Analyser/Lung cancer data set/res_NOWAC_lung_MM_Methyl_techcovar_residuals_storagetime_age_std_exclsubj_sampleID.rds")

##-------------------------------------------------------------------------------------------------------------------------------## 
## Unconditional logistic regression

# define model family and corresponding p-value notatio
model.family <-  "binomial"
pval.name <-  "Pr(>|z|)"

# build formula, covariates of interest
covs$case <- covs$case_ctrl=="case"

covars <- c("age.sample", "CSI")
covars <- c("age.sample", "smoking_status_cat", "B", "CD4T", "CD8T", "NK", "Neutrophils", "Eosinophils")
covars <- c("B", "CD4T", "CD8T", "NK", "Neutrophils") #Monocytes out
covars <- c("B", "CD4T", "CD8T", "Monocytes", "Neutrophils") # NK out
cov <- covs[apply(!is.na(covs[,covars]), 1, all),]
dim(cov)

f <- sprintf("%s ~ CpG + %s", "case", paste0(covars, collapse=" + ")) 
f <- sprintf("%s ~ CpG", "case")

#teste på en CpG
DNAmres_std =t(DNAmres_std)  # CpGs as columns
DNAmres_std2 <- DNAmres_std[rownames(DNAmres_std) %in% cov$sampleID,]
DNAmres_std2 <- DNAmres_std2[cov$sampleID,]
all(rownames(DNAmres_std)==covs$sampleID)

CpG <- DNAmres_std2[,9]
testmodel<- glm(as.formula(f), family="binomial", data=cov)
coef(summary(testmodel))

# Build model matrix (without CpG) and extract variable names
X <- model.matrix(as.formula(sprintf("~ %s", paste0(covars, collapse=" + "))), data=cov)
var.names <- c("CpG", colnames(X))
var.names <- c("CpG", "(Intercept)")

#results array
y <- DNAmres_std2#[,c(23423:23443)]

results_logreg <- array(NA, c(length(var.names), ncol(y), 6),
                        dimnames=list(var.names,
                                      colnames(y),
                                      c("nobs", "model.r2", "var.r2",
                                        "coef", "coef.se", "pval")))
for (i in 1:ncol(y)) {
  CpG <- y[,i]
  model <- try(glm(as.formula(f), family=model.family, data=cov),
               silent=TRUE)
  if (!inherits(model, "try-error")) {
    coefs <- try(coef(summary(model)), silent=TRUE)
    if (inherits(coefs, "try-error")) {
      next
    }
    
    # Compute sum of squares for all variables, and RSS
    ss <- effects(model)[setdiff(var.names, "(Intercept)")]**2
    rss <- sum(residuals(model, type="response")**2)
    
    # Standardize sum of squares to unit sum
    ss <- ss / (sum(ss) + rss)
    
    for (v in var.names) if (v %in% rownames(coefs)) {
      results_logreg[v,i,"nobs"] <- nobs(model)
      results_logreg[v,i,"model.r2"] <- sum(ss)
      results_logreg[v,i,"var.r2"] <- ss[v]
      results_logreg[v,i,"coef"] <- coefs[v,"Estimate"]
      results_logreg[v,i,"coef.se"] <- coefs[v,"Std. Error"]
      results_logreg[v,i,"pval"] <- coefs[v, pval.name]
    }
  }
}

saveRDS(results_logreg,file="/home/tno010/Desktop/NOWAC DNAm Therese/Analyser/output/res_logReg_NOWAC_LC_meth_betares_techcov_storagetime_age_std_CpG_WBCadj_notmonoeos_n260_260416.rds")

res.CpG <- as.data.frame(results_logreg["CpG",,])

res.CpG <-res.CpG[which(!is.na(res.CpG$pval)),]
res.CpG$Name<-rownames(res.CpG)
str(res.CpG)# X obs. of  7 variables





##-------------------------------------------------------------------------------------------------------------------------------## 
## Unconditional logistic regression

# Remove unmatched pairs
cases<-covs[covs$case=="case",]
controls<-covs[covs$case=="ctrl",]
cases2 <- cases[!cases$pair.no %in% controls$pair.no,] # 18, 40, 43
controls2 <- controls[!controls$pair.no %in% cases$pair.no,] # 68
exclude <- c(cases2$pair.no, controls2$pair.no)
covs <- covs[!covs$pair.no %in% as.character(exclude),]
rm(cases, cases2, controls, controls2, exclude, InclSubj)

# Align Methylation data to covariates
DNAmres_std <- t(DNAmres_std)
dim(DNAmres_std) # 256 485512
DNAmres2 <- DNAmres_std[rownames(DNAmres_std) %in% covs$sampleID,]
DNAmres2 <- DNAmres2[covs$sampleID,]
all(rownames(DNAmres2)==covs$sampleID)

# test on one CpG
CpG <- DNAmres2[,9]
mod2 <- clogit(covs$case=="case" ~ CpG + strata(pairs), data=covs)
summary(mod2)
mod2$n #nobs
summary(mod2)$rsq[1]  #Rsquared
#var.r2 CpG if first
coef(summary(mod2))[1,1] #coef CpG if first
coef(summary(mod2))[1,3] #coef.se CpG if first
coef(summary(mod2))[1,5] #pval CpG if first

covars <- c("age.sample", "smoking_status_cat", "B", "CD4T", "CD8T", "NK", "Neutrophils", "Eosinophils")
covars <- c("B", "CD4T", "CD8T", "NK", "Neutrophils") #Monocytes out
covars <- c("B", "CD4T", "CD8T", "Monocytes", "Neutrophils") # NK out
cov <- covs[apply(!is.na(covs[,covars]), 1, all),]
dim(cov)

(f <- sprintf("%s ~ CpG + %s", "case", paste0(covars, collapse=" + ")) )
#f <- sprintf("%s ~ CpG", "case")

# Build model matrix (without CpG) and extract variable names
X <- model.matrix(covs$case=="TRUE" ~ CpG + strata(pairs), data=covs)
X <- model.matrix(as.formula(sprintf("~ %s", paste0(covars, collapse=" + "))), data=covs)
var.names <- c("CpG", colnames(X))
var.names <- c("CpG", "(Intercept)")

#results array
y <- DNAmres2#[,c(23423:23430)]

results_condlogreg <- array(NA, c(length(var.names), ncol(y), 5),
                            dimnames=list(var.names,
                                          colnames(y),
                                          c("nobs", "model.r2",
                                            "coef", "coef.se", "pval")))
for (i in 1:ncol(y)) {
  CpG <- y[,i]
  model <- try(clogit(covs$case=="case" ~ CpG + strata(pairs), method="exact", data=cov),
               silent=TRUE)
  if (!inherits(model, "try-error")) {
    coefs <- try(coef(summary(model)), silent=TRUE)
    if (inherits(coefs, "try-error")) {
      next
    }
    
    for (v in var.names) if (v %in% rownames(coefs)) {
      results_condlogreg[v,i,"nobs"] <- model$n
      results_condlogreg[v,i,"model.r2"] <- summary(model)$rsq[1]
      results_condlogreg[v,i,"coef"] <- coefs[1,1]
      results_condlogreg[v,i,"coef.se"] <- coefs[1,3]
      results_condlogreg[v,i,"pval"] <- coefs[1,5]
    }
  }
}

saveRDS(results_condlogreg,file="/home/tno010/Desktop/NOWAC DNAm Therese/Analyser/output/res_cond_logReg_NOWAC_LC_meth_betares_techcov_storagetime_age_std_CpG_WBCadj_notmonoeos_n256_260416.rds")

res.CpG <- as.data.frame(results_condlogreg["CpG",,])

res.CpG <-res.CpG[which(!is.na(res.CpG$pval)),]
res.CpG$Name<-rownames(res.CpG)
str(res.CpG)# X obs. of  7 variables

results_logreg <- readRDS("/home/tno010/Desktop/DNAm GE Lung cancer case-control 2016/NOWAC DNAm Therese/Analyser/output/res_logReg_NOWAC_LC_meth_betares_techcov_storagetime_age_std_CpGonly_n260_210416.rds")
res.CpG <- as.data.frame(results_logreg["CpG",,])
res.CpG$Name<-rownames(res.CpG)
##-------------------------------------------------------------------------------------------------------------------------------## 
## Annotation and top 100 list 

### Load annotationhs
annot<-load("/home/tno010/Desktop/DNAm past information from Imperial BC LC/OMICs Imperial College/Data Analysis involving NOWAC Imperial Florence/Rskripts/annot450K.R")
str(annot)# 485577 12

#### Load Cross-hybridisation probes
add<-read.table("/home/tno010/Desktop/DNAm past information from Imperial BC LC/OMICs Imperial College/R scripts Florence Guida/add_annot_price.txt", sep = "\t", header=TRUE)
add$ID<-as.character(add$ID)
str(add)# 485512 22

# Probes to delete
todel<-add[which(add$XY_Hits!="XY_NO" & add$Autosomal_Hits!="A_NO"), "ID"]
length(todel)# 11101

#### Process RESULTS
# remove chr X and Y
res.CpG2<-merge(res.CpG, annot450K, by ="Name")
dim(res.CpG2)# 465878     18
res.CpG2<-subset(res.CpG2, res.CpG2$CHR!="X" & res.CpG2$CHR != "Y")
dim(res.CpG2)# 454977     18

# Delete cross-hybridised probes
addmargins(table(res.CpG2$Name %in% todel, exclude=NULL))# 10789 T
res.CpG2<- res.CpG2[!(res.CpG2$Name %in% todel),]
dim(res.CpG2)# 444188     18

# Calculate bonferroni threshold
n<-0.05/nrow(res.CpG2)
res.CpG2$bonfthreshold<-n#1.08E-07
res.CpG2$bonf<-ifelse(res.CpG2$pval<res.CpG2$bonfthreshold, "sign", "not_sign")# 52 signif
table(res.CpG2$bonf)

# Calculate fdr threshold
library(fdrtool)
fdr<-fdrtool(res.CpG2[!is.na(res.CpG2$pval),]$pval, statistic = "pvalue")
res.CpG2$qval<-ifelse(!is.na(res.CpG2$pval), fdr$qval, NA)
res.CpG2$fdr<-ifelse(res.CpG2$qval<0.05, "sign", "not_sign")# 305 signif
table(res.CpG2$fdr)

# make top 100 significant CpG
res.CpG2<-res.CpG2[order(res.CpG2$pval),]
top_100<-res.CpG2[1:10000,]
head(top_100)

write.csv(top_100, file="/home/tno010/Desktop/NOWAC DNAm Therese/Analyser/output/Top10000CpG_logReg_NOWAC_LC_meth_betares_techcov_storagetime_age_std_CpGonly_n260_210416_250516.csv")

signCpGs <- top_100$Name[top_100$bonf=="sign"]
saveRDS(signCpGs, file="/home/tno010/Desktop/NOWAC DNAm Therese/Analyser/output/Significant_CpGs_betares_storagetime_age_std_WBCadj_notNKeos220416.rds")

res.CpG2$pos<-substring(res.CpG2$UCSC_CpG_Islands_Name, 7,14)
PCs <-res.CpG2[which(res.CpG2$UCSC_RefGene_Name=="PC"),]
plot(PCs$pval)
plot(density(PCs$pval, na.rm=T),col='blue', xlab='p-value', main="Distribution p-values for CpGs on the PC gene",
     ylab='Density', xlim=c(-0.5,1.2), ylim=c(0,1.7))

plot(PCs$pval, PCs$pos, xlim=c(0,0.2))
write.csv(top_100, file="/home/tno010/Desktop/NOWAC DNAm Therese/Analyser/output/Top10000CpG_logReg_NOWAC_LC_meth_betares_techcov_storagetime_age_std_CpGonly_n260_210416_250516.csv")
