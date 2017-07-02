##-------------------------------------------------------------------------------------------------------------------------------## 
#                Creating smoking variables in NOWAC
##-------------------------------------------------------------------------------------------------------------------------------## 
#  Smoking status decided by Florence and Therese, see separate document for reasoning. BSMOKE
#  Intensity: Based on SAS code from Nicolle and Tonje with intervals trying to harminoze data from different questionnaires
rm(list=ls())

covs <- readRDS("/home/tno010/Desktop/NOWAC DNAm Therese/Analyser/output/covs_8jul+bloodq+update.rds")  # 263 women, 199 vars

##-------------------------------------------------------------------------------------------------------------------------------## 
## Smoking status variables

# Creating three CFN categories ourselves
covs$smoking.status <- NA
covs$smoking.status[covs$BROYK==0] <- "Current" 
covs$smoking.status[covs$smoking.status=="Current" & covs$SIGROYK==1 & !is.na(covs$SIGROYK) & !is.na(covs$smoking.status)]<-"Never" 
covs$smoking.status[covs$BROYK==1 & (covs$EVERROK==1 | covs$SIGROYK==1)]<-"Never" 
covs$smoking.status[is.na(covs$smoking.status) & (covs$EVERROK==0 | covs$SIGROYK==0)]<-"Former" 
covs$smoking.status[is.na(covs$BROYK)]<-"Never"
covs$smoking.status[covs$smoking.status=="Former" & covs$BROYK==1 & covs$yEVERROK==1 &!is.na(covs$smoking.status) & covs$ZSIGROYK==1 &!is.na(covs$ZSIGROYK)]<-"Never"
covs$smoking_status_cat <- factor(covs$smoking.status, labels =c("Current", "Former", "Never"))
covs$smoking_status_cat <- factor(covs$smoking.status, levels = c("Never", "Former", "Current"))

### Create a variable for people with inconsistent information "liars"
covs$problem[covs$smoking.status=="Current" & covs$SIGROYK==1 & !is.na(covs$SIGROYK) & !is.na(covs$smoking.status)]<-1
covs$problem[covs$smoking.status=="Current" & !is.na(covs$smoking.status) 
             & ((!is.na(covs$yROYKSTOP) & covs$yROYKSTOP<=covs$age.sample) | (!is.na(covs$ZROYKSTOP) & covs$ZROYKSTOP<=covs$age.sample))]<-2
covs$problem[covs$smoking.status=="Never" & covs$BROYK==1 &!is.na(covs$smoking.status) & covs$ZSIGROYK==0 &!is.na(covs$ZSIGROYK)]<-1
covs$problem[covs$smoking.status=="Former" & covs$BROYK==1 &!is.na(covs$smoking.status) & covs$ZSIGROYK==1 &!is.na(covs$ZSIGROYK)]<-1
#covs$smoking.status[is.na(covs$BROYK)]<-"Never" #conditional considering breast set

#saveRDS(covs[,c("labnr", "smoking.status", "problem")], file="/home/tno010/Desktop/NOWAC DNAm Therese/Analyser/output/Lung_smoking.status.rds")
t[,c("labnr","BROYK","BROYKANTGAR", "age.sample","EVERROK","yEVERROK", "ROYKNAA", "yROYKNAA", "ZROYKNAA", "SIGROYK","ySIGROYK", "ZSIGROYK","SIGALDER", "YSIGALDER", "ZSIGALDER", "roykald", "yROYKSTOP", "ZROYKSTOP", "ZROKSIST5", "ZROYKAVTIL","packyear", "roykaar")]

# Comparing manual categories
covs$smoking_status_cat<-factor(covs$smoking.status, labels =c("Current", "Former", "Never"))
table(covs$Smoking_status_T, covs$smoking.status)
addmargins(table(covs$Francesca_classification, covs$smoking.status))

##      Ever - never smokers
covs$evernever <- ifelse(covs$smoking_status_cat=="Never", "Never", "Ever");covs$evernever <-factor(covs$evernever, levels=c("Never", "Ever"))

## Experiencing passive smoking as children
covs$passiv_child <- ifelse(!is.na(covs$ROKBARN), as.numeric(covs$ROKBARN), NA)
covs$passiv_child <- ifelse(!is.na(covs$yROKBARN), as.numeric(covs$yROKBARN), covs$passiv_child)
covs$passiv_child <- ifelse(!is.na(covs$ZBARNROYK), as.numeric(covs$ZBARNROYK), covs$passiv_child)
covs$passiv_child <- ifelse((is.na(covs$passiv_child) & covs$smoking_status_cat=="Never"), 0, covs$passiv_child)  # 2 women
table(covs$smoking_status_cat, covs$passiv_child, exclude=NULL)

## Experiencing passive smoking as adults
covs$passiv_adult <- ifelse(!is.na(covs$ROKBOR), as.numeric(covs$ROKBARN), NA)
covs$passiv_adult <- ifelse(!is.na(covs$yROKBOR), as.numeric(covs$yROKBOR), covs$passiv_adult)
covs$passiv_adult <- ifelse(is.na(covs$passiv_adult), 1, covs$passiv_adult)
covs$passiv_adult <- ifelse(covs$passiv_adult==3, 0, covs$passiv_adult)

covs$passiv_ad_int <- ifelse(!is.na(covs$ROKBORNO), as.numeric(covs$ROKBORNO), NA)
covs$passiv_ad_int <- ifelse(!is.na(covs$yROKBORNO), as.numeric(covs$yROKBORNO), covs$passiv_ad_int)
covs$passiv_ad_int <- ifelse(covs$passiv_ad_int==98, NA, covs$passiv_ad_int)  # 98 is not sure, put as NA
covs$passiv_ad_int <- ifelse((covs$passiv_adult==1 & is.na(covs$passiv_ad_int)), 0, covs$passiv_ad_int) 
addmargins(table(covs$passiv_adult, covs$passiv_ad_int, exclude=NULL)) # 19 NAs that have smokers at home, but without intensity

covs$passiv_ad_int_cat <- as.factor(cut(covs$passiv_ad_int, breaks=c(0, 15, 30, 60)))
covs$passiv_ad_int_cat <- ifelse(is.na(covs$passiv_ad_int_cat), "no", covs$passiv_ad_int_cat)
covs$passiv_ad_int_cat <- factor(covs$passiv_ad_int_cat, labels = c("low", "medium", "high", "none"))
table(covs$passiv_ad_int_cat, exclude=NULL)

## go down to where the finished ROK is loaded as RDS

## Harmonizing intensity questions in ROK matrix
ROK <- matrix(NA, nrow=nrow(covs), ncol=12, dimnames=list(rownames(covs), c("ROK1", "ROK2", "ROK3", "ROK4", "ROK5", "ROK6", "ROK7", "ROK8", "ROK9", "ROK10", "ROK11", "ROK12")))

serie0<-c(1:10)
ROK[,1] <- ifelse(covs$serienr %in% serie0 & !is.na(covs$ROYKANT1014), covs$ROYKANT1014, ROK[,1])
ROK[,2] <- ifelse(covs$serienr %in% serie0 & !is.na(covs$ROYKANT1519), covs$ROYKANT1519, ROK[,2])
ROK[,3] <- ifelse(covs$serienr %in% serie0 & !is.na(covs$roykant2024), covs$roykant2024, ROK[,3]) 
ROK[,4] <- ifelse(covs$serienr %in% serie0 & !is.na(covs$roykant2529), covs$roykant2529, ROK[,4]) 
ROK[,5] <- ifelse(covs$serienr %in% serie0 & !is.na(covs$roykant3034), covs$roykant3034, ROK[,5]) 
ROK[,6] <- ifelse(covs$serienr %in% serie0 & !is.na(covs$roykant3539), covs$roykant3539, ROK[,6]) 
ROK[,7] <- ifelse(covs$serienr %in% serie0 & !is.na(covs$roykant4044), covs$roykant4044, ROK[,7])
ROK[,8] <- ifelse(covs$serienr %in% serie0 & !is.na(covs$roykant4549), covs$roykant4549, ROK[,8])
ROK[,9] <- ifelse(covs$serienr %in% serie0, NA, ROK[,9])
ROK[,10]  <- ifelse(covs$serienr %in% serie0, NA, ROK[,10])
ROK[,11] <- ifelse(covs$serienr %in% serie0, NA, ROK[,11])
ROK[,12] <- ifelse(covs$serienr %in% serie0, NA, ROK[,12])

serie1<-c(11,12,13,14,15,16,19,20,21,22,23,24)
ROK[,1] <- ifelse(covs$serienr %in% serie1 & !is.na(covs$roykant1019), covs$roykant1019, ROK[,1])
ROK[,2] <- ifelse(covs$serienr %in% serie1 & !is.na(covs$roykant1019), covs$roykant1019, ROK[,2])
ROK[,3] <- ifelse(covs$serienr %in% serie1 & !is.na(covs$ROYKANT2029), covs$ROYKANT2029, ROK[,3]) 
ROK[,4] <- ifelse(covs$serienr %in% serie1 & !is.na(covs$ROYKANT2029), covs$ROYKANT2029, ROK[,4]) 
ROK[,5] <- ifelse(covs$serienr %in% serie1 & !is.na(covs$ROYKANT3039), covs$ROYKANT3039, ROK[,5]) 
ROK[,6] <- ifelse(covs$serienr %in% serie1 & !is.na(covs$ROYKANT3039), covs$ROYKANT3039, ROK[,6]) 
ROK[,7] <- ifelse(covs$serienr %in% serie1 & !is.na(covs$ROYKANT4049), covs$ROYKANT4049, ROK[,7])
ROK[,8] <- ifelse(covs$serienr %in% serie1 & !is.na(covs$ROYKANT4049), covs$ROYKANT4049, ROK[,8])
ROK[,9] <- ifelse(covs$serienr %in% serie1 & !is.na(covs$roykant5059), covs$roykant5059, ROK[,9])
ROK[,10] <- ifelse(covs$serienr %in% serie1 & !is.na(covs$roykant5059), covs$roykant5059, ROK[,10])
ROK[,11] <- ifelse(covs$serienr %in% serie1, NA, ROK[,11])
ROK[,12] <- ifelse(covs$serienr %in% serie1, NA, ROK[,12])

serie2 <- c(35:45)
ROK[,1] <- ifelse(covs$serienr %in% serie2 | covs$zserienr %in% serie2 & !is.na(covs$ROYKANT1014), covs$ROYKANT1014, ROK[,1])
ROK[,2] <- ifelse(covs$serienr %in% serie2 | covs$zserienr %in% serie2 & !is.na(covs$ROYKANT1519), covs$ROYKANT1519, ROK[,2])
ROK[,3] <- ifelse(covs$serienr %in% serie2 | covs$zserienr %in% serie2 & !is.na(covs$ROYKANT2029), covs$ROYKANT2029, ROK[,3])
ROK[,4] <- ifelse(covs$serienr %in% serie2 | covs$zserienr %in% serie2 & !is.na(covs$ROYKANT2029), covs$ROYKANT2029, ROK[,4])
ROK[,5] <- ifelse(covs$serienr %in% serie2 | covs$zserienr %in% serie2 & !is.na(covs$ROYKANT3039), covs$ROYKANT3039, ROK[,5])
ROK[,6] <- ifelse(covs$serienr %in% serie2 | covs$zserienr %in% serie2 & !is.na(covs$ROYKANT3039), covs$ROYKANT3039, ROK[,6])
ROK[,7] <- ifelse(covs$serienr %in% serie2 | covs$zserienr %in% serie2 & !is.na(covs$ROYKANT4049), covs$ROYKANT4049, ROK[,7])
ROK[,8] <- ifelse(covs$serienr %in% serie2 | covs$zserienr %in% serie2 & !is.na(covs$ROYKANT4049), covs$ROYKANT4049, ROK[,8])
ROK[,9] <- ifelse(covs$serienr %in% serie2 | covs$zserienr %in% serie2 & !is.na(covs$ROYKANT50MM), covs$ROYKANT50MM, ROK[,9])
ROK[,10] <- ifelse(covs$serienr %in% serie2 | covs$zserienr %in% serie2 & !is.na(covs$ROYKANT50MM), covs$ROYKANT50MM, ROK[,10])
ROK[,11] <- ifelse(covs$serienr %in% serie2 | covs$zserienr %in% serie2, NA, ROK[,11])
ROK[,12] <- ifelse(covs$serienr %in% serie2 | covs$zserienr %in% serie2, NA, ROK[,12])

serie3 <- 26
ROK[,1] <- ifelse(covs$yserienr %in% serie3 & !is.na(covs$ROYKANT1014), NA, ROK[,1])
ROK[,2] <- ifelse(covs$yserienr %in% serie3 & !is.na(covs$ROYKANT1519), covs$ROYKANT1519, ROK[,2])
ROK[,3] <- ifelse(covs$yserienr %in% serie3 & !is.na(covs$roykant2024), covs$roykant2024, ROK[,3]) 
ROK[,4] <- ifelse(covs$yserienr %in% serie3 & !is.na(covs$roykant2529), covs$roykant2529, ROK[,4]) 
ROK[,5] <- ifelse(covs$yserienr %in% serie3 & !is.na(covs$roykant3034), covs$roykant3034, ROK[,5]) 
ROK[,6] <- ifelse(covs$yserienr %in% serie3 & !is.na(covs$roykant3539), covs$roykant3539, ROK[,6]) 
ROK[,7] <- ifelse(covs$yserienr %in% serie3 & !is.na(covs$roykant4044), covs$roykant4044, ROK[,7])
ROK[,8] <- ifelse(covs$yserienr %in% serie3 & !is.na(covs$roykant4549), covs$roykant4549, ROK[,8])
ROK[,9] <- ifelse(covs$yserienr %in% serie3, NA, ROK[,9])
ROK[,10]  <- ifelse(covs$yserienr %in% serie3, NA, ROK[,10])
ROK[,11] <- ifelse(covs$yserienr %in% serie3, NA, ROK[,11])
ROK[,12] <- ifelse(covs$yserienr %in% serie3, NA, ROK[,12])

serie4 <- c(25,27:29)
age94 <- (1994-as.numeric(format(covs$birthyr, "%Y")))
age98 <- (1998-as.numeric(format(covs$birthyr, "%Y")))
ROK[,5] <- ifelse(covs$yserienr %in% serie4 & !is.na(covs$yROYKAR1) & (31 <= age94 & age94 <=35), covs$yROYKAR1, ROK[,5])
ROK[,6] <- ifelse(covs$yserienr %in% serie4 & !is.na(covs$yROYKAR1) & (36<= age94 & age94 <=40), covs$yROYKAR1, ROK[,6])
ROK[,7] <- ifelse(covs$yserienr %in% serie4 & !is.na(covs$yROYKAR1) & (41<= age94 & age94 <=45), covs$yROYKAR1, ROK[,7])
ROK[,8] <- ifelse(covs$yserienr %in% serie4 & !is.na(covs$yROYKAR1) & (46<= age94 & age94 <=50), covs$yROYKAR1, ROK[,8])
ROK[,9] <- ifelse(covs$yserienr %in% serie4 & !is.na(covs$yROYKAR1) & (51<= age94 & age94 <=55), covs$yROYKAR1, ROK[,9])
ROK[,10] <- ifelse(covs$yserienr %in% serie4 & !is.na(covs$yROYKAR1) & (56<= age94 & age94 <=60), covs$yROYKAR1, ROK[,10])
ROK[,5] <- ifelse(covs$yserienr %in% serie4 & !is.na(covs$yROYKAR2) & (31<= age98 & age98 <=35), covs$yROYKAR2, ROK[,5])
ROK[,6] <- ifelse(covs$yserienr %in% serie4 & !is.na(covs$yROYKAR2) & (36<= age98 & age98 <=40), covs$yROYKAR2, ROK[,6])
ROK[,7] <- ifelse(covs$yserienr %in% serie4 & !is.na(covs$yROYKAR2) & (41<= age98 & age98 <=45), covs$yROYKAR2, ROK[,7])
ROK[,8] <- ifelse(covs$yserienr %in% serie4 & !is.na(covs$yROYKAR2) & (46<= age98 & age98 <=50), covs$yROYKAR2, ROK[,8])
ROK[,9] <- ifelse(covs$yserienr %in% serie4 & !is.na(covs$yROYKAR2) & (51<= age98 & age98 <=55), covs$yROYKAR2, ROK[,9])
ROK[,10] <- ifelse(covs$yserienr %in% serie4 & !is.na(covs$yROYKAR2) & (56<= age98 & age98 <=60), covs$yROYKAR2, ROK[,10])

serie5 <- c(32,33)
age01 <- (2001-as.numeric(format(covs$birthyr, "%Y")))
ROK[,7] <- ifelse(covs$yserienr %in% serie5 & !is.na(covs$yROYKAR1) & (41<= age01 & age01 <=45), covs$yROYKAR1, ROK[,7])
ROK[,8] <- ifelse(covs$yserienr %in% serie5 & !is.na(covs$yROYKAR1) & (46<= age01 & age01 <=50), covs$yROYKAR1, ROK[,8])
ROK[,9] <- ifelse(covs$yserienr %in% serie5 & !is.na(covs$yROYKAR1) & (51<= age01 & age01 <=55), covs$yROYKAR1, ROK[,9])
ROK[,10] <- ifelse(covs$yserienr %in% serie5 & !is.na(covs$yROYKAR1) & (56<= age01 & age01 <=60), covs$yROYKAR1, ROK[,10])

serie6 <- c(39,42)
ROK[,7] <- ifelse(covs$zserienr %in% serie6 & !is.na(covs$ZROKSIST5) & (42<= covs$age.quez & covs$age.quez <=46), covs$ZROKSIST5, ROK[,7])
ROK[,8] <- ifelse(covs$zserienr %in% serie6 & !is.na(covs$ZROKSIST5) & (47<= covs$age.quez & covs$age.quez <=51), covs$ZROKSIST5, ROK[,8])
ROK[,9] <- ifelse(covs$zserienr %in% serie6 & !is.na(covs$ZROKSIST5) & (52<= covs$age.quez & covs$age.quez <=56), covs$ZROKSIST5, ROK[,9])
ROK[,10] <- ifelse(covs$zserienr %in% serie6 & !is.na(covs$ZROKSIST5) & (57<= covs$age.quez & covs$age.quez <=61), covs$ZROKSIST5, ROK[,10])

rm(age01, age94, age98, serie0, serie1, serie2, serie3, serie4, serie5, serie6)

ROK <- as.data.frame(ROK)
ROK$IntNA<-rep(0,nrow(ROK))
for (i in 1:nrow(ROK[,c(1:12)]))
{
  ROK$IntNA[i]<-length(which(is.na(ROK[i,c(1:12)]) | ROK[i,c(1:12)]==0))
}
Int.allNAs <- which(ROK$IntNA==12)  # 75 women

ROK$Intervals <- 12- ROK$IntNA
ROK$smokeyears <- ROK$Intervals*5
table(ROK$smokeyears, exclude=NULL)

ROK$Startald.NA<-rep(0,nrow(ROK))
for (i in 1:nrow(ROK[,c(1:12)]))
{
  ROK$Startald.NA[i] <- min(which(ROK[i,c(1:12)] > 0))
}

ROK$Stopald.NA<-rep(0,nrow(ROK))
for (i in 1:nrow(ROK[,c(1:12)]))
{
  ROK$Stopald.NA[i] <- max(which(ROK[i,c(1:12)] > 0))
}
ROK[which(ROK$Startald.NA==Inf),]$Startald.NA <- NA 
ROK[which(ROK$Stopald.NA==-Inf),]$Stopald.NA <- NA
ROK$Startald.NA[ROK$smoking_status_cat=="Never"] <- NA

#ROK[is.na(ROK)] <- 0
ROK[Int.allNAs,] <- NA

map.intensity <- function(x, include.zero) {
  mapping <- c(0, 2.5, 7, 12, 17, 22, 27) # number of sigarettes representing interval in questionnaire
  mapping[as.factor(x)]
}

ROK$int.ROK1 <- map.intensity(ROK$ROK1, include.zero=TRUE)
ROK$int.ROK2 <- map.intensity(ROK$ROK2, include.zero=TRUE)
ROK$int.ROK3 <- map.intensity(ROK$ROK3, include.zero=TRUE)
ROK$int.ROK4 <- map.intensity(ROK$ROK4, include.zero=TRUE)
ROK$int.ROK5 <- map.intensity(ROK$ROK5, include.zero=TRUE)
ROK$int.ROK6 <- map.intensity(ROK$ROK6, include.zero=TRUE)
ROK$int.ROK7 <- map.intensity(ROK$ROK7, include.zero=TRUE)
ROK$int.ROK8 <- map.intensity(ROK$ROK8, include.zero=TRUE)
ROK$int.ROK9 <- map.intensity(ROK$ROK9, include.zero=TRUE)
ROK$int.ROK10 <- map.intensity(ROK$ROK10, include.zero=TRUE)
ROK$int.ROK11 <- map.intensity(ROK$ROK11, include.zero=TRUE)
ROK$int.ROK12 <- map.intensity(ROK$ROK12, include.zero=TRUE)

ROK$int.ROK1 <- ifelse(ROK$int.ROK1 > 0, ROK$int.ROK1, NA)
ROK$int.ROK2 <- ifelse(ROK$int.ROK2 > 0, ROK$int.ROK2, NA)
ROK$int.ROK3 <- ifelse(ROK$int.ROK3 > 0, ROK$int.ROK3, NA)
ROK$int.ROK4 <- ifelse(ROK$int.ROK4 > 0, ROK$int.ROK4, NA)
ROK$int.ROK5 <- ifelse(ROK$int.ROK5 > 0, ROK$int.ROK5, NA)
ROK$int.ROK6 <- ifelse(ROK$int.ROK6 > 0, ROK$int.ROK6, NA)
ROK$int.ROK7 <- ifelse(ROK$int.ROK7 > 0, ROK$int.ROK7, NA)
ROK$int.ROK8 <- ifelse(ROK$int.ROK8 > 0, ROK$int.ROK8, NA)
ROK$int.ROK9 <- ifelse(ROK$int.ROK9 > 0, ROK$int.ROK9, NA)
ROK$int.ROK10 <- ifelse(ROK$int.ROK10 > 0, ROK$int.ROK10, NA)

ROK$mean.int <- rowMeans(ROK[,c("int.ROK1", "int.ROK2", "int.ROK3", "int.ROK4", "int.ROK5", "int.ROK6",
                                "int.ROK7", "int.ROK8", "int.ROK9", "int.ROK10")], na.rm = TRUE) #, "int.ROK11", "int.ROK12"
summary(ROK$int.ROK3)
summary(ROK$mean.int)

#saveRDS(ROK, file="/home/tno010/Desktop/NOWAC DNAm Therese/Analyser/output/Lung_smoke_intensity_ROK_matrix_030316.rds")

##

ROK <- readRDS("/home/tno010/Desktop/NOWAC DNAm Therese/Analyser/output/Lung_smoke_intensity_ROK_matrix_030316.rds")  
covs<-merge(covs, ROK, by.x="row.names", by.y="row.names") 

covs[covs$smoking_status_cat=="Former",c("int.ROK1", "int.ROK2", "int.ROK3", "int.ROK4", "int.ROK5", "int.ROK6", "int.ROK7", "int.ROK8", "int.ROK9", "int.ROK10", "int.ROK11", "int.ROK12")]
covs[covs$smoking_status_cat=="Current",c("int.ROK1", "int.ROK2", "int.ROK3", "int.ROK4", "int.ROK5", "int.ROK6", "int.ROK7", "int.ROK8", "int.ROK9", "int.ROK10", "int.ROK11", "int.ROK12")]

t <- covs[which(is.na(covs$mean.int)),]  # 75 women
t <- covs[which(is.na(covs$mean.int) & covs$BROYK==1),] # 69 women
t[,c("labnr", "smoking.status","age.sample", "age.que1", "age.quey", "age.quez", "ROYKNAA", "yROYKNAA", "ZROYKNAA", "age.start", "roykald","SIGALDER", "YSIGALDER", "ZSIGALDER", "roykant1019", "ROYKANT2029","ROYKANT1014", "ROYKANT1519")]
t[,c("labnr", "smoking.status","age.sample", "age.start", "roykald", "age.stop", "stopald", "yROYKSTOP", "ZROYKSTOP", "Intervals", "Startald.NA", "Stopald.NA", "roykaar", "duration.max", "duration.adj")]

covs[which(is.na(covs$mean.int) & covs$BROYK==1 & (covs$ROYKNAA==0 | covs$yROYKNAA==0 | covs$ZROYKNAA==0)),]  # no one
hist(round(covs$mean.int), main="Mean intensity cigarette smoking", xlab="cigarettes per day")

covs[which((is.na(covs$mean.int) | covs$mean.int>0) & covs$smoking_status_cat=="Never"),]$mean.int <- 0  # no one seem to not have smoked considerably, 2 Never women have intensities of 1
covs[which(is.na(covs$mean.int) & covs$BROYK==0 & (covs$ROYKNAA==0 | covs$yROYKNAA==0 | covs$ZROYKNAA==0)),]$mean.int <- NA  # 2 Current women
covs[which(is.nan(covs$mean.int)),]$mean.int <- NA # 2 current, 2 formers missing info intensity
#covs[which(covs$mean.int==0 & covs$smoking_status_cat=="Former"),]$mean.int  <- NA #no one?

summary(covs$mean.int)
addmargins(table(round(covs[covs$smoking_status_cat=="Never",]$mean.int), exclude=NULL))
addmargins(table(round(covs[covs$smoking_status_cat=="Former",]$mean.int), exclude=NULL))
addmargins(table(round(covs[covs$smoking_status_cat=="Current",]$mean.int), exclude=NULL))

##          Smoke start age
addmargins(table(covs$smoking_status_cat, covs$Startald.NA, exclude=NULL))
t <- covs[which(is.na(covs$Startald.NA) & !covs$smoking_status_cat=="Never"),]  # 4 current women, 2 former women

covs$smokeald <- with(covs, ifelse(Startald.NA==1, 12.5, ifelse(Startald.NA==2, 17.5, ifelse(Startald.NA==3, 22.5, ifelse(Startald.NA==4, 27.5, ifelse(Startald.NA==5, 32.5, ifelse(Startald.NA==6, 37.5, 
                                                                                                                                                                                    ifelse(Startald.NA==7, 42.5, ifelse(Startald.NA==8, 47.5, ifelse(Startald.NA==9, 52.5, ifelse(Startald.NA==10, 57.5, ifelse(Startald.NA==11, 62.5, ifelse(Startald.NA==12, 67.5, NA)))))))))))))

summary(covs$smokeald)
addmargins(table(covs$smoking_status_cat, covs$smokeald, exclude=NULL))

##      Age at start
covs$age.start<-ifelse((!is.na(covs$SIGALDER) ), covs$SIGALDER, covs$smokeald) # & covs$SIGALDER >= 15
covs$age.start<-ifelse((!is.na(covs$YSIGALDER) ), covs$YSIGALDER, covs$age.start) # & covs$YSIGALDER >= 15
covs$age.start<-ifelse((!is.na(covs$ZSIGALDER) ), covs$ZSIGALDER, covs$age.start) # & covs$ZSIGALDER >= 15
covs$age.start<-ifelse(covs$smoking_status_cat=="Never", covs$age.sample, covs$age.start)
#covs$age.start<-ifelse(covs$age.start<15, 15, covs$age.start)

addmargins(table(covs[!covs$smoking_status_cat=="Never",]$age.start, covs[!covs$smoking_status_cat=="Never",]$roykald, exclude=NULL))
addmargins(table(covs$smoking_status_cat, covs$age.start, exclude=NULL))
addmargins(table(covs$smoking_status_cat, covs$roykald, exclude=NULL))
addmargins(table(covs$roykald, covs$age.start, exclude=NULL))
hist(covs$age.start, main="Age at start", xlab="years")
plot(covs$age.start, covs$roykald, col=as.factor(covs$smoking_status_cat), pch=16); abline(a=0, b=1)

t <- covs[which(!covs$age.start==covs$roykald & !covs$smoking_status_cat=="Never"),] # 50 stk
t <- covs[which(covs$age.start>=covs$age.stop),] # 2 stk

##      Age at stop
covs$stopald <- with(covs, ifelse(Stopald.NA==1, 12.5, ifelse(Stopald.NA==2, 17.5, ifelse(Stopald.NA==3, 22.5, ifelse(Stopald.NA==4, 27.5, ifelse(Stopald.NA==5, 32.5, ifelse(Stopald.NA==6, 37.5, 
                                                                                                                                                                              ifelse(Stopald.NA==7, 42.5, ifelse(Stopald.NA==8, 47.5, ifelse(Stopald.NA==9, 52.5, ifelse(Stopald.NA==10, 57.5, ifelse(Stopald.NA==11, 62.5, ifelse(Stopald.NA==12, 67.5, NA)))))))))))))
covs$age.stop<-ifelse(covs$smoking_status_cat=="Former" & !is.na(covs$yROYKSTOP) & (covs$age.sample>=covs$yROYKSTOP) & ((covs$yROYKSTOP+3)>(covs$stopald)),
                      covs$yROYKSTOP, covs$stopald)
covs$age.stop<-ifelse(covs$smoking_status_cat=="Former" & !is.na(covs$ZROYKSTOP) & (covs$age.sample>=covs$ZROYKSTOP) & ((covs$ZROYKSTOP+3)>(covs$stopald)),
                      covs$ZROYKSTOP, covs$age.stop)
covs$age.stop<-ifelse((!covs$smoking_status_cat=="Former"), NA, covs$age.stop)

addmargins(table(covs$stopald, covs$age.stop, exclude=NULL))
addmargins(table(covs$smoking_status_cat, covs$age.stop, exclude=NULL))

#covs$age.stop2 <- (covs$age.start + covs$Intervals*5)
addmargins(table(covs$age.stop2, covs$stopald, exclude=NULL))
hist(covs$age.stop, main="Age at stop", xlab="years")

t <- covs[which(covs$Stopald.NA==covs$Startald.NA),] # 6 stk

#covs$roykstop <- ifelse(covs$smoking_status_cat=="Former", (covs$roykald + covs$roykaar), NA)
covs$duration.max <- ifelse(covs$smoking_status_cat=="Never", 0, NA)
covs$duration.max <- ifelse(covs$smoking_status_cat=="Former", covs$age.stop - covs$age.start, covs$duration.max)
covs$duration.max <- ifelse(covs$smoking_status_cat=="Current", covs$age.sample - covs$age.start, covs$duration.max)
#covs[which(covs$duration.max<0),]$duration.max <- NA

addmargins(table(covs$smoking_status_cat, covs$roykstop, exclude=NULL))
addmargins(table(covs$smoking_status_cat, covs$duration.max, exclude=NULL))
addmargins(table(round(covs[covs$smoking_status_cat=="Never",]$duration.max), exclude=NULL))
addmargins(table(round(covs[covs$smoking_status_cat=="Former",]$duration.max), exclude=NULL))
addmargins(table(round(covs[covs$smoking_status_cat=="Current",]$duration.max), exclude=NULL))

plot(covs$age.stop, covs$stopald, col=as.factor(covs$smoking_status_cat), pch=16); abline(a=0, b=1)
plot(covs$duration.max, covs$roykaar, col=as.factor(covs$smoking_status_cat), pch=16); abline(a=0, b=1)

covs$Intervals.actu <- covs$Intervals
covs$Intervals.teor <- (covs$Stopald.NA - covs$Startald.NA)+1

table(covs$Intervals.actu, covs$Intervals.teor)
t <- covs[which(!covs$age.stop==covs$roykstop & !covs$smoking_status_cat=="Never"),] # 66 stk
t <- covs[which(covs$smoking_status_cat=="Former" & covs$duration.max<10),] # 9 stk
t <- covs[which(!covs$Intervals.actu==covs$Intervals.teor & covs$smoking_status_cat=="Former"),] # 17 stk
t <- covs[which(!covs$Intervals.actu==covs$Intervals.teor & covs$smoking_status_cat=="Current"),] # 17 stk

covs$diff.int <- ifelse(!(covs$Startald.NA==covs$Stopald.NA),  (covs$Intervals.teor - covs$Intervals.actu), -0.5)
covs$duration.adj <- ifelse((!covs$smoking_status_cat=="Never" & !is.na(covs$diff.int)), (covs$duration.max - (covs$diff.int*5)), covs$duration.max)
covs$duration.adj <- ifelse((covs$diff.int<0 & covs$duration.max>1), covs$duration.max, covs$duration.adj)
covs$duration.adj <- ifelse(covs$smoking_status_cat=="Never", 0, covs$duration.adj)

covs$duration.int <- (covs$Intervals.actu*5)
plot(covs$duration.int, covs$duration.max, pch=16)
addmargins(table(round(covs[covs$smoking_status_cat=="Never",]$duration.adj), exclude=NULL))
addmargins(table(round(covs[covs$smoking_status_cat=="Former",]$duration.adj), exclude=NULL))
addmargins(table(round(covs[covs$smoking_status_cat=="Current",]$duration.adj), exclude=NULL))

t <- covs[which(!covs$smoking_status_cat=="Never" & covs$age.stop>covs$stopald+3),]
t <- covs[which(covs$smoking_status_cat=="Current" & covs$duration.max>30 & covs$roykaar<25),] # 8 stk
t <- covs[which(covs$smoking_status_cat=="Former" & covs$duration.max>30 & covs$ZROYKSTOP<25),] # 8 stk

plot(covs$duration.max, covs$duration.adj, col=as.factor(covs$smoking_status_cat), pch=16); abline(a=0, b=1)
plot(covs$duration.adj, covs$duration.int, col=as.factor(covs$smoking_status_cat), pch=16); abline(a=0, b=1)
plot(covs$roykaar, covs$duration.adj, col=as.factor(covs$smoking_status_cat), pch=16); abline(a=0, b=1)

hist(covs$duration.max, main=paste("Duration of smoking", " \n ", "- Disregarding breaks"), xlab="years")
hist(covs$duration.adj, main=paste("Duration of smoking", " \n ", "- Including breaks"), xlab="years")

##      Time since quitting for formers
covs$tsq_smoking<-ifelse(covs$smoking_status_cat=="Former" & !is.na(covs$age.stop) & (covs$age.sample>=covs$age.stop),
                         covs$age.sample-covs$age.stop, NA)
covs$tsq_smoking<-ifelse(covs$smoking_status_cat=="Never", covs$age.sample, covs$tsq_smoking)
covs$tsq_smoking<-ifelse(covs$smoking_status_cat=="Current", 0, covs$tsq_smoking)

addmargins(table(covs$smoking_status_cat, covs$tsq_smoking, exclude=NULL))
t <- covs[covs$tsq_smoking<2 & covs$smoking_status_cat=="Former",]
hist(covs$tsq_smoking, main="Time since quitting smoking", xlab="years")
t <- covs[which(covs$tsq_smoking<20 & covs$duration.max<10),] # 10 stk

##      Fasanellis two categories for formers, over or under 10 years TSQ
covs$fasanelli_formersmoke_twocat <- ifelse(covs$smoking_status_cat=="Former" & covs$tsq_smoking < 10, "Former <10 years", NA)
covs$fasanelli_formersmoke_twocat <- ifelse(covs$smoking_status_cat=="Former" & covs$tsq_smoking >= 10, "Former >10 years", covs$fasanelli_formersmoke_twocat)
covs$fasanelli_formersmoke_twocat <- ifelse(covs$smoking_status_cat=="Current", "Current", covs$fasanelli_formersmoke_twocat)
covs$fasanelli_formersmoke_twocat <- ifelse(covs$smoking_status_cat=="Never", "Never", covs$fasanelli_formersmoke_twocat)
covs$fasanelli_formersmoke_twocat <- factor(covs$fasanelli_formersmoke_twocat, levels = c("Current", "Former <10 years", "Former >10 years", "Never"))

##      Intensity metrics
## Cumulative mean intensity

addmargins(table(covs$diff, exclude=NULL))
covs$cum.mean.int <- (rowSums(covs[,c("int.ROK1", "int.ROK2", "int.ROK3", "int.ROK4", "int.ROK5", "int.ROK6",
                                      "int.ROK7", "int.ROK8", "int.ROK9", "int.ROK10")], na.rm = TRUE) *5) / (covs$Intervals.actu*5)
covs$cum.mean.int2 <- (rowSums(covs[,c("int.ROK1", "int.ROK2", "int.ROK3", "int.ROK4", "int.ROK5", "int.ROK6",
                                       "int.ROK7", "int.ROK8", "int.ROK9", "int.ROK10")], na.rm = TRUE) *5) / (covs$Intervals.teor*5)
covs[which((is.na(covs$cum.mean.int) | covs$cum.mean.int>0) & covs$smoking_status_cat=="Never"),]$cum.mean.int <- 0 
covs[which((is.na(covs$cum.mean.int2) | covs$cum.mean.int2>0) & covs$smoking_status_cat=="Never"),]$cum.mean.int2 <- 0 

summary(covs$cum.mean.int)
addmargins(table(covs$smoking_status_cat, round(covs$cum.mean.int), exclude= NULL))
hist(covs$cum.mean.int2, main=paste("Cumulative mean intensity", " \n ", "- Including breaks"), xlab="cigarettes per day")

## Cumulative intensity
covs$cum.int <- ifelse(covs$smoking_status_cat=="Never", 0, (rowSums(covs[,c("int.ROK1", "int.ROK2", "int.ROK3", "int.ROK4", "int.ROK5", "int.ROK6",
                                                                             "int.ROK7", "int.ROK8", "int.ROK9", "int.ROK10")], na.rm = TRUE) *5))
boxplot(covs$smoking_status_cat, covs$cum.int, exclude= NULL)
by(covs$cum.int, covs$smoking_status_cat, summary)

## Pack years
covs$packyrsTF <- ((rowSums(covs[,c("int.ROK1", "int.ROK2", "int.ROK3", "int.ROK4", "int.ROK5", "int.ROK6",
                                    "int.ROK7", "int.ROK8", "int.ROK9", "int.ROK10")], na.rm = TRUE) *5)/20) 
covs$packyrsTF <-ifelse(covs$smoking_status_cat=="Never", 0, covs$packyrsTF)
addmargins(table(covs$smoking_status_cat, round(covs$packyrsTF), exclude= NULL))
by(covs$packyrsTF, covs$smoking_status_cat, summary)
plot(covs$packyear, covs$packyrsTF, pch=16)
plot(covs$cum.int/20, covs$packyrsTF, pch=16)

#  Imputing intensity from last questionnaire to blood sample    NOT FINISHED
addmargins(table(covs$smoking_status_cat, covs$diff.s.maxq))
addmargins(table(round(covs[covs$smoking_status_cat=="Former" & covs$tsq_smoking<covs$diff.s.maxq,]$cum.mean.int2), covs[covs$smoking_status_cat=="Former" & covs$tsq_smoking<covs$diff.s.maxq,]$diff.s.maxq))
addmargins(table(round(covs[covs$smoking_status_cat=="Current",]$cum.mean.int2), covs[covs$smoking_status_cat=="Current",]$diff.s.maxq))
addmargins(table(round(covs[covs$smoking_status_cat=="Current",]$mean.int), covs[covs$smoking_status_cat=="Current",]$diff.s.maxq))

covs$last.int <- apply(covs[,c("int.ROK1", "int.ROK2", "int.ROK3", "int.ROK4", "int.ROK5", "int.ROK6",
                               "int.ROK7", "int.ROK8", "int.ROK9", "int.ROK10")], 1, function(x) rev(x[is.finite(x)])[1] )
addmargins(table(round(covs[covs$smoking_status_cat=="Current",]$last.int), covs[covs$smoking_status_cat=="Current",]$diff.s.maxq))

covs$diff.s.stop <- ifelse(covs$smoking_status_cat=="Current",(covs$age.sample - (covs$Stopald.NA*5)), NA)
addmargins(table(covs$diff.s.stop, covs$smoking_status_cat, exclude=NULL))
addmargins(table(covs$diff.s.maxq, covs$smoking_status_cat, exclude=NULL))
addmargins(table(covs[covs$smoking_status_cat=="Current",]$diff.s.maxq, covs[covs$smoking_status_cat=="Current",]$diff.s.stop, exclude=NULL))

covs$last_int_age <- ((covs$Stopald.NA+1)*5 +4)
covs$gap.lastq.b <- ifelse(covs$smoking_status_cat=="Current", (covs$age.sample - covs$last_int_age), NA)
addmargins(table(covs$Stopald.NA, covs$last_int_age, exclude= NULL))
addmargins(table(covs$gap.lastq.b, covs$smoking_status_cat, exclude= NULL))

t <- covs[which(covs$smoking_status_cat=="Current" & covs$gap.lastq.b>3),]
t <- covs[which(covs$smoking_status_cat=="Current" & (covs$age.max.q < covs$last_int_age),]

currents <- covs[covs$smoking_status_cat=="Current",c("int.ROK1", "int.ROK2", "int.ROK3", "int.ROK4", "int.ROK5", "int.ROK6", "int.ROK7", "int.ROK8", "int.ROK9", "int.ROK10")]

ROK$imp.gap.lastq.b <- ifelse(covs$smoking_status_cat=="Current", covs$last.int*
                                
                                covs$mean.int.imp <- rowMeans(ROK[,c("int.ROK1", "int.ROK2", "int.ROK3", "int.ROK4", "int.ROK5", "int.ROK6",
                                                                     "int.ROK7", "int.ROK8", "int.ROK9", "int.ROK10")], na.rm = TRUE) #, "int.ROK11", "int.ROK12"
                              
                              ##
                              addmargins(table(covs$smoking_status_cat, round(covs$duration.max), exclude= NULL))
                              addmargins(table(covs$smoking_status_cat, round(covs$tsq_smoking), exclude= NULL))
                              addmargins(table(covs$smoking_status_cat, round(covs$cum.mean.int), exclude= NULL))
                              addmargins(table(covs$smoking_status_cat, round(covs$mean.int), exclude= NULL))
                              addmargins(table(covs$roykaar, round(covs$duration.max), exclude= NULL))
                              addmargins(table(covs$roykstop, round(covs$age.stop), exclude= NULL))
                              addmargins(table(covs$roykald, round(covs$age.start), exclude= NULL))
                              
                              plot(covs$age.stop, covs$tsq_smoking, pch=16, xlab="Age at stop", ylab="Time since quitting"); text(52,39, paste("r = ", round(cor(covs$age.stop, covs$tsq_smoking, method="s", use="pairwise.complete"), 2)))
                              
                              ### Cumulative smoking index
                              
                              # Comprehensive smoking index
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
                              
                              morecovs <- readRDS("/home/tno010/Desktop/NOWAC DNAm Therese/Analyser/Reproductive factors/lung_BMI_education.rds")
                              covs <- merge(covs, morecovs, by.x="labnr", by.y="labnr") 
                              
                              saveRDS(covs, file="/home/tno010/Desktop/NOWAC DNAm Therese/Analyser/output/Lung_covs_040416.rds")  # 261 women, 260 vars
                              
                              covs <- readRDS("/home/tno010/Desktop/NOWAC DNAm Therese/Analyser/output/Lung_covs_040416.rds")  
                              
                              # save only the smoking variables
                              saveRDS(covs[,c("labnr", "sampleID", "smoking_status_cat", "age.start", "age.stop", "passiv_child", "passiv_adult", "passiv_ad_int", "passiv_ad_int_cat", "CSI", "duration.int", "duration.max", "duration.adj", "tsq_smoking", "cum.mean.int", "packyrsTF", "cum.int", "mean.int")],
                                      file="/home/tno010/Desktop/NOWAC DNAm Therese/Analyser/output/Lung_calc_smoking_vars_080316.rds") 
                              
                              
                              
                              t[,c("labnr", "smoking.status","age.sample", "age.que1", "age.quey", "age.quez", "age.start", "age.stop", "Intervals.actu", "Intervals.teor", "Startald.NA", "Stopald.NA", "roykaar", "duration.max", "duration.adj")]
                              t[,c("labnr", "smoking.status","age.sample", "age.que1", "age.quey", "age.quez", "ROYKNAA", "yROYKNAA", "ZROYKNAA", "age.start", "roykald","SIGALDER", "YSIGALDER", "ZSIGALDER", "roykant1019", "ROYKANT2029","ROYKANT1014", "ROYKANT1519")]
                              t[,c("labnr", "smoking.status","age.sample", "age.start", "roykald", "age.stop", "stopald", "yROYKSTOP", "ZROYKSTOP", "Intervals.actu", "Intervals.teor", "Startald.NA", "Stopald.NA", "roykaar", "duration.max", "duration.adj", "diff.int")]
                              
                              
                              