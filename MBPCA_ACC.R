# MBPCA for Accelorometer data
# written by Chen Yue on 05/01/2014

# load data
load('C:/Users/Chen/Dropbox/Research/PCA_working_group/code/Data_Jiawei/NewDataNewSub.RData')
load("NewDataNewSub.RData")
valid_id <- NULL
for(i in 1:dim(ActivityCount)[1]){
  if(sum(is.na(ActivityCount[i,]))==0){
    valid_id <- c(valid_id, i)
  }
}

Count_Valid <- ActivityCount[valid_id,]
Info_Valid <- SubjectInfo[valid_id,]
J <- table(Info_Valid$BLSA_id)
valid_id <- NULL
for(i in 1:dim(ActivityCount)[1]){
  if((sum(is.na(ActivityCount[i,]))==0)&(!(as.character(SubjectInfo$BLSA_id[i]) %in% names(which(J==1))))){
    valid_id <- c(valid_id, i)
  }
}
Count_Valid <- ActivityCount[valid_id,]
Info_Valid <- SubjectInfo[valid_id,]
J <- table(Info_Valid$BLSA_id)
Count_Valid[Count_Valid>0] <- 1
source('MBPCA.R')
fit <- MBPCA(X=Count_Valid, J=J, k_between=5, k_within=5, n_iter=30)
save(fit, file="MBPCA_AccDat.rda")