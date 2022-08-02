library(FortMyrmidon) ####R bindings
library(Rcpp)
library(circular)
library(R.utils)

dir_data <- "/media/ll16598/Seagate Desktop Drive/"
dir_SI <- "/home/ll16598/Documents/SICC_METADATA/SI_ID/13P/"
dir_chall <- "/home/ll16598/Documents/SICC_METADATA/Pathogen_ID/"

myrmidon_file1 <- paste(dir_data,"13P_IV31_110522_HDN/13P_IV31_110522_HDN_manual_orient_test.myrmidon",sep='')
Metadata_exp <- fmExperimentOpen(myrmidon_file1)

#FORPATHOGEN CHALLENGED COLONIES
CHALL_KVL_ <- paste(dir_chall, ".myrmidon",sep='') #Sham SI IDs
CHALL_KVL <- fmExperimentOpen(CHALL_KVL_)
CHALL_BB_ <- paste(dir_chall, ".myrmidon",sep='') #Sham SI IDs
CHALL_BB <- fmExperimentOpen(CHALL_BB_)

#PATHOGEN CHALLENGED IDS
BB_tag_statistics <- fmQueryComputeTagStatistics(CHALL_BB) #do for myrmidon file also
BB_IDs <- BB_tag_statistics$tagDecimalValue
for ( i in 1:nrow(BB_tag_statistics)) {  #####loop over each tag
  #if ( tag_statistics[i,"count"] >= 0.001*max(tag_statistics[,"count"],na.rm=T) ) { ### optional: here I decide to create an antID only if the tag detection rate was more than 1/1000 * the best tag detection rate. You can choose your own criterion
  a <- BB$createAnt(); ###this actually creates an antID, i.e. associates a decimal antID number to that particular tagID
  identification <- CHALL_BB$addIdentification(a$ID,BB_tag_statistics[i,"tagDecimalValue"],fmTimeSinceEver(),fmTimeForever())
  print(identification)
}
BB_CHALL_ants <- CHALL_BB$ants #BEAUVERIA IDS

KVL_tag_statistics <- fmQueryComputeTagStatistics(KVL_CHALL) #do for myrmidon file also
KVL_IDs <- KVL_tag_statistics$tagDecimalValue
for ( i in 1:nrow(KVL_tag_statistics)) {  #####loop over each tag
  #if ( tag_statistics[i,"count"] >= 0.001*max(tag_statistics[,"count"],na.rm=T) ) { ### optional: here I decide to create an antID only if the tag detection rate was more than 1/1000 * the best tag detection rate. You can choose your own criterion
  a <- KVL_CHALL$createAnt(); ###this actually creates an antID, i.e. associates a decimal antID number to that particular tagID
  identification <- KVL_CHALL$addIdentification(a$ID,KVL_tag_statistics[i,"tagDecimalValue"],fmTimeSinceEver(),fmTimeForever())
  print(identification)
}
KVL_CHALL_ants <- KVL_CHALL$ants

#FOR BEAD CHALLENGED COLONIES
BLUE_ <- paste(dir_chall, ".myrmidon",sep='') #Sham SI IDs
BLUE <- fmExperimentOpen(BLUE_)
YELLOW_ <- paste(dir_chall, ".myrmidon",sep='') #Sham SI IDs
YELLOW <- fmExperimentOpen(YELLOW_)

#BEAD CHALLENGE IDS
BLUE_tag_statistics <- fmQueryComputeTagStatistics(BLUE) #do for myrmidon file also
BLUE_IDs <- BLUE_tag_statistics$tagDecimalValue
for ( i in 1:nrow(BLUE_tag_statistics)) {  #####loop over each tag
  #if ( tag_statistics[i,"count"] >= 0.001*max(tag_statistics[,"count"],na.rm=T) ) { ### optional: here I decide to create an antID only if the tag detection rate was more than 1/1000 * the best tag detection rate. You can choose your own criterion
  a <- BLUE$createAnt(); ###this actually creates an antID, i.e. associates a decimal antID number to that particular tagID
  identification <- BLUE$addIdentification(a$ID,BLUE_tag_statistics[i,"tagDecimalValue"],fmTimeSinceEver(),fmTimeForever())
  print(identification)
}
BLUE_ants <- BLUE$ants #BLUE ANT IDS

YELLOW_tag_statistics <- fmQueryComputeTagStatistics(YELLOW) #do for myrmidon file also
YELLOW_IDs <- YELLOW_tag_statistics$tagDecimalValue
for ( i in 1:nrow(YELLOW_tag_statistics)) {  #####loop over each tag
  #if ( tag_statistics[i,"count"] >= 0.001*max(tag_statistics[,"count"],na.rm=T) ) { ### optional: here I decide to create an antID only if the tag detection rate was more than 1/1000 * the best tag detection rate. You can choose your own criterion
  a <- YELLOW$createAnt(); ###this actually creates an antID, i.e. associates a decimal antID number to that particular tagID
  identification <- YELLOW$addIdentification(a$ID,YELLOW_tag_statistics[i,"tagDecimalValue"],fmTimeSinceEver(),fmTimeForever())
  print(identification)
}
YELLOW_ants <- YELLOW$ants

#READ SI FILES
SIBB_ <- paste(dir_SI, "13P_SIBB.myrmidon",sep='') #Beauveria SI IDs
SIBB <- fmExperimentOpen(SIBB_)
SIKVL_ <- paste(dir_SI, "13P_SIKVL.myrmidon",sep='') #Metarhizium SI IDs
SIKVL <- fmExperimentOpen(SIKVL_)
SISH_ <- paste(dir_SI, "13P_SISH.myrmidon",sep='') #Sham SI IDs
SISH <- fmExperimentOpen(SISH_)



#CREATE SIBB ANTS
SIBB_tag_statistics <- fmQueryComputeTagStatistics(SIBB) #do for myrmidon file also
SIBB_IDs <- SIBB_tag_statistics$tagDecimalValue
for ( i in 1:nrow(SIBB_tag_statistics)) {  #####loop over each tag
  #if ( tag_statistics[i,"count"] >= 0.001*max(tag_statistics[,"count"],na.rm=T) ) { ### optional: here I decide to create an antID only if the tag detection rate was more than 1/1000 * the best tag detection rate. You can choose your own criterion
    a <- SIBB$createAnt(); ###this actually creates an antID, i.e. associates a decimal antID number to that particular tagID
    identification <- SIBB$addIdentification(a$ID,SIBB_tag_statistics[i,"tagDecimalValue"],fmTimeSinceEver(),fmTimeForever())
    print(identification)
  }
SIBB_ants <- SIBB$ants

#CREATE SIKVL ANTS
SIKVL_tag_statistics <- fmQueryComputeTagStatistics(SIKVL) #do for myrmidon file also
SIKVL_IDs <- SIKVL_tag_statistics$tagDecimalValue
for ( i in 1:nrow(SIKVL_tag_statistics)) {  #####loop over each tag
  #if ( tag_statistics[i,"count"] >= 0.001*max(tag_statistics[,"count"],na.rm=T) ) { ### optional: here I decide to create an antID only if the tag detection rate was more than 1/1000 * the best tag detection rate. You can choose your own criterion
  a <- SIKVL$createAnt(); ###this actually creates an antID, i.e. associates a decimal antID number to that particular tagID
  identification <- SIKVL$addIdentification(a$ID,SIKVL_tag_statistics[i,"tagDecimalValue"],fmTimeSinceEver(),fmTimeForever())
  print(identification)
}
SIKVL_ants <- SIKVL$ants

#CREATE SISH ANTS
SISH_tag_statistics <- fmQueryComputeTagStatistics(SISH) #do for myrmidon file also
SISH_IDs <- SISH_tag_statistics$tagDecimalValue
for ( i in 1:nrow(SISH_tag_statistics)) {  #####loop over each tag
  #if ( tag_statistics[i,"count"] >= 0.001*max(tag_statistics[,"count"],na.rm=T) ) { ### optional: here I decide to create an antID only if the tag detection rate was more than 1/1000 * the best tag detection rate. You can choose your own criterion
  a <- SISH$createAnt(); ###this actually creates an antID, i.e. associates a decimal antID number to that particular tagID
  identification <- SISH$addIdentification(a$ID,SISH_tag_statistics[i,"tagDecimalValue"],fmTimeSinceEver(),fmTimeForever())
  print(identification)
}
SISH_ants <- SISH$ants

#CREATE ANTS FOR MAIN EXP FILE
Metadata_ants <- Metadata_exp$ants
t <- fmTimeCreate(offset = 0)#SET TIME TO 1970

##this changes metadata if tag IDs match for SI ants 02/08/22
for (a in Metadata_exp$ants){   ####LOOP OVER ANTS
  ###get corresponding ant from metadata
  # e <- a$createAnt()
  for (b in SIBB$ants){
    if (a$identifications[[1]]$tagValue == b$identifications[[1]]$tagValue){
      a$setValue("SIBB",TRUE, t)}}
  for (k in SIKVL$ants){
    if (a$identifications[[1]]$tagValue == k$identifications[[1]]$tagValue){
      a$setValue("SIKVL",TRUE, t)}}
  for (s in SISHAM$ants){
    if (a$identifications[[1]]$tagValue == s$identifications[[1]]$tagValue){
      a$setValue("SISH",TRUE, t)}}
}
# changes metadata for challenge ants
for (a in Metadata_exp$ants){   ####LOOP OVER ANTS
  ###get corresponding ant from metadata
  # e <- a$createAnt()
  for (B in BB_CHALL$ants){
    if (a$identifications[[1]]$tagValue == B$identifications[[1]]$tagValue){
      a$setValue("BB_CHALL",TRUE, t)}}
  for (K in KVL_CHALL$ants){
    if (a$identifications[[1]]$tagValue == K$identifications[[1]]$tagValue){
      a$setValue("KVL_CHALL",TRUE, t)}}
}
#CONTROL BEAD ANT METADATA
for (a in Metadata_exp$ants){   ####LOOP OVER ANTS
  ###get corresponding ant from metadata
  # e <- a$createAnt()
  for (N in BLUE$ants){
    if (a$identifications[[1]]$tagValue == N$identifications[[1]]$tagValue){
      a$setValue("BLUE",TRUE, t)}}
  for (Y in YELLOW$ants){
    if (a$identifications[[1]]$tagValue == Y$identifications[[1]]$tagValue){
      a$setValue("YELLOW",TRUE, t)}}
}

Metadata_exp$save(paste0(sub("\\..*", "", myrmidon_file1),"_withMetaData.myrmidon"))  
