library(FortMyrmidon) ####R bindings
library(Rcpp)
library(circular)
library(R.utils)
rm(list=ls())
gc()
dir_data <- "/media/ll16598/SeagateDesktopDrive/SICC_DATA/"
dir_SI <- "/media/ll16598/One Touch/SICC_METADATA/"
dir_chall <- "/media/ll16598/One Touch/SICC_METADATA/CHALL_ID/"

myrmidon_file1 <- paste("/media/ll16598/SeagateDesktopDrive/SICC_DATA/2S_IV58_090222_WST/LLEXPP2S_IV58_090222_WST_KEY_DEATHS_AutoOriented.myrmidon",sep='')
Metadata_exp <- fmExperimentOpen(myrmidon_file1)

#FOR BEAD CHALLENGED COLONIES
BLUE_ <- paste(dir_chall, "2S_NB.myrmidon",sep='') #Sham SI IDs
BLUE <- fmExperimentOpen(BLUE_)
YELLOW_ <- paste(dir_chall, "2S_LY.myrmidon",sep='') #Sham SI IDs
YELLOW <- fmExperimentOpen(YELLOW_)

#SIBB_ <- paste("/media/ll16598/One Touch/LUKE/SI_ID/LLEXPP1S_010222_HDN/LLEXPP1S_IV55_010222_SIBB.myrmidon",sep='') #Beauveria SI IDs
#SIKVL_ <- paste("/media/ll16598/One Touch/LUKE/SI_ID/LLEXPP1S_010222_HDN/LLEXPP1S_IV55_010222_SIKVL1.myrmidon",sep='') #Beauveria SI IDs
#SISH_ <- paste("/media/ll16598/One Touch/LUKE/SI_ID/LLEXPP1S_010222_HDN/LLEXPP1S_IV55_010222_SIS.myrmidon",sep='') #Beauveria SI IDs

#READ SI FILES
SIBB_ <- paste(dir_SI, "2S_SIBB.myrmidon",sep='') #Beauveria SI IDs
SIBB <- fmExperimentOpen(SIBB_)
SIKVL_ <- paste(dir_SI, "2S_SIKVL.myrmidon",sep='') #Metarhizium SI IDs
SIKVL <- fmExperimentOpen(SIKVL_)
SISH_ <- paste(dir_SI, "2S_SISH.myrmidon",sep='') #Sham SI IDs
SISH <- fmExperimentOpen(SISH_)


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
BLUE_CHALL_LIST <- NULL
for (ant in BLUE_ants){
  tag_val <- ant$identifications[[1]]$tagValue
  print(tag_val)
  BLUE_CHALL_LIST <- rbind(BLUE_CHALL_LIST, data.frame(hex_id = tag_val, stringsAsFactors = F))
}
setDT(BLUE_CHALL_LIST)
write.table(BLUE_CHALL_LIST, file = paste0(dir_chall, "/BLUE_CHALL_list",basename(BLUE_), ".txt", sep=""), append = FALSE,
            row.names = FALSE, col.names = TRUE,quote = FALSE)

YELLOW_tag_statistics <- fmQueryComputeTagStatistics(YELLOW) #do for myrmidon file also
YELLOW_IDs <- YELLOW_tag_statistics$tagDecimalValue
for ( i in 1:nrow(YELLOW_tag_statistics)) {  #####loop over each tag
  #if ( tag_statistics[i,"count"] >= 0.001*max(tag_statistics[,"count"],na.rm=T) ) { ### optional: here I decide to create an antID only if the tag detection rate was more than 1/1000 * the best tag detection rate. You can choose your own criterion
  a <- YELLOW$createAnt(); ###this actually creates an antID, i.e. associates a decimal antID number to that particular tagID
  identification <- YELLOW$addIdentification(a$ID,YELLOW_tag_statistics[i,"tagDecimalValue"],fmTimeSinceEver(),fmTimeForever())
  print(identification)
}
YELLOW_ants <- YELLOW$ants
YELLOW_CHALL_LIST <- NULL
for (ant in YELLOW_ants){
  tag_val <- ant$identifications[[1]]$tagValue
  print(tag_val)
  YELLOW_CHALL_LIST <- rbind(YELLOW_CHALL_LIST, data.frame(hex_id = tag_val, stringsAsFactors = F))
}
setDT(YELLOW_CHALL_LIST)
write.table(YELLOW_CHALL_LIST, file = paste0(dir_chall, "/YELLOW_CHALL_list",basename(YELLOW_), ".txt", sep=""), append = FALSE,
            row.names = FALSE, col.names = TRUE,quote = FALSE)


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
SIBB_LIST <- NULL
for (ant in SIBB_ants){
  tag_val <- ant$identifications[[1]]$tagValue
  print(tag_val)
  SIBB_LIST <- rbind(SIBB_LIST, data.frame(hex_id = tag_val, stringsAsFactors = F))
}
setDT(SIBB_LIST)
write.table(SIBB_LIST, file = paste0(dir_SI, "/SIBB_list",basename(SIBB_), ".txt", sep=""), append = FALSE,
            row.names = FALSE, col.names = TRUE,quote = FALSE)
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
SIKVL_LIST <- NULL
for (ant in SIKVL_ants){
  tag_val <- ant$identifications[[1]]$tagValue
  print(tag_val)
  SIKVL_LIST <- rbind(SIKVL_LIST, data.frame(hex_id = tag_val, stringsAsFactors = F))
}
setDT(SIKVL_LIST)
write.table(SIKVL_LIST, file = paste0(dir_SI, "/SIKVL_list",basename(SIKVL_), ".txt", sep=""), append = FALSE,
            row.names = FALSE, col.names = TRUE,quote = FALSE)

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
SISH_LIST <- NULL
for (ant in SISH_ants){
  tag_val <- ant$identifications[[1]]$tagValue
  print(tag_val)
  SISH_LIST <- rbind(SISH_LIST, data.frame(hex_id = tag_val, stringsAsFactors = F))
}


write.table(SISH_LIST, file = paste0(dir_SI, "/SISH_list",basename(SISH_), ".txt", sep=""), append = FALSE,
            row.names = FALSE, col.names = TRUE,quote = FALSE)
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
  for (s in SISH$ants){
    if (a$identifications[[1]]$tagValue == s$identifications[[1]]$tagValue){
      a$setValue("SISH",TRUE, t)}}
}

#CONTROL BEAD ANT METADATA
for (a in Metadata_exp$ants){   ####LOOP OVER ANTS
  ###get corresponding ant from metadata
  # e <- a$createAnt()
  for (N in BLUE$ants){
    if (a$identifications[[1]]$tagValue == N$identifications[[1]]$tagValue){
      a$setValue("CHALL_BLUE",TRUE, t)}}
  for (Y in YELLOW$ants){
    if (a$identifications[[1]]$tagValue == Y$identifications[[1]]$tagValue){
      a$setValue("CHALL_YELLOW",TRUE, t)}}
}

Metadata_exp$save(paste0(sub("\\..*", "", myrmidon_file1),"_withMetaData.myrmidon"))  
