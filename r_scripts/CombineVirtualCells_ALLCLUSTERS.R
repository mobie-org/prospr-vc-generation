# Copy of CombineVirtualCells.R to process ALL virtual cells for automatic curation

# Hernando M. Vergara. August 2019

# Virtual cells have been selected manually (bad ones are removed) by looking at the images generated
# with PullClustersFromRegions.R

# This script will load, for each region, the dataframe containing the virtual cells, and will select those
# that are remaining in the folder of pictures. It might rename them to include information for the region

##############

args = commandArgs(trailingOnly=TRUE)
# test if number of arguments are correct: if not, return an error
if (length(args)<1) {
  stop("Script called incorrectly.
       Please provide:
       Output folder of ProSPr_6dpf_SuperVoxelPixCount\n", call.=FALSE)
}

MainDir = args[1]
setwd(MainDir)
# get a list of the directories
filesInDir = list.files(path = ".", pattern = 'RegPro*')
# for naming SV
generalCounter = 1
PreviousShortRegionName = FALSE
# for storing dfs
DFlist = list()

# go through each region
for(folder in filesInDir[1:length(filesInDir)]){
  setwd(file.path(MainDir, folder))
  print(folder)
  #########get the SV dataframe
  inputfilename <- list.files(path = '.', pattern = '*_VirtualCells_ALL')[1]
  load(inputfilename)


  clustersHQ <- clustersDF

  #add new fields with easier IDs
  #rename old name
  names(clustersHQ)[names(clustersHQ) == "clustersID"] <- "HCidentifier"
  #add a general number to the Clusters ID for easier recognition
  name_vec <- rep(NA,nrow(clustersHQ))
  for(i in c(1:nrow(clustersHQ))){
    name_vec[i] <- sprintf("%04.f", generalCounter)
    generalCounter <- generalCounter + 1
  }
  clustersHQ$clusterNumber <- name_vec
  #add a field with the region they are coming from
  clustersHQ$OriginRegion <- strsplit(folder,'RegPro_')[[1]][2]

  #add a unique ID for each cluster with the region origin
  # for those that the names to give are the same (e.g. Head), put them one after the other
  #CHANGE THIS FOR EACH REGION
  if(folder == "RegPro_All-MHCL4-Positive"){ShortRegionName <- 'MHCL4'}
  if(folder == "RegPro_CrypticSegment"){ShortRegionName <- 'CS'}
  if(folder == "RegPro_Head-Lmo4-Negative"){ShortRegionName <- 'Head'}
  if(folder == "RegPro_Head-Lmo4-Positive"){ShortRegionName <- 'Head'}
  if(folder == "RegPro_PNS"){ShortRegionName <- 'PNS'}
  if(folder == "RegPro_Pygidium"){ShortRegionName <- 'Pyg'}
  if(folder == "RegPro_RestOfAnimal"){ShortRegionName <- 'ROA'}
  if(folder == "RegPro_Stomodeum"){ShortRegionName <- 'Stom'}
  if(folder == "RegPro_VNC"){ShortRegionName <- 'VNC'}

  name_vec <- rep(NA,nrow(clustersHQ))
  # restart counter if region changes name
  if(ShortRegionName!=PreviousShortRegionName){
    SRNcounter <- 0
  }
  # name
  for(i in c(1:nrow(clustersHQ))){
    SRNcounter <- SRNcounter + 1
    name_vec[i] <- paste(ShortRegionName,"_",sprintf("%04.f", SRNcounter),sep = '')
  }
  clustersHQ$clustersID <- name_vec

  # add it to the dataframes list
  DFlist <- append(DFlist, list(clustersHQ))

  # for not start numbering from 0 if the region name is the same
  PreviousShortRegionName <- ShortRegionName
}

# merge and save
ClustDF <- data.frame(DFlist[1])
for(i in c(2:length(DFlist))){
  ClustDF <- rbind(ClustDF, data.frame(DFlist[i]))
}
setwd(MainDir)
VCprop <- ClustDF
save(VCprop, file = "VirtualCells_ALL_Properties")
