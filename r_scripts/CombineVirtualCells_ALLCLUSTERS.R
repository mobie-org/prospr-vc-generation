# Copy of CombineVirtualCells.R to process ALL virtual cells for automatic curation

# Hernando M. Vergara. August 2019

# Virtual cells have been selected manually (bad ones are removed) by looking at the images generated
# with PullClustersFromRegions.R

# This script will load, for each region, the dataframe containing the virtual cells, and will select those
# that are remaining in the folder of pictures. It might rename them to include information for the region

##############

MainDir = '/Users/herny/Desktop/EMBL/ProSPr/PlatyBrowser/VirtualCells/GenerationOfVirtualCells/npix3/'
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


####### make a matrix with expression information for each Virtual Cell ####
print('reading data...')
#load the npix3 matrix
SVprofile=read.table(paste(MainDir, 'SuperVoxels_npix3_GeneExpression.txt', sep=''),sep=",",header=T)
# Hard coded location for a diccionary for gene list change
DictFile = '/Users/herny/Desktop/EMBL/ProSPr/PlatyBrowser/VirtualCells/GenerationOfVirtualCells/OlegGenesDiccionary.txt'
# File with genes to remove
listtoremoveFile = '/Users/herny/Desktop/EMBL/ProSPr/PlatyBrowser/VirtualCells/GenerationOfVirtualCells/GenesToRemove.txt'
# load file
listtoremove = read.table(listtoremoveFile,header=F)[,1]


####Rename Columns for Oleg's genes#####
dicc = read.table(DictFile, sep=',',header=F,row.names=1)
for(i in 1:ncol(SVprofile)){
  if(colnames(SVprofile)[i] %in% rownames(dicc)){
    colnames(SVprofile)[i] <- paste(colnames(SVprofile)[i],dicc[colnames(SVprofile)[i],1],sep="--")
  }
}


####Rename genes if they were wrong:
colnames(SVprofile)[which(colnames(SVprofile) == 'Hox5')] <- 'Hox4'

#Make a cluster matrix
#consider clusters as entities and assign gene expression to them
#initialize matrix
print('calculating expression...')
AllClust = matrix(NA,nrow = nrow(ClustDF), ncol = ncol(SVprofile))
colnames(AllClust) = colnames(SVprofile)
rownames(AllClust) = ClustDF$clustersID
for(i in c(1:nrow(AllClust))){
  SuperVoxel_List = strsplit(as.character(ClustDF$supervoxelsID[i]),",")[[1]]
  #Now we calculate the coverage of each gene in the group
  AllClust[i,] = colMeans(SVprofile[SuperVoxel_List,])
}


#Make zero those values that are smaller a certain threshold, to remove the noise
ToZeroThr <- 0.5
AllClust[AllClust < ToZeroThr] <- 0
#Remove empty rows and columns
AllClust = AllClust[,colSums(AllClust)>0]
AllClust = AllClust[rowSums(AllClust)>0,]

#Save
VCexp <- AllClust
save(VCexp, file = "VirtualCells_ALL_ExpressionProfile")
print('done')