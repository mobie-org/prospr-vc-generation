# Modification of PullClustersFromRegions.R to get ALL virtual cells for automatic curation.


# Get the raw Virtual Cells from each region (created with 'Clusters_Extraction_6dpf_npix3.R')
# to filter them by size and save plots of each of them for manual curation.

# This has to be implemented because I did not solve the rgl xfvb issue in the slurm cluster


library(ggplot2)
library(png)
library(rgl)

args = commandArgs(trailingOnly=TRUE)
# test if number of arguments are correct: if not, return an error
if (length(args)<1) {
  stop("Script called incorrectly.
       Please provide:
       Output folder of ProSPr_6dpf_SuperVoxelPixCount\n", call.=FALSE)
}

HelperFunctionsFile = './scripts/Clusters_Extraction_functions.R'
source(HelperFunctionsFile)
MainDir = args[1]
setwd(MainDir)
#########get the supervoxels coordinate file
SVcoordinates=read.table('./SuperVoxels_npix3_SpatialInfo.txt',sep=",",header=T)

# get a list of the directories
filesInDir = list.files(path = ".", pattern = 'RegPro*')

for(folder in filesInDir[1:length(filesInDir)]){
  setwd(file.path(MainDir, folder))
  print(folder)
  #########get the SV matrix
  inputfilename <- list.files(path = '.', pattern = '*_SupervoxelsExpression_NoCorrRemoved')[1]
  load(inputfilename)
  initMat <- Regionvoxels_reduced

  #### get the coordinates for the region
  initCoord <- SVcoordinates[as.character(rownames(initMat)),]

  ### get the raw VCs
  inputfilename <- list.files(path = '.', pattern = '*_VirtualCells_Raw')[1]
  load(inputfilename)

  ### create dataframe
  print('Creating dataframe')
  clustersDF <- data.frame(
    clusterID <- names(finalclusters),
    supervoxelsID <- as.character(lapply(finalclusters,function(x){paste(x,collapse=",")})),
    correlation <- as.numeric(lapply(finalclusters,function(x){MyCorrFunc(x,initMat)})),
    spatialCoherence <- as.numeric(lapply(finalclusters,function(x){SpatCoheFunc(x,initCoord)})),
    clusterSize <- as.numeric(lapply(finalclusters,length))
  )
  names(clustersDF) <- c("clustersID","supervoxelsID","correlation","spatialCoherence","clusterSize")

  save(clustersDF, file = paste(folder,'_VirtualCells_ALL',sep=''))
  print('data saved')
}
