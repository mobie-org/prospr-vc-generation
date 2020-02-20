# The purpose of this is to create a table that looks like this:

# Cellular model ID          SuperVoxel Coordinates                      Gene1          Gene2    ....
# 000001                     ([x1,y1,z1], [x2,y2,z2],...)                0.6              0.1        ....
# 000002                     ([x1,y1,z1], [x2,y2,z2],...)                0                 1           ....
# 000003                      ([x1,y1,z1], [x2,y2,z2],...)              0.4              0.9           ....

# THEY ARE ACTUALLY Y,X,Z


#working directory
args = commandArgs(trailingOnly=TRUE)
# test if number of arguments are correct: if not, return an error
if (length(args)<1) {
  stop("Script called incorrectly.
       Please provide:
       Output folder of ProSPr_6dpf_SuperVoxelPixCount\n", call.=FALSE)
}

MainDir = args[1]

setwd(MainDir) #***

#SPECIFIC DATASETS

#load clusters data frame (this contains information about the SV in the clusters or cells):
x <- load("VirtualCells_ALL_Properties") #***
Clust_DF <- get(x)
rm(x)


#load 3D coordinates:
coord3d=read.table(paste(MainDir, 'SuperVoxels_npix3_SpatialInfo.txt', sep=''),sep=",",header=T)


# get the coordinates of a SuperVoxel
getSVcoord <- function(SVID,coord3d_file){
  coord <- coord3d_file[toString(SVID),]
  return (coord)
}

###########Get a list of supervoxels in a group of cells#########
#Returns a list of SV grouped by color
getListOfSV <- function(cellnames,cellGroups,cellsSVs){
  #group cells based on the number of clusters
  if(!is.null(cellGroups)){
    colors.names <- unique(cellGroups[,1])
  }else{
    colors.names <- c("red")
  }
  SVL <- vector("list", length(colors.names))
  names(SVL) <- colors.names

  if(!is.null(cellGroups)){
    for (CellName in cellnames){
      col <- cellGroups[CellName,1]
      SVL[[col]] = append(SVL[[col]],strsplit(as.character(cellsSVs[cellsSVs$clustersID==CellName,]$supervoxelsID),",")[[1]])
    }
  }else{
    for (CellName in cellnames){
      col <- "red"
      SVL[[col]] = append(SVL[[col]],strsplit(as.character(cellsSVs[cellsSVs$clustersID==CellName,]$supervoxelsID),",")[[1]])
    }
  }
  return (SVL)
}

Clust_DF$SVcoords = NA

for (i in 1:nrow(Clust_DF)){
  CellModel = Clust_DF[i,]
  SVtest = getListOfSV(CellModel$clustersID,cellGroups = NULL,Clust_DF)$red

  coordsList <- vector("list", length(SVtest))
  coordsList <- setNames(coordsList, SVtest)
  for (SV in SVtest){
    coords <- paste("(",toString(getSVcoord(SV, coord3d)),")", collapse = "")
    coordsList[[SV]] <- coords
  }
  ListCollapsed <- paste(coordsList, collapse = " ; ")
  Clust_DF[i,]$SVcoords = ListCollapsed
  print(i)
}

NewDF <- Clust_DF[c("clustersID", "SVcoords")]

write.table(NewDF, file='CellModels_ALL_coordinates.tsv', quote=FALSE, sep='\t', col.names = F, row.names = F)
