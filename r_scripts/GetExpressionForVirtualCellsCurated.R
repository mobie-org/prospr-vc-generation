## This is for creating expression for every virtual cell out of Valentina's final curation ####

######
MainDir = '/Users/herny/Desktop/EMBL/ProSPr/PlatyBrowser/VirtualCells/GenerationOfVirtualCells/npix3/'
setwd(MainDir)

## Get the list from Valentyna and transform it into a dataframe:
# ClustDF, with columns clustersID and supervoxelsID
CuratedFile = 'CellModels_ALL_coordinates_clust_curated.tsv'
curated_vcs = read.table(file=CuratedFile, sep='\t')
#load 3D coordinates:
coord3d=read.table(paste(MainDir, 'SuperVoxels_npix3_SpatialInfo.txt', sep=''),sep=",",header=T)

# parse it
# get the supervoxelID out of the x,y,z coordinates
getSVid <- function(SVcoord,coord3d_file){
  # SVcoord is in the form ( y, x, z )
  parsed = as.numeric(strsplit(gsub("[\\(\\)]", "", SVcoord), ",")[[1]])
  x = parsed[2]
  y = parsed[1]
  z = parsed[3]
  SVidx = which(coord3d_file$Xpos==x & coord3d_file$Zpos==z & coord3d_file$Ypos==y)
  if(length(SVidx)>1){
    print('More than one SV found')
  }
  SVID = rownames(coord3d_file)[SVidx]
  return (SVidx)
}

# get a list of coordinates
getCoordinatesList <- function(coordlist){
  clist = strsplit(as.character(coordlist),";")[[1]]
  return(clist)
}

SVsuperlist = vector()
for(i in 1:nrow(curated_vcs)){
  SVlist <- vector()
  VCname = curated_vcs[i,1]
  rawlist = curated_vcs[i,2]
  SVcoord_list = getCoordinatesList(rawlist)
  for(coord in SVcoord_list){
    SVid = getSVid(coord, coord3d)
    SVlist <- c(SVlist, SVid)
  }
  SVsuperlist <- c(SVsuperlist, paste(SVlist, collapse = ","))
  print(i)
}

### create dataframe
print('Creating dataframe')
ClustDF <- data.frame(
  clusterID <- curated_vcs[,1],
  supervoxelsID <- SVsuperlist
)
names(ClustDF) <- c("clustersID","supervoxelsID")

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


####Rename Columns for Oleg's genes
dicc = read.table(DictFile, sep=',',header=F,row.names=1)
for(i in 1:ncol(SVprofile)){
  if(colnames(SVprofile)[i] %in% rownames(dicc)){
    colnames(SVprofile)[i] <- paste(colnames(SVprofile)[i],dicc[colnames(SVprofile)[i],1],sep="--")
  }
}


####Rename genes if they were wrong:
colnames(SVprofile)[which(colnames(SVprofile) == 'Hox5')] <- 'Hox4'
colnames(SVprofile)[which(colnames(SVprofile) == 'irx')] <- 'irx6'


#############
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
  AllClust[i,] = colMeans(SVprofile[as.numeric(SuperVoxel_List),])
  print(i)
}


#Make zero those values that are smaller a certain threshold, to remove the noise
ToZeroThr <- 0.5
AllClust[AllClust < ToZeroThr] <- 0
#Remove empty rows and columns
AllClust = AllClust[,colSums(AllClust)>0]
AllClust = AllClust[rowSums(AllClust)>0,]

#Save
write.table(AllClust, file='CellModels_ALL_profile_clust_curated.tsv', quote=FALSE, sep='\t')
print('done')
