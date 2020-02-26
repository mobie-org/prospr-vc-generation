#Region Profiler for 6dpf
#Profile a particular region in the animal, segmented using TrackEM for example

# variables read from command line.

# Region_Profiler.R
# Output directory of ProSPr_6dpf_SuperVoxelPixCount,
# Name of animal region (can be 'All')
# [Optional gene name] (here, "-Hb9" would take whatever is not Hb9 for example)

# Usage example
# Rscript Region_Profiler.R npix3/ "Pygidium" "Pax6"

#########initiate#######
# Loading it from the script directory, before the directory is set to OutputDir
DictFile = './helper_files/OlegGenesDiccionary.txt'
# load file
dicc = read.table(DictFile, sep=',',header=F,row.names=1)
# File with genes to remove
listtoremoveFile = './helper_files/GenesToRemove.txt'
# load file
listtoremove = read.table(listtoremoveFile,header=F)[,1]

# get inputs
args = commandArgs(trailingOnly=TRUE)
# test if number of arguments are correct: if not, return an error
if (length(args)<2) {
  stop("Script called incorrectly.
       Please provide:
       Output directory of ProSPr_6dpf_SuperVoxelPixCount,
       Name of animal region
       [optional gene name]\n", call.=FALSE)
}
matdir = args[1]
NameOfRegion = args[2]
OutputDirName = paste('RegPro_',NameOfRegion,sep='')
selectBasedOnGene=F
if (length(args)==3){
  selectBasedOnGene=T
  GeneInfo = args[3]
  # parse the information to create appropriate variables
  GIsp = strsplit(GeneInfo,'')[[1]]
  if(GIsp[1]=='-'){
    isneg=T
    GeneToSelect = paste(GIsp[2:length(GIsp)],collapse='')
    OutputDirName = paste(OutputDirName,GeneToSelect,'Negative',sep='-')
  }else{
    isneg=F
    GeneToSelect = paste(GIsp[1:length(GIsp)],collapse='')
    OutputDirName = paste(OutputDirName,GeneToSelect,'Positive',sep='-')
    }
}

setwd(matdir)
#get the number of npix
npix = tail(strsplit(matdir,'/')[[1]],n=1)[1]
# load matrices
SVprofile=read.table(paste('SuperVoxels_',npix,'_GeneExpression.txt', sep=''),sep=",",header=T)
# This one might not be neccessary as it is redundant
#SVcoordinates=read.table(paste('SuperVoxels_',npix,'_SpatialInfo.txt', sep=''),sep=",",header=T)
#load neighbour matrix
load(paste('./SuperVoxels_',npix,'_neighbour_matrix', sep=''))

cat('Data loaded\n')

#check that the gene name and the region match
if(NameOfRegion!='All'){
  if(!(NameOfRegion %in% colnames(SVprofile))){
    stop('Name of the region not found in the data')
  }
}
if(selectBasedOnGene){
  if(!(GeneToSelect %in% colnames(SVprofile))){
    stop('Name of the gene not found in the data')
  }
}

#create directory and move there
dir.create(OutputDirName)
setwd(OutputDirName)
####Rename Columns for Oleg's genes#####
for(i in 1:ncol(SVprofile)){
  if(colnames(SVprofile)[i] %in% rownames(dicc)){
    colnames(SVprofile)[i] <- paste(colnames(SVprofile)[i],dicc[colnames(SVprofile)[i],1],sep="--")
  }
}

###########Region voxels########
cat(paste('\nSelecting the supervoxels corresponding to the ', NameOfRegion, sep=''))
#isolate specific voxels (rows corresponding to voxels in the Region of interest)
if(NameOfRegion!='All'){
  Regionvoxels = SVprofile[SVprofile[,NameOfRegion]>=0.8*max(SVprofile),]
}else{
  Regionvoxels = SVprofile
}
#rownames(Regionvoxels)=which(SVprofile[,nameofregion]>=0.8*max(SVprofile)) #keep the index of voxels

#Select the SV based on gene information
if(selectBasedOnGene==T){
  cat(paste('\nSelecting the SV based on gene (', GeneToSelect , ')information\n', sep = ''))
  if(isneg){
    Regionvoxels <- Regionvoxels[which(Regionvoxels[,GeneToSelect]<(max(Regionvoxels)/2)),]
  }else{
    Regionvoxels <- Regionvoxels[which(Regionvoxels[,GeneToSelect]>=(max(Regionvoxels)/2)),]
  }
}

cat('\nCleaning the data')
#remove those columns that have no expression:
Regionvoxels = Regionvoxels[,colSums(Regionvoxels)>max(Regionvoxels)]
#remove genes in the list
Regionvoxels = Regionvoxels[,!(colnames(Regionvoxels) %in% listtoremove)]
#remove empty rows
Regionvoxels = Regionvoxels[rowSums(Regionvoxels)>max(Regionvoxels),]

#Regioncoordinates = SVcoordinates[rownames(Regionvoxels),]

#This might not even be necessary
#save(Regionvoxels,file = paste(OutputDirName,'_SupervoxelsExpression',sep=""))
#The coordinates should not be necessary to compute or save
#save(Regioncoordinates,file = paste(OutputDirName,'_SupervoxelsCoordinates',sep=""))

#cat('\nRegion-specific tables saved\n')

########removal of low correlation######
cat('Removing low-correlation supervoxels\n')
###New way of doing it with the neighbour matrix
#Calculate the correlation of each voxel just with the neighbours
NeighCorr <- function(SVID){ #SuperVoxel ID as a parameter
  #Neighbour size
  npixsp = strsplit(npix,'')[[1]]
  SVSize = strtoi(paste(npixsp[5:length(npixsp)],collapse=''))
  NeSize = SVSize^3 - 1
  #Values of the SV
  SVvals <- Regionvoxels[as.character(SVID),]
  #submatrix with only the neighbours
  NeighList <- neighbour_mat[as.character(SVID),]
  #remove those that are not in the dataset
  NeighList_In <- NeighList[NeighList %in% rownames(Regionvoxels)]
  submat <- Regionvoxels[NeighList_In,,drop=F]
  if (nrow(submat)>=1){ #if there are any neighbours
    #remove columns that have 0s
    submat_clean <- submat[,(colSums(submat)+SVvals)>0,drop=F]
    SVvals_clean <- SVvals[(colSums(submat)+SVvals)>0]
    #calculate the correlations
    corvec <- cor(t(submat_clean),t(t(SVvals_clean)))[,1]
  }else{corvec <- rep(NA,NeSize)}
  #fill with NAs till the end to get a matrix as an output
  if(length(corvec)<NeSize){
    corvec <- append(corvec,rep(NA,NeSize-length(corvec)))
  }
  #return the list of correlations
  corvec
  #SVvals_clean
}

CorrMat <- t(apply(as.matrix(rownames(Regionvoxels)),1,NeighCorr))
rownames(CorrMat) <- rownames(Regionvoxels)
save(CorrMat,file = paste(OutputDirName,"_SV_NeighbourCorrelationMatrix",sep=""))
cat('\nNeighbour correlation matrix saved\n')

########

cat('Making the data small enough for R to handle\n')

#mean vs std
MeansVec <- apply(CorrMat,1,function(x){mean(x,na.rm=T)})
SdsVec <- apply(CorrMat,1,function(x){sd(x,na.rm=T)})

CorrResDf <- data.frame(MeansVec,SdsVec)
#remove missing values
CorrResDf <- CorrResDf[complete.cases(CorrResDf),]
CorrResDf$Label <- rep("In",nrow(CorrResDf))

# TODO: keep removing SV until the data is small enough to be processed
MeanCO <- min(CorrResDf$MeansVec) + 0.01
while(length(rownames(CorrResDf[CorrResDf$Label == "In",]))^2>.Machine$integer.max){
  MeanCO = MeanCO + 0.01
  CorrResDf[CorrResDf$MeansVec < MeanCO, ]$Label <- "Out"
  cat(paste('\nTesting value ',MeanCO,sep=''))
}
cat(paste('\n\nValue to remove supervoxels based on neighbour correlation selected to ', MeanCO, '\n', sep=''))
Regionvoxels_reduced <- Regionvoxels[rownames(CorrResDf[CorrResDf$Label == "In",]),]
cat(paste('Number of supervoxels: ', nrow(Regionvoxels_reduced),sep=''))
#Save matrix
save(Regionvoxels_reduced,file = paste(OutputDirName,"_SupervoxelsExpression_NoCorrRemoved",sep=""))
cat('\nData saved, plotting and wrapping up\n')
###############plot

#library(ggplot2)

#p1 <- ggplot(CorrResDf, aes(x=MeansVec, y=SdsVec))
#p1 +
#  geom_point(alpha = 0.1,size = 2, aes(colour = Label)) +
#  #geom_point(alpha = 0.1,size = 2) +
#  scale_color_manual(values = c("black","blue")) +
#  #geom_path(data=ell, aes(x=x, y=y), size=1, linetype=1, color="red") +
#  stat_density_2d(aes(fill = ..level.., alpha = ..level..),size = 0.2, geom = 'polygon') +
#  scale_fill_gradient(low = "green", high = "red") +
#  scale_alpha(range = c(0.2, 0.5), guide = FALSE) +
#  #stat_density2d(geom="tile", aes(fill = ..density..), contour = FALSE) +
#  #geom_smooth() +
#  #scale_colour_gradientn(colours = topo.colors(10)) +
#  coord_cartesian(xlim = c(0, 1)) +
#  #scale_y_log10() +
#  ggtitle("Neighbour Correlation Analysis") +
#  theme(plot.title = element_text(size=18,vjust=1.1, face="bold")) +
#  ylab("Neighbour correlation Standard Deviation") +
#  xlab("Neighbour correlation Mean") +
#  theme(axis.title=element_text(size=14)) +
#  theme(legend.justification=c(0,0), legend.position=c(.8,0.65),
#        #legend.key.size = unit(0.6, "cm"),
#        legend.text = element_text(size = 14, colour = "black", angle = 0),
#        legend.title = element_text(size = 14),
#        legend.box = "horizontal") +
#  theme(legend.key = element_blank()) +
#  theme(legend.background = element_rect(fill=alpha('grey', 0.0)))

#ggsave(paste(NameOfRegion,"_SV_NeighbourCorrelationSelected.pdf",sep=""), width=10, height=6, dpi=700, useDingbats=FALSE)


#cat('Region_Profiler.R Done\n')

#******From this point you can use the cluster to get the cells of the region*******
