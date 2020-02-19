# variables read from command line.

# ProSPr_6dpf_SuperVoxelPixCount.R
# Curated_MEDs directory
# TrackEM segmentation directory
# number of pixels to use for the supervoxels
# and output directory

# Usage example
# Rscript ProSPr_6dpf_SuperVoxelPixCount.R 4Curated_MEDs_Good/ TrackEM/ 3 .


####### Check the variables#####

# #use the dapi shell to generate the voxel space and to restrict the number of voxels in the calculations
# #the interior of the animal is 255 and the outside is 0.
# this is hardcoded for now
library(tiff)
dapishell = 'helper_files/DapiShellInv.tif'
ds = readTIFF(dapishell,all=T,as.is=F)

# get inputs
args = commandArgs(trailingOnly=TRUE)
# test if number of arguments are correct: if not, return an error
if (length(args)<4) {
  stop("Script not called correctly.
       Please provide Curated_MEDs directory,
       TrackEM segmentation directory,
       number of pixels to use for the supervoxels,
       and output directory\n", call.=FALSE)
}
bspmdir = args[1]
trackemdir = args[2]
npix = strtoi(args[3])
outputFolder = args[4]

setwd(outputFolder)

dir.create(paste(outputFolder,'/npix',npix,sep = ''))

############# Generate Voxel Space ######

zs = length(ds)
xs = ncol(ds[1][[1]])
ys = nrow(ds[1][[1]])

#Generate voxel space
startpoint = ceiling(npix/2)
zcentlist = seq(startpoint,zs,npix)
xcentlist = seq(startpoint,xs,npix)
ycentlist = seq(startpoint,ys,npix)
totnsupvox = length(zcentlist)*length(xcentlist)*length(ycentlist) #total number of supervoxels
svmatrix = matrix(NA,totnsupvox,3)
#supervoxels IDs
rownames(svmatrix) <- c(1:totnsupvox)
count = 0
#Create a neighbouring 3D matrix at the same time, storing the ID of each SV in a 3D matrix
neighbour_list <- as.list(rep(NA,length(zcentlist))) #as z-stacks
names(neighbour_list) <- zcentlist
rep_mat <- matrix(NA,length(ycentlist),length(xcentlist)) #as x-y images
rownames(rep_mat) <- ycentlist
colnames(rep_mat) <- xcentlist
for(i in c(1:length(neighbour_list))){
  neighbour_list[[i]] <- rep_mat
}

#fill the matrix and list
for(r in zcentlist){
  for(s in xcentlist){
    for(t in ycentlist){
      count=count+1
      svmatrix[count,] = c(t,s,r) #same order as before, as it has not changed.
      neighbour_list[[as.character(r)]][as.character(t),as.character(s)] <- as.character(count)
    }
  }
}

# Restrict supervoxels to those with dapi signal
#Measure intensity in supervoxel
measmat = matrix(0,totnsupvox,1)
for(j in 1:totnsupvox){
  #y,x,z
  ycent = svmatrix[j,1]
  xcent = svmatrix[j,2]
  zcent = svmatrix[j,3]

  zstart = zcent-floor(npix/2)
  zstop = zcent+floor(npix/2)
  ystart = ycent-floor(npix/2)
  ystop = ycent+floor(npix/2)
  xstart = xcent-floor(npix/2)
  xstop = xcent+floor(npix/2)
  if(zstop>zs){zstop=zs}
  if(ystop>ys){ystop=ys}
  if(xstop>xs){xstop=xs}

  tosum = 0
  for(i in zstart:zstop){
    tosum = tosum + sum(ds[i][[1]][ystart:ystop,xstart:xstop])
  }
  measmat[j,1] = tosum
}
svmatrix = cbind(svmatrix,measmat)

svmatrix = svmatrix[svmatrix[,4]>0,][,1:3]
colnames(svmatrix) = c("Ypos","Xpos","Zpos")

cat('SuperVoxel space generated\n')

############# Get neighbours ########
#Create a table with the neighbours of each SV
#function to use apply, which makes it waaaaaay faster
getNeighbours <- function(SVposition){
  #get the spatial position
  SVy <- SVposition["Ypos"]
  SVx <- SVposition["Xpos"]
  SVz <- SVposition["Zpos"]

  #get the possible search spaces for each dimension
  #default cases (the own supervoxel)
  ZSS <- c(SVz)
  YSS <- c(as.character(SVy))
  XSS <- c(as.character(SVx))
  #add the previous position if they are not in the edge
  if (SVz != zcentlist[1]){ZSS <- append(ZSS,SVz-npix)}
  if (SVz != tail(zcentlist,1)){ZSS <- append(ZSS,SVz+npix)}
  if (SVy != ycentlist[1]){YSS <- append(YSS,as.character(SVy-npix))}
  if (SVy != tail(ycentlist,1)){YSS <- append(YSS,as.character(SVy+npix))}
  if (SVx != xcentlist[1]){XSS <- append(XSS,as.character(SVx-npix))}
  if (SVx != tail(xcentlist,1)){XSS <- append(XSS,as.character(SVx+npix))}

  SVNeigh_vec <- vector()
  #for each possible neighbour value in Z
  for (z in ZSS){
    #like this should be faster:
    SVNeigh_vec <- append(SVNeigh_vec, as.character(neighbour_list[[as.character(z)]][YSS,XSS]))
  }
  #add NAs to the end of the vector if it is not complete (because the SV being in the edge)
   if(length(SVNeigh_vec)<27){
     SVNeigh_vec <- append(SVNeigh_vec,rep(NA,27-length(SVNeigh_vec)))
   }

  SVNeigh_vec
}
neighbour_mat <- t(apply(svmatrix,1,getNeighbours))
#column 1 is the own SV
neighbour_mat <- neighbour_mat[,c(2:ncol(neighbour_mat))]

save(neighbour_mat,file = paste("./npix",npix,"/SuperVoxels_npix",npix,"_neighbour_matrix",sep=""))

cat('Neighbour matrix done and saved\n')

############# Measure genes ########
bspmfiles = vector()
#get the tif files for the genes
bspmfolders = list.files(path = bspmdir)
for (gene in bspmfolders){
  MEDmap = list.files(path = paste(bspmdir,gene,sep=''),pattern = "*_MEDs.tif")
  if (length(MEDmap)==1){
    bspmfiles = append(bspmfiles,paste(bspmdir,gene,MEDmap,sep='/'))
  }
}
#get the tif files for the TrackEM segmented animal regions
bspmfolders = list.files(path = trackemdir)
for (gene in bspmfolders){
  MEDmap = list.files(path = paste(trackemdir,gene,sep=''),pattern = "*_MEDs.tif")
  if (length(MEDmap)==1){
    bspmfiles = append(bspmfiles,paste(trackemdir,gene,MEDmap,sep='/'))
  }
}


#for loop for genename
persigarr = vector()
for (impath in bspmfiles){
  # read image
  bspm = readTIFF(impath,all=T,as.is=T) #not-binarized MED

  #Binarize the MEDs, as their pixel value is the MED index
  for(k in 1:length(bspm)){
    bspm[k][[1]] <- ifelse(bspm[k][[1]]>0,1,bspm[k][[1]])
  }

  #Measure intensity in supervoxel
  measmat = matrix(0,nrow(svmatrix),1)
  #extract gene name:
  colnames(measmat) = tail(strsplit(impath,'/')[[1]],n=2)[1]

  for(j in 1:length(measmat)){
    #y,x,z
    #restrict this search to the voxels occupied by the dapi
    ycent = svmatrix[j,1]
    xcent = svmatrix[j,2]
    zcent = svmatrix[j,3]

    zstart = zcent-floor(npix/2)
    zstop = zcent+floor(npix/2)
    ystart = ycent-floor(npix/2)
    ystop = ycent+floor(npix/2)
    xstart = xcent-floor(npix/2)
    xstop = xcent+floor(npix/2)
    if(zstop>zs){zstop=zs}
    if(ystop>ys){ystop=ys}
    if(xstop>xs){xstop=xs}

    tosum = 0
    for(k in zstart:zstop){
      tosum = tosum + sum(bspm[k][[1]][ystart:ystop,xstart:xstop])
    }
    measmat[j,1] = tosum
  }

  svmatrix = cbind(svmatrix,measmat)

  cat(paste(colnames(measmat)," Done\n"))
}

##############

#normalize it from 0 to 1
svmatrix[, 4:ncol(svmatrix)] = svmatrix[, 4:ncol(svmatrix)]/(npix^3)

###############

#Save matrix into separate files for the spatial and the expression information
write.table(svmatrix[, 1:3],file=paste("./npix",npix,"/SuperVoxels_npix",npix,"_SpatialInfo.txt",sep=''),quote=F,sep=',',row.names=T)
write.table(svmatrix[, 4:ncol(svmatrix)],file=paste("./npix",npix,"/SuperVoxels_npix",npix,"_GeneExpression.txt",sep=''),quote=F,sep=',',row.names=T)

cat('Matrix saved, I am done\n')
