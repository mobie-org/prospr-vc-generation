#print stuff
print("Grouping the Supervoxels into cells in the cluster")

# Get arguments
# Usage example
# Rscript Clusters_Extraction_6dpf_npix3.R FileCreatedWithRegion_Profiler.R
# Rscript Clusters_Extraction_6dpf_npix3.R RegPro...SupervoxelsExpression_NoCorrRemoved

# helper functions is assumed to be in the scripts folder
HelperFunctionsFile = './scripts/Clusters_Extraction_functions.R'
#load functions
source(HelperFunctionsFile)
# get inputs
args = commandArgs(trailingOnly=TRUE)
# test if number of arguments are correct: if not, return an error
if (length(args)<1) {
  stop("Script called incorrectly.
       Please provide:
       Output expression file of Region_Profiler\n", call.=FALSE)
}
ExpFile = args[1]

library(ggplot2)
library(png)
library(rgl)

############################################
#########Define your directory
setwd(dirname(ExpFile))
#########get a matrix
inputfilename <- basename(ExpFile)
load(inputfilename)
#########get the supervoxels coordinate file
SVcoordinates=read.table('../SuperVoxels_npix3_SpatialInfo.txt',sep=",",header=T)
#########Change the name of the matrix and coordinates (I don't know now why this is for)
initMat <- Regionvoxels_reduced
initCoord <- SVcoordinates[as.character(rownames(initMat)),]
#########Define your size range for your clusters
#For 6dpf npix3:
clsize <- c(20,80)  #take on account the SV size, that each cell is covered by 250-1000 pixels, and the bilaterality
############################################

print("data loaded, ready to process")

#Specify a list, and iterate over the elements of the list until the subdivision is achieved
#Store the good clusters in another list
tempclusters <- list()
finalclusters <- list()
#Initiate it with all the elements in one vector
tempclusters[[1]] <- rownames(initMat)
names(tempclusters) <- ""

#names are arranged so the hierachy is preserved: e.g: ABBAAABAABA...
while(length(tempclusters)>0){
  #select the voxels to cluster
  vox2clust <- tempclusters[[1]]
  #get the name
  clustname <- names(tempclusters[1])
  #remove them from the list
  tempclusters[[1]] <- NULL
  #get the matrix that corresponds to them
  subMat <- initMat[vox2clust,]
  #divide them
  clust <- splitMat(subMat)
  #check for size
  for(i in c(1:length(clust))){
    #get name (split) of the cluster:
    subclustname <- paste(clustname,names(clust)[i],sep="")
    #if size fits, save it as final, else save it back in tempclusters

    #This decision here can be a little bit more elaborate to take on account both sizes, their correlation... etc:
    #A more elaborated method is implemented:

    #if size is bigger than maxsize, send it again to the queue
    if (length(clust[[i]])>clsize[2]){
      tempclusters[[length(tempclusters)+1]] <- clust[[i]]
      names(tempclusters)[length(tempclusters)] <- subclustname  #not + 1 because it's updated
    }
    #if size is smaller than minsize, save the cluster (it will probably be removed, but just to know)
    if (length(clust[[i]])<clsize[1]){
      finalclusters[[length(finalclusters)+1]] <- clust[[i]]
      names(finalclusters)[length(finalclusters)] <- subclustname #not + 1 because it's updated
    }

    #if size is in between, check the values of correlation, spatial coherence and size
    #then subdivide it and check the same values for the subclusters
    #if there is an improvement, save the two subclusters separately, otherwise save the initial subdivision
    if (length(clust[[i]])>=clsize[1] & length(clust[[i]])<=clsize[2]){
      IniClust_Corr <- MyCorrFunc(clust[[i]],initMat)
      IniClust_SpCoh <- SpatCoheFunc(clust[[i]],initCoord)
      IniClust_Size <- length(clust[[i]])
      #subdivide and measure
      sub_clust <- splitMat(initMat[clust[[i]],])
      SubCl_1_Corr <-  MyCorrFunc(sub_clust[[1]],initMat)
      SubCl_1_SpCoh <-  SpatCoheFunc(sub_clust[[1]],initCoord)
      SubCl_1_Size <- length(sub_clust[[1]])
      SubCl_2_Corr <-  MyCorrFunc(sub_clust[[2]],initMat)
      SubCl_2_SpCoh <-  SpatCoheFunc(sub_clust[[2]],initCoord)
      SubCl_2_Size <- length(sub_clust[[2]])
      #Correct the values in case the length is = 1
      if(SubCl_1_Size<2){
        SubCl_1_Corr <- 1
        SubCl_1_SpCoh <- 0
      }
      if(SubCl_2_Size<2){
        SubCl_2_Corr <- 1
        SubCl_2_SpCoh <- 0
      }
      #Compare
      #if both correlation and spatial coherence improve, save them separate
      if(SubCl_1_Corr > IniClust_Corr & SubCl_2_Corr > IniClust_Corr & SubCl_1_SpCoh < IniClust_SpCoh & SubCl_2_SpCoh < IniClust_SpCoh){
        finalclusters[[length(finalclusters)+1]] <- sub_clust[[1]]
        names(finalclusters)[length(finalclusters)] <- paste(subclustname,"_1",sep="") #not + 1 because it's updated
        finalclusters[[length(finalclusters)+1]] <- sub_clust[[2]]
        names(finalclusters)[length(finalclusters)] <- paste(subclustname,"_2",sep="")
      }else{
        finalclusters[[length(finalclusters)+1]] <- clust[[i]]
        names(finalclusters)[length(finalclusters)] <- subclustname #not + 1 because it's updated
      }
      print(subclustname)
    }
  }
}



print("processing done, saving the data")

regionname <- strsplit(inputfilename,'_')[[1]][2]

save(finalclusters, file = paste(regionname,'_VirtualCells_Raw',sep=''))

print("data saved")

#######
# EVERYTHING THAT FOLLOWS WILL BE COMMENTED OUT FOR NOW
# CAN BE UNCOMMENTED IF WE NEED PICS
#######
#print("Cleaning the data and processing it for curation\n")

##create a dataframe with the final clusters:
###Cluster identifier: ABBAB...
###Supervoxels IDs: 1256,1252,...
###Inner-cluster correlation value: 0.86
###Spatial coherence / bilaterality: 11
###Cluster size: 14
#clustersDF <- data.frame(
#  clusterID <- names(finalclusters),
#  supervoxelsID <- as.character(lapply(finalclusters,function(x){paste(x,collapse=",")})),
#  correlation <- as.numeric(lapply(finalclusters,function(x){MyCorrFunc(x,initMat)})),
#  spatialCoherence <- as.numeric(lapply(finalclusters,function(x){SpatCoheFunc(x,initCoord)})),
#  clusterSize <- as.numeric(lapply(finalclusters,length))
#)
#names(clustersDF) <- c("clustersID","supervoxelsID","correlation","spatialCoherence","clusterSize")

##########
#sizelimit <- 15
##plot nicely
#p <- ggplot(clustersDF[clustersDF$clusterSize>sizelimit,], aes(x = correlation, y = spatialCoherence, color=clusterSize)) +
#  geom_point(size = 3, alpha = 0.8) +
#  #scale_colour_gradientn(colours = topo.colors(5)) +
#  scale_color_gradient(low = 'red', high = 'green') +
#  ggtitle(paste(regionname," clusters characteristics",sep="")) +
#  theme(axis.title=element_text(size=14)) +
#  theme(plot.title = element_text(size=18,vjust=1.1, face="bold")) +
#  theme(axis.title=element_text(size=14)) +
#  theme(legend.justification=c(0,0), legend.position=c(.85,0.65),
#        legend.key.size = unit(1, "cm"),
#        legend.text = element_text(size = 14, colour = "black", angle = 0),
#        legend.title = element_text(size = 14),
#        legend.box = "horizontal") +
#  theme(legend.key = element_blank()) +
#  theme(legend.background = element_rect(fill=alpha('grey', 0.0)))
#print(p)
#ggsave(paste(regionname,"_Clusters_characteristics.pdf",sep=""), width=10, height=10, dpi=700, useDingbats=FALSE)

#p <- ggplot(clustersDF[clustersDF$clusterSize>sizelimit,], aes(x = clusterSize, y = spatialCoherence, color=correlation)) +
#  geom_point(size = 3, alpha = 0.5) +
#  #scale_colour_gradientn(colours = topo.colors(5)) +
#  scale_color_gradient(low = 'red', high = 'green') +
#  ggtitle(paste(regionname," clusters characteristics",sep="")) +
#  theme(axis.title=element_text(size=14)) +
#  theme(plot.title = element_text(size=18,vjust=1.1, face="bold")) +
#  theme(axis.title=element_text(size=14)) +
#  theme(legend.justification=c(0,0), legend.position=c(.85,0.65),
#        legend.key.size = unit(1, "cm"),
#        legend.text = element_text(size = 14, colour = "black", angle = 0),
#        legend.title = element_text(size = 14),
#        legend.box = "horizontal") +
#  theme(legend.key = element_blank()) +
#  theme(legend.background = element_rect(fill=alpha('grey', 0.0)))
#print(p)
#ggsave(paste(regionname,"_Clusters_characteristics_SizeVSspatial.pdf",sep=""), width=10, height=10, dpi=700, useDingbats=FALSE)


###Inspect All Cells Manually######
#cat("\n\nSaving images\n\n")
##Create 3D visualizations for every cell, and show where in the plot are they.
#dir.create("Cell_Inspection")
##load('ClustersHQ')
#setwd("Cell_Inspection")
##########
##clusters above size limit
#clustersAS <- clustersDF[clustersDF$clusterSize>sizelimit,]
#for(i in c(1:nrow(clustersAS))){

#  #see it in 3D
#  selected_cluster = clustersAS[i,"clustersID"]

#  cat(paste("Processing cluster",selected_cluster,"\n",sep = ' '))

#  #plot voxels
#  cols=rep("grey92",nrow(initCoord))
#  names(cols) = rownames(initCoord)
#  cols[strsplit(as.character(clustersAS[clustersAS$clustersID==selected_cluster,"supervoxelsID"]),',')[[1]]]='red'
#  plot3d(initCoord$Xpos,initCoord$Ypos,initCoord$Zpos,size=6,col='grey92',alpha=0.05)#,box=NA)
#  spheres3d(initCoord$Xpos[cols!='grey92'],initCoord$Ypos[cols!='grey92'],initCoord$Zpos[cols!='grey92'],radius=2,col=cols[cols!='grey92'],alpha=0.5)
#  aspect3d("iso")
#  par3d("windowRect"= c(0,0,800,800))
#  bg3d("black")
#  rgl.pop(type = "bboxdeco") #remove the box
#  rgl.viewpoint( theta = 0, phi = 0, fov = 60, zoom = 1)
#  rgl.snapshot(filename="SVviewTemp.png")


#  #Make the plot (run same parameters as above)
#  p <- ggplot(clustersAS, aes(x = correlation, y = spatialCoherence, color=clusterSize)) +
#    #p <- ggplot(clustersAS, aes(x = correlation*(clusterSize^size_correction), y = spatialCoherence/(clusterSize^size_correction), color=clusterSize)) +
#    #p <- ggplot(clustersAS, aes(x = correlation, y = spatialCoherence/(clusterSize^size_correction), color=clusterSize)) +
#    geom_point(size = 3, alpha = 0.8) +
#    #scale_colour_gradientn(colours = topo.colors(5)) +
#    scale_color_gradient(low = 'red', high = 'green') +
#    ggtitle(clustersAS[i,"clustersID"]) +
#    theme(axis.title=element_text(size=14)) +
#    theme(plot.title = element_text(size=18,vjust=1.1, face="bold")) +
#    theme(axis.title=element_text(size=14)) +
#    theme(legend.position="none") +
#    #   theme(legend.justification=c(0,0), legend.position=c(.85,0.65),
#    #         legend.key.size = unit(1, "cm"),
#    #         legend.text = element_text(size = 14, colour = "black", angle = 0),
#    #         legend.title = element_text(size = 14),
#    #         legend.box = "horizontal") +
#    #   theme(legend.key = element_blank()) +
#    #   theme(legend.background = element_rect(fill=alpha('grey', 0.0))) +
#    #  geom_abline(intercept = eq_int, slope = eq_sl) + #completely arbitrary shit
#    #  geom_point(aes(x = clustersAS[i,"correlation"]*(clustersAS[i,"clusterSize"]^size_correction), y = clustersAS[i,"spatialCoherence"]/(clustersAS[i,"clusterSize"]^size_correction)),shape = 21, colour = "black", size = 4, stroke = 2)
#    geom_point(aes(x = clustersAS[i,"correlation"], y = clustersAS[i,"spatialCoherence"]),shape = 21, colour = "black", size = 4, stroke = 2)
#  #print(p)
#  ggsave("DotViewTemp.png",width=5, height = 5, dpi=150)

#  # create a subplot and save the two images
#  SVview = readPNG(source = "SVviewTemp.png")
#  DotView = readPNG(source = "DotViewTemp.png")
#  # setup plot
#  #plot.new()
#  dev.new()
#  par(mar=rep(0,4)) # no margins
#  # layout the plots into a matrix w/ 2 columns, by row
#  layout(matrix(1:2, ncol=2, byrow=TRUE))
#  # do the plotting
#  plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
#  rasterImage(SVview,0,0,1,1)
#  plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
#  rasterImage(DotView,0,0,1,1)

#  # write to PDF
#  dev.print(png, paste(regionname,"_Clusters_",clustersAS[i,"clustersID"],".png",sep=""), width=2000,height=1000)
#  dev.off()

#  # delete temporary images
#  file.remove("SVviewTemp.png")
#  file.remove("DotViewTemp.png")


#}



#cat('Images saved\n')

###########
##Create a dicctionary to Sort the cells by their spatial coherence value
#txtout <- list()
#for(i in c(1:nrow(clustersAS))){
#  txtout[[i]] <- (paste(clustersAS[i,"clustersID"],sprintf("%05.1f",clustersAS[i,"spatialCoherence"]), sep=","))
#}

#lapply(txtout, write, "000SortedCells.txt", append=TRUE, ncolumns=2000)

########
#cat("Clusters Extraction Done\n")
