#!/bin/bash

# For now I assume it's run on cluster
# We need this specific version of R, might not work with later ones
module load R/3.4.0-foss-2016b

# The pathes to MEDs and TrackEM hardcoded for now
# TODO: add as an optional argument
# Number of pixels (supervoxel length) will be hardcoded to 3
# since further processing heavily relies on this
Rscript ./scripts/ProSPr_6dpf_SuperVoxelPixCount.R /g/arendt/PrImR/ProSPr6/4SPM_Binarization/CuratedMEDs_Good/ /g/arendt/PrImR/ProSPr6/TrackEM/ 3 $1

# The previous script outputs everything into a npix3 folder
npix_folder="${1}/npix3/"

# The list of regions from TrackEM, let it be hardcoded
declare -a ListOfRegions=("CrypticSegment" "Head Lmo4" "Head -Lmo4" "PNS"
                          "Pygidium" "RestOfAnimal" "Stomodeum" "VNC" "All MHCL4")

# For each region take a separate process, then wait till they all finish
for region in "${ListOfRegions[@]}"; do
  Rscript ./scripts/Region_Profiler.R $npix_folder $region &
done

wait

# The same here
for file in $(find $npix_folder -name '*SupervoxelsExpression_NoCorrRemoved'); do
  Rscript ./scripts/Clusters_Extraction_6dpf_npix3.R $file &
done

wait

Rscript ./scripts/PullClustersFromRegions_ALLCLUSTERS.R $npix_folder

Rscript ./scripts/CombineVirtualCells_ALLCLUSTERS.R $npix_folder

Rscript ./scripts/ExtractCellularModelsForFiji.R $npix_folder

# TODO: environment for it
./scripts/vc_coord_to_volume.py ."${npix_folder}CellModels_ALL_coordinates.tsv" "${npix_folder}volume_prospr_space_clust"

Rscript ./scripts/scripts/GetExpressionForVirtualCellsCurated.R $npix_folder
