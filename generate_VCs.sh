#!/bin/bash

# If pathes to MEDs and body parts segmentation folder are not provided
# take the default ones

if [ -z "$2" ]
then
        MED_DIR="/g/arendt/PrImR/ProSPr6/4SPM_Binarization/CuratedMEDs_Good/"
else
        MED_DIR=$2
fi

if [ -z "$3" ]
then
        SEGM_DIR="/g/arendt/PrImR/ProSPr6/TrackEM/"
else
        SEGM_DIR=$3
fi


# Number of pixels (supervoxel length) will be hardcoded to 3
# since further processing heavily relies on this
Rscript ./scripts/ProSPr_6dpf_SuperVoxelPixCount.R $MED_DIR $SEGM_DIR 3 $1

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

#ImageJ might crash here, but will probably do the job properly anyway
python ./scripts/vc_coord_to_volume.py "${npix_folder}CellModels_ALL_coordinates.tsv" "${npix_folder}volume_prospr_space_clust"

Rscript ./scripts/GetExpressionForVirtualCellsCurated.R $npix_folder

# add an empty cell in the first row, otherwise the gene names are shifted by one column left
sed -i -e 1's/.*/\t&/' "${npix_folder}CellModels_ALL_profile_clust_curated.tsv"

# if the name dict file isn't present, download it
if [ ! -f "helper_files/new_name_lut.json" ]; then
    wget -P helper_files/ https://github.com/platybrowser/platybrowser-backend/raw/master/misc/new_name_lut.json
fi

python ./scripts/update_table_gene_names.py "${npix_folder}CellModels_ALL_profile_clust_curated.tsv" helper_files/new_name_lut.json
