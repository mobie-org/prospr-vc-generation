# Generate Virtual Cells from ProSPr

## Requirements:

- Curated ProSPr output (a.k.a. MEDs)
- Manual segmentation of the animal into regions
- Environment as described in environment.txt (TODO)

## Usage:

Due to R peculiarities the bash script has to be run from the directory it is located in as follows:

`./generate_VCs.sh /PATH/TO/OUTPUT_DIR /PATH/TO/MEDS_DIR /PATH/TO/SEGM_DIR`

For example:

`./generate_VCs.sh /g/arendt/EM_6dpf_segmentation/GenerationOfVirtualCells/ProSPr_VirtualCells/ /g/arendt/PrImR/ProSPr6/4SPM_Binarization/CuratedMEDs_Good/ /g/arendt/PrImR/ProSPr6/TrackEM/`

Note: please, make sure the directory names __always__ end with a slash `/` since otherwise R doesn't parse them as directories and simply ignores them.
  
