# Synthesized datasets
The datasets in this folder are slightly different from the datasets used in the paper, as they have been synthesized from the real data using the `synthpop` package for privacy reasons. The results may therefore slightly differ from the results in the paper. The datasets are provided as uncompressed R data structures (`.rds` files).

## Cam-CAN volume
Synthesized data frame with 647 rows (participants) and 70 columns (ROIs). Values indicate grey matter volume.

## Cam-CAN wm
White matter tractography. Synthesized data frame with 646 rows (participants) and 42 columns (ROIs). Values indicate fractional anisotropy.

## Cam-CAN fun
List of synthesized resting state functional connectivity matrices of 5 participants for 90 ROIs. These come from 261 observations per participant. Synthetisation was done by sampling from a Wishart distribution with the original matrix as Sigma parameter and 261 degrees of freedom.