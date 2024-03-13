# 2bRAD-workshop
These instructions were developed as an introduction for undergraduates researching corals from St. Croix (STX), USVI. These materials were made for participants to reference when learning to process next-generation sequences (2bRAD specifically), then implement population genetics and genotype-environment association analyses.

Notes in TACC_Workshop.pdf are beginner lessons in using TACC (Texas Advanced Computing Center) and shell scripting. Scripts in STX.sh contain the 2bRAD pipeline for processing sequences in the STX coral project. If implementing these scripts in other projects, please change the project IDs, genomes, and emails referenced throughout.

The Data folder contains all environmental and sample data. Internal folders contain genotype data organized by species:\
AAGA: *Agaricia agaricites*\
MCAV: *Montastrea cavernosa*\
OFAV: *Orbicella faveolata*\
PAST: *Porites astreoides*\
PSTR: *Pseudodiploria strigosa*\
SSID: *Siderastrea siderea*

These species folders contains the Identity-by-State genetic distance matrix, sample list, and lineage assignments. These files can be input to the following R scripts:

STX_popstructure.R examines population structure and cryptic genetic lineages within each species.
STX_spatial.R explores the distribution of cryptic genetic lineages across the map.
STX_environment.R finds genotype-environment associations using gradient forest.

Additionally, environmental heterogeneity across the St. Croix seascape can be summarized into ecoregions using STX_envPCA.R
And lastly, cryptic lineages within each species can be combined into a community abundance matrix, and associations between cryptic communities and ecoregions can be explored in Ecoregions.R

