# habitat_resotration
This repository contains the code for spatially explicit metapopulation models used in the manuscript "Habitat restoration in spatially explicit metacommunity models".

Each file is named afer the interaciton type (e.g. pairwise_mutualism) and the restoration method (e.g. random) simulated by the model coded in that file.

The code in all files has the following structure:

1) HABITAT DESTRUCTION
   - input parameters for habitat destruction stage
   - function "dynamicsD" which computes the regional abundance as a function of habitat loss during habitat destruction
   
2) HABITAT RESTORATION
   - input parameters for habitat restoration stage
   - function "dynamicsR" which computes the regional abundance as a function of habitat loss during habitat restoration
   - code to replicate the habitat resotration simulation (using funciton "dynamiscR") n times
   - calculation of restoration efficiency
   - storage of results for post-processing
