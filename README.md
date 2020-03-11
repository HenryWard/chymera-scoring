# chymera-scoring

This repository contains all scripts used to score data from the paper "Genetic interaction mapping and exon-resolution functional genomics with a hybrid Cas9–Cas12a platform", 
available to read at *insert link here*. 

To run the script "score_paralogs.R" and generate scores for the paralog CHyMErA library, perform the following steps: 

1. Make sure that you have R version > 3.6.3 installed, along with the packages *ggplot2* and *ggthemes*.
2. Change the four parameters at the top of the script to match your local machine. 
    1. `setwd` takes the path of the *chymera-scoring* folder. 
	2. `output_folder` is the path to a desired output folder. This will be created automatically if it does not exist.
	3. `which_cell_line` is one of "hap1", "rpe1", or "hap1_torin" for scoring each screen, respectively.
3. Run the script, e.g. with `source("score_paralogs.R")`, to score the type of data specified in step 2.3. 

To run the script "score_dual_targeting.R" and generate scores for the dual-targeting CHyMErA libraries, perform the following steps: 

1. Make sure that you have R version > 3.6.3 installed, along with the packages *ggplot2* and *ggthemes*.
2. Change the four parameters at the top of the script to match your local machine. 
    1. `setwd` takes the path of the *chymera-scoring* folder. 
	2. `output_folder` is the path to a desired output folder. This will be created automatically if it does not exist.
	3. `guide_type` is one of "single", "dual", "paralog_single" or "paralog_dual" for scoring single-targeting guides or dual-targeting guides
	    from the dual-targeting library ("single" and "dual") or the paralog library ("paralog_single" and "paralog_dual"). 
3. Run the script, e.g. with `source("score_dual_targeting.R")`, to score the type of data specified in step 2.3. 

To run the script "score_exon.R" and generate scores for the exon deletion CHyMErA library, perform the following steps: 

1. Make sure that you have R version > 3.6.3 installed, along with the packages *ggplot2*, *ggthemes*, *dplyr* and *scales*.
2. Change the two parameters at the top of the script to match your local machine. 
    1. `setwd` takes the path of the *chymera-scoring* folder. 
	2. `output_folder` is the path to a desired output folder. This will be created automatically if it does not exist.
3. Run the script, e.g. with `source("score_exon.R")`.