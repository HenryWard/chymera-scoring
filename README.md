# chymera-scoring

This repository contains all scripts used to score data from the paper "Genetic interaction mapping and exon-resolution functional genomics with a hybrid Cas9–Cas12a platform", 
available to read at *insert link here*. 

To run the script "score_paralogs.R" and generate scores for the paralog CHyMErA libraries, perform the following steps: 

1. Download the original screening data from https://crispr.ccbr.utoronto.ca/chymera/#. 
2. Make sure that you have R version > 3.6.3 installed, along with the packages *ggplot2* and *ggthemes*.
3. Change the four parameters at the top of the script to match your local machine. 
    1. `setwd` takes the path of the *chymera-scoring* folder. 
	2. `input_file` is the path to the screening data downloaded in step 1. 
	3. `output_folder` is the path to a desired output folder. This will be created automatically if it does not exist.
	4. `which_cell_line` is one of "hap1", "rpe1", or "hap1_torin" for scoring each screen, respectively.
4. Run the script, e.g. with `source("paralog_scoring.R")`, to score the type of data specified in step 3.4. 