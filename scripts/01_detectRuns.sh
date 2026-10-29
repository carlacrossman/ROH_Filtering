#!/bin/bash

module load r/4.4.0

PREFIX=$(basename $1)

sed 's/\t/ /g' 00_data/"$PREFIX"_chr.ped > 00_data/"$PREFIX"_chr.peds

##################################################################################
# detect RUNs
# based on https://cran.r-project.org/web/packages/detectRUNS/readme/README.html
##################################################################################

R --vanilla --args "$PREFIX" << EOF 

if (!requireNamespace("detectRUNS", quietly = TRUE)) {
  install.packages("detectRUNS")
}
library(detectRUNS)

args = commandArgs(TRUE)
PREFIX = args[1]

#1) detectRUNS (sliding-windows method)

    genotypeFile <- paste0("00_data/", PREFIX,"_chr.peds")
    mapFile <- paste0("00_data/", PREFIX,"_chr.map")
    
    # calculating runs with sliding window approach
    runs <- slidingRuns <- slidingRUNS.run(
      genotypeFile = genotypeFile, 
      mapFile = mapFile, 
      windowSize = 30, 
      threshold = 0.05,
      minSNP = 20, 
      ROHet = FALSE, 
      maxOppWindow = 2, 
      maxMissWindow = 2,
      maxGap = 500000, 
      minLengthBps = 100000, 
      minDensity = 1/100, # SNP/kbps
      maxOppRun = NULL,
      maxMissRun = NULL)
    
    write.table(runs, file = paste0("02_detectRUNS/", PREFIX, ".detectRuns2.runs"), quote= F, row.names = F)

#2) summarise results

    summaryList <- summaryRuns(runs = runs, mapFile = mapFile, genotypeFile = genotypeFile, Class = 6, snpInRuns = TRUE)

    save(summaryList, file = paste0("02_detectRUNS/", PREFIX, ".detectRuns2.summary.Rdata"))


#3) plot results

    pdf(paste0("02_detectRUNS/", PREFIX, ".detectRuns2.pdf"))
    plot_Runs(runs = runs)
    dev.off()

EOF
