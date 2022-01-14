library(here)
library(tidyverse)
library(CARNIVAL)

# capture arguments
args <- commandArgs(trailingOnly = TRUE)

# set wd 
setwd(here())

# list index
index <- as.numeric(args[1])

# read input list
inputList <- readRDS(here("data/carnival_input_list.rds"))

# move to tempDir
tempDir <- paste0( here( "tmp/", index ) )
dir.create(path = tempDir, recursive = TRUE, showWarnings = FALSE)
setwd(tempDir)

# execute carnival
carnivalRes <- runCARNIVAL(inputObj = inputList[[index]]$inputs,
                           measObj = inputList[[index]]$measurements, 
                           weightObj = inputList[[index]]$weights,
                           netObj = inputList[[index]]$pkn,
                           solver = "cplex",
                           solverPath = here("raw/cplex"),
                           timelimit = 3600, 
                           threads = 2)

# create results dir
if( ! dir.exists(here("data/carnival_results")) ) dir.create(here("data/carnival_results"))

# write results
saveRDS(carnivalRes, here("data/carnival_results/", paste0(index, ".rds" )))

