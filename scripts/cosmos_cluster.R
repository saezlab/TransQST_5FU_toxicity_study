library(here)
library(dplyr)
library(tibble)

library(readr)
library(purrr)
library(cosmos)

# capture arguments
args <- commandArgs(trailingOnly = TRUE)

# set wd
setwd(here())

# list index
index <- as.numeric(args[1])
rep <- as.numeric(args[2])

# read input list
diff_expr_genes_list <- read_rds(here("data/carnival_inputs_cluster/diff_expr_genes.rds"))
metabolomics_data_list <- read_rds(here("data/carnival_inputs_cluster/metabolomics_data.rds"))
tf_data_list <- read_rds(here("data/carnival_inputs_cluster/tf_data.rds"))
tf_regulon <- read_rds(here("data/carnival_inputs_cluster/tf_regulon.rds"))
meta_network <- read_rds(here("data/carnival_inputs_cluster/meta_network.rds"))
expressed_genes <- read_rds(here("data/carnival_inputs_cluster/expressed_genes.rds"))


# move to tempDir
tempDir <- paste0( here( "tmp/", index,"_", rep ) )
dir.create(path = tempDir, recursive = TRUE, showWarnings = FALSE)
setwd(tempDir)


preproc_opts <- default_CARNIVAL_options()
preproc_opts$solverPath = here("cplex/cplex") ## <-- adapt for cluster!
preproc_opts$solver = "cplex"
preproc_opts$threads = 1
preproc_opts$timelimit = 3600
preproc_opts$mipGAP = 0.1


# repetitions:
# shuffle the inputs, measurements and network to induce difference in multiple runs
set.seed(132*rep)
tf_data_rep = tf_data_list[[index]]
tf_data_rep = tf_data_rep[sample(1:length(tf_data_rep))]

metabolomics_data_rep = metabolomics_data_list[[index]]
metabolomics_data_rep = metabolomics_data_rep[sample(1:length(metabolomics_data_rep))]

meta_network_rep = meta_network[sample(1:nrow(meta_network)),]

# preprocess COSMOS inputs
prep_inputs <- preprocess_COSMOS(meta_network = meta_network_rep,
				  tf_regulon = tf_regulon,
				  signaling_data = tf_data_rep,
				  metabolic_data = metabolomics_data_rep,
				  diff_expression_data = diff_expr_genes_list[[index]],
				  diff_exp_threshold = 0.5,
				  expressed_genes = expressed_genes[[index]],
				  maximum_network_depth = 7,
				  remove_unexpressed_nodes = TRUE,
				  filter_tf_gene_interaction_by_optimization = TRUE,
				  CARNIVAL_options = preproc_opts
				  )

# create results dir
if( ! dir.exists(here("data/cosmos_preprocess_results")) ) dir.create(here("data/cosmos_preprocess_results"))
# write results
write_rds(prep_inputs, here("data/cosmos_preprocess_results/", paste0(index,"_",rep, ".rds" )))

signaling_data <- prep_inputs$signaling_data_bin
meta_network <- prep_inputs$meta_network
diff_expression_data <- prep_inputs$diff_expression_data_bin
metabolic_data <- prep_inputs$metabolic_data

run_opts <- default_CARNIVAL_options()
run_opts$solverPath = here("cplex/cplex") ## <-- adapt for cluster!
run_opts$solver = "cplex"
run_opts$threads = 1
run_opts$timelimit = 3600*15
run_opts$mipGAP = 0.1

res_network = run_COSMOS_signaling_to_metabolism(meta_network = meta_network,
												 metabolic_data = metabolic_data,
												 signaling_data = signaling_data,
												 CARNIVAL_options = run_opts,
												 test_run = FALSE)


# create results dir
if( ! dir.exists(here("data/cosmos_results")) ) dir.create(here("data/cosmos_results"))

# write results
write_rds(res_network, here("data/cosmos_results/", paste0(index,"_",rep, ".rds" )))

