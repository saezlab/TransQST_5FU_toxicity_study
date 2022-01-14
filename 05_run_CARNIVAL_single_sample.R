# Run Carnival on a single sample


library(readr)
library(igraph)
library(CARNIVAL)
library(biomaRt)

source("./revision_COSMOS_functions.R")
source("./cosmos.R")
library(tidyverse)
library(dorothea)

###
# Threshold used to find significantly modified genes:
# We will use the "stat" variable. If we plot the histogram for non-significant
# stat, we get the cutoff at +-4.65 for colon and 4.66 for SI.
# run:  rna_seq_si %>% filter(padj >0.05) %>% pull(stat) %>% max()
my_stat_RNA <- "stat"
my_threshold_RNA <- 4.66



# Load inputs
if(FALSE){
	dorothea <- format_dorothea()
	write_rds(dorothea,"./data/carnival_inputs/dorothea_tf_regulons_omnipath.rds")
}else{
	dorothea <- read_rds("./data/carnival_inputs/dorothea_tf_regulons_omnipath.rds")
}
meta_network_colon <-	read_rds("./data/carnival_inputs/meta_pkn_colon.rds")
meta_network_si <-	read_rds("./data/carnival_inputs/meta_pkn_si.rds")
metabolomics <- read_rds("./data/carnival_inputs/human_metabolomics_data.rds")
tf_activity <- read_rds("./data/carnival_inputs/human_tf_data.rds")
rna_seq_colon <- read_rds("./data/5-FU in vitro human/transcriptomics/colon/colon_transcriptomics_diffExp_filtered_data.rds")
rna_seq_si <- read_rds("./data/5-FU in vitro human/transcriptomics/SI/SI_transcriptomics_diffExp_filtered_data.rds")
expressed_genes_SI <- read_rds("./data/carnival_inputs/expressed_genes_human_SI.rds")
expressed_genes_colon <- read_rds("./data/carnival_inputs/expressed_genes_human_colon.rds")

# adjust log2FC naming
rna_seq_si <- rna_seq_si %>%
	mutate(entrezID = gene_ensembl_to_entrez(ensembl_gene_id)) %>%
	mutate(ID = paste0("X",entrezID))

rna_seq_colon <- rna_seq_colon %>%
	mutate(entrezID = gene_ensembl_to_entrez(ensembl_gene_id)) %>%
	mutate(ID = paste0("X",entrezID))

# entrez and ensembl IDs are different and does not map 1-1, because of e.g. splicing.
# The RNA data is used to constrain the PKN: remove nodes which did not respond.
# First, we remove nodes which are not expressed.
rna_seq_si <- rna_seq_si %>%
	left_join(dplyr::rename(expressed_genes_SI, ensembl_gene_id = "gene_ensembl_id"),by = c("entrezID","ensembl_gene_id"))

rna_seq_colon <- rna_seq_colon %>%
	left_join(dplyr::rename(expressed_genes_colon, ensembl_gene_id = "gene_ensembl_id"),by = c("entrezID","ensembl_gene_id"))


### Select case
# colon 1000uM, time == 72

tf_activity_inp = tf_activity %>%
	filter(organ=="colon", time==72, concentration==1000) %>%
	dplyr::select(cosmosID,activity) %>% deframe() %>%
	as.matrix() %>% t() %>% as.data.frame()


metabolomics_inp = metabolomics %>%
	filter(organ=="colon",sample_id=="colon_72 h_5-FU_1000 uM vs vehicle") %>%
	dplyr::select(cosmosID,`log2 fc`) %>% deframe()	%>%
	as.matrix() %>% t() %>% as.data.frame()

rna_seq_inp = rna_seq_colon %>% filter( time==72, concentration==1000) %>%
	filter(keep) %>%
	group_by(ID) %>%
	summarise(stat = mean(stat))



#Preprocess pkn
names(tf_activity_inp)[which(!names(tf_activity_inp) %in% c(meta_network_colon$source,meta_network_colon$target))]

tf_activity_inp <- filter_inputs(tf_activity_inp, meta_network_colon)

metabolomics_inp <- filter_inputs(metabolomics_inp, meta_network_colon)
metabolomics_inp = metabolomics_inp[,which(names(metabolomics_inp) %in% c(meta_network_colon$source,meta_network_colon$target))]

meta_network_colon_inp <- downstream_neighbours(meta_network_colon, 8, names(tf_activity_inp))

###############

meta_network_colon_inp <- filter_TF_sign(meta_network = meta_network_colon_inp,
							   ttop_RNA = rna_seq_inp,
							   inputs = tf_activity_inp,
							   TF_targets = dorothea,
							   my_stat = my_stat_RNA,
							   my_threshold = my_threshold_RNA)

meta_network_colon_inp <- meta_network_colon_inp[complete.cases(meta_network_colon_inp),]
###############

CARNIVAL_Result <- runCARNIVAL(inputObj = sign(tf_activity_inp),
							   measObj = metabolomics_inp,
							   netObj = meta_network_colon_inp,
							   solverPath = "/Applications/CPLEX_Studio128/cplex/bin/x86-64_osx/cplex",
							   solver = "cplex",
							   timelimit = 2400,
							   mipGAP = 0.2)

write_rds(CARNIVAL_Result,paste0("./carnival_temp_results/","carnival_step1_results_",format(Sys.time(), "%b_%d_%Hh_%Mm"),".rds"))

save.image("./carnival_temp_results/carni_forwardrun_res_fullomni_fullfilter.RData")

TF_signs <- get_TF_sign_from_CARNI(CARNIVAL_Result, dorothea)

meta_network_colon_inp <- filter_TF_sign(meta_network = meta_network_colon_inp,
							   ttop_RNA = rna_seq_inp,
							   inputs = TF_signs,
							   TF_targets = dorothea,
							   my_stat = my_stat_RNA,
							   my_threshold = my_threshold_RNA)

CARNIVAL_Result_rerun <- runCARNIVAL(inputObj = sign(tf_activity_inp),
									 measObj = metabolomics_inp,
									 netObj = meta_network_colon_inp,
									 solverPath = "/Applications/CPLEX_Studio128/cplex/bin/x86-64_osx/cplex",
									 solver = "cplex",
									 timelimit = 1200,
									 mipGAP = 0.2)

write_rds(CARNIVAL_Result,paste0("./carnival_temp_results/","carnival_step2_results_",format(Sys.time(), "%b_%d_%Hh_%Mm"),".rds"))

save.image("./carnival_temp_results/carni_forwardrun_res_fullomni_fullfilter.RData")

meta_network_colon <-	read_rds("./data/carnival_inputs/meta_pkn_colon.rds")

### FOR RELIEF COSMOS CHEM TO SIG

# input_network <- meta_network[meta_network$source %in% names(metabolomics),]
#
# meta_network$edgeID <- paste0(meta_network$source, meta_network$interaction, meta_network$target)
# input_network$edgeID <- paste0(input_network$source, input_network$interaction, input_network$target)
#
# meta_network <- meta_network[!(meta_network$edgeID %in% input_network$edgeID),]
#
# meta_network <- meta_network[,-4]
# input_network <- input_network[,-4]
#
# input_network_chem_to_sig <- input_network
# input_network_chem_to_sig$source <- gsub("___[a-z]____", "",input_network_chem_to_sig$source)
#
# ##For later
# input_network_conector <- input_network
# input_network_conector$target <- input_network_chem_to_sig$source
# input_network_conector$interaction <- 1
# input_network_conector <- unique(input_network_conector)
# ##
#
# input_network_chem_to_comp <- input_network
# input_network_chem_to_comp$target <- input_network_chem_to_comp$source
# input_network_chem_to_comp$source <- input_network_chem_to_sig$source
# input_network_chem_to_comp$interaction <- 1
# input_network_chem_to_comp <- unique(input_network_chem_to_comp)
#
# input_network_chem_to_sig <- unique(as.data.frame(rbind(input_network_chem_to_sig,input_network_chem_to_comp)))
#
# meta_network <- as.data.frame(rbind(input_network_chem_to_sig, meta_network))
#
# metabolomics_reduced <- as.data.frame(t(metabolomics))
# metabolomics_reduced$metab <- gsub("___[a-z]____", "",row.names(metabolomics_reduced))
# metabolomics_reduced <- unique(metabolomics_reduced)
# row.names(metabolomics_reduced) <- metabolomics_reduced$metab
# metabolomics_reduced <- metabolomics_reduced[,-2,drop = F]
#
# metabolomics_reduced <- as.data.frame(t(metabolomics_reduced))
#
# meta_network <- downstream_neighbours(meta_network, 6, names(metabolomics_reduced))

### FOR RELIEF COSMOS CHEM TO SIG (END)

meta_network <- downstream_neighbours(meta_network_colon, 6, names(metabolomics_inp))

######
######
######

meta_network <- filter_TF_sign(meta_network = meta_network,
							   ttop_RNA = rna_seq_inp,
							   inputs = tf_activity_inp,
							   TF_targets = dorothea,
							   my_stat = my_stat_RNA,
							   my_threshold = my_threshold_RNA)

meta_network <- filter_TF_sign(meta_network = meta_network,
							   ttop_RNA = rna_seq_inp,
							   inputs = TF_signs,
							   TF_targets = dorothea,
							   my_stat = my_stat_RNA,
							   my_threshold = my_threshold_RNA)

####
####
####

# CARNIVAL_Result_2 <- runCARNIVAL(inputObj = sign(metabolomics_reduced),
#                                  measObj = tf_activity,
#                                  netObj = meta_network,
#                                  solverPath = "~/Documents/cplex",
#                                  solver = "cplex",
#                                  timelimit = 7200,
#                                  mipGAP = 0.2)

meta_network <- meta_network[complete.cases(meta_network),]


CARNIVAL_Result_2 <- runCARNIVAL(inputObj = sign(metabolomics_inp),
								 measObj = tf_activity_inp,
								 netObj = meta_network,
								 solverPath = "/Applications/CPLEX_Studio128/cplex/bin/x86-64_osx/cplex",
								 solver = "cplex",
								 timelimit = 2400,
								 mipGAP = 0.2)

write_rds(CARNIVAL_Result_2,paste0("./carnival_temp_results/","carnival_step3_results_",format(Sys.time(), "%b_%d_%Hh_%Mm"),".rds"))

save.image("./carnival_temp_results/carni_doublerun_res_fullomni_met_expfiltered.RData")

#####

metabolite_map <- metabolomics %>% dplyr::select(`possible annotation based on accurate mass`,pubchemID) %>%
	dplyr::rename(name = `possible annotation based on accurate mass`) %>% unique()

cosmos_forward <- format_COSMOS_res(cosmos_res = CARNIVAL_Result_rerun,
									metab_mapping = metabolite_map,
									measured_nodes = c(names(tf_activity_inp),names(metabolomics_inp)),
									omnipath_ptm = read_rds("./data/carnival_inputs/omnipath_ptm.rds"))

cosmos_backward <- format_COSMOS_res(cosmos_res = CARNIVAL_Result_2,
									 metab_mapping = metabolite_map,
									 measured_nodes = c(names(tf_activity_inp),names(metabolomics_inp)),
									 omnipath_ptm = read_rds("./data/carnival_inputs/omnipath_ptm.rds"))


write_csv(cosmos_forward[[1]], "data/results/carnival_downstream/sif_first_colon_1000uM_72h.csv")
write_csv(cosmos_forward[[2]], "data/results/carnival_downstream/att_first_colon_1000uM_72h.csv")

write_csv(cosmos_backward[[1]], "data/results/carnival_downstream/sif_second_colon_1000uM_72h.csv")
write_csv(cosmos_backward[[2]], "data/results/carnival_downstream/att_second_colon_1000uM_72h.csv")

full_sif <- as.data.frame(unique(rbind(cosmos_forward[[1]],cosmos_backward[[1]])))
full_att <- as.data.frame(unique(rbind(cosmos_forward[[2]][,-6],cosmos_backward[[2]][,-6])))
full_att <- unique(full_att[,c(1,6,7,8)])
full_att[duplicated(full_att$Nodes),]

full_att$Nodes <- gsub(",","_",full_att$Nodes)
full_sif$Node1 <- gsub(",","_",full_sif$Node1)
full_sif$Node2 <- gsub(",","_",full_sif$Node2)
full_sif$Weight <- as.numeric(full_sif$Weight)

full_sif <- full_sif %>%
	group_by(Node1, Node2, Sign) %>%
	summarise(Weight= sum(Weight))

full_sif[full_sif$Weight > 100, "Weight"] <- 100

full_att <- full_att[full_att$Nodes %in% full_sif$Node1 | full_att$Nodes %in% full_sif$Node2,]

write_csv(full_sif, "data/results/carnival_downstream/subnet_combined_sif_full_colon_1000uM_72h.csv")
write_csv(full_att, "data/results/carnival_downstream/subnet_combined_att_full_colon_1000uM_72h.csv")

display_node_neighboorhood("ATF2",sif = full_sif, att = full_att, n = 3)
