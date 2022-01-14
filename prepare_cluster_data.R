# prepare the input list for COSMOS on the cluster

# Run Carnival on a single sample

library(cosmos)
library(tidyverse)
library(tidyr)


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
meta_network_colon <- read_rds("./data/carnival_inputs/meta_pkn_colon.rds")
meta_network_si <-	read_rds("./data/carnival_inputs/meta_pkn_si.rds")
meta_pkn <- readr::read_csv(file = "./data/meta_prior_knowledge/meta_network_carnival_ready_exch_solved_fullomni_metfiltered.csv")
metabolomics <- read_rds("./data/carnival_inputs/human_metabolomics_data.rds")
tf_activity <- read_rds("./data/carnival_inputs/human_tf_data.rds")
rna_seq_colon <- read_rds("./data/5-FU in vitro human/transcriptomics/colon/colon_transcriptomics_diffExp_filtered_data.rds")
rna_seq_si <- read_rds("./data/5-FU in vitro human/transcriptomics/SI/SI_transcriptomics_diffExp_filtered_data.rds")
expressed_genes_SI <- read_rds("./data/carnival_inputs/expressed_genes_human_SI.rds")
expressed_genes_colon <- read_rds("./data/carnival_inputs/expressed_genes_human_colon.rds")



# prepare meta network -------------------
# we keep the folllowing samples only, because they have significant metabolomics:
keep_sample  = c("colon_1000uM_24h", "colon_1000uM_48h",
				 "colon_1000uM_72h", "colon_100uM_24h",
				 "si_1000uM_24h",    "si_1000uM_48h",
				 "si_1000uM_72h"   )

# nothing to do:
readr::write_rds(meta_pkn,"./data/carnival_inputs_cluster/meta_network.rds")

# prepare tf_regulon -------------
tf_regulon = cosmos::load_tf_regulon_dorothea()
tf_regulon <- tf_regulon %>% select(tf,sign,target)
write_rds(tf_regulon,"./data/carnival_inputs_cluster/tf_regulon.rds")


# prepare signaling/TF data -----------------
# TF activity is already filtered for Abs(activity) > 1.96 and TF in AB confidence

tf_data <- tf_activity %>% select(sample,activity,cosmosID) %>%
	dplyr::mutate(sample = gsub("SI","si",sample)) %>%
	dplyr::mutate(activity = sign(activity)) %>%
	group_by(sample) %>%
	arrange(sample) %>%
	group_split() %>%
	setNames(.,unlist(lapply(.,function(x)x$sample[[1]]))) %>%
	lapply(.,function(df){
		signaling_data = df$activity
		names(signaling_data) = df$cosmosID
		return(signaling_data)
	})

# the elements of the list are in order:
stopifnot(all(order(names(tf_data)) == sort(order(names(tf_data)))))


tf_data <- tf_data[names(tf_data) %in% keep_sample]

write_rds(tf_data,"./data/carnival_inputs_cluster/tf_data.rds")



# prepare metabolomics data -----------------
# generate the same sampleID as above and filter noise.
# use pval = .05
# there are different ions that have the same name. Why? no idea...
#  -- fold change is consistent?t
#
# we derive a score: -log10(qvalue)*sign(FC)
# this score ranges [-20 to 20], but COSMOS cares about data in [-1 0 1]
# so we take the sign to convert to {-1,1}
#
# filter for metabolites in compartment e
metabolomics_data <- metabolomics %>%
	filter(grepl(pattern = "5-FU",sample_id)) %>%
	# format sampleID
	separate(col = sample_id,into = c("organ","time",NA,"condition"),sep = "_") %>%
	mutate(condition = stringr::str_extract(condition, "[0-1]+ uM")) %>%
	mutate(condition = gsub(" ", "",condition)) %>%
	mutate(time = gsub(" ","",time)) %>%
	mutate(sample = paste(organ,condition,time,sep="_")) %>%
	dplyr::select(sample,cosmosID,`log2 fc`,`q-value (FDR)`) %>%
	# filter for significance
	filter(`q-value (FDR)` < 0.05) %>%
	# calculate score
	mutate(score = -log10(`q-value (FDR)`)*sign(`log2 fc`)) %>%
	mutate(binarised_score = sign(score)) %>%
	# remove duplicated metabolites with inconsistent sign
	# decide if in sample all metabolite with same ID has same score:
	group_by(sample,cosmosID) %>%
	nest() %>%
	mutate(consistency_data = purrr::map(data,function(df){
		df %>% mutate(consistent_data = (all(binarised_score==1) | all(binarised_score==-1)))
	})) %>% select(-data) %>% unnest(consistency_data) %>%
	# filter out inconsistent metabolites
	dplyr::filter(consistent_data) %>%
	group_by(sample,cosmosID) %>%
	summarise(binarised_score  = mean(binarised_score)) %>%
	# format to list
	select(sample,binarised_score,cosmosID) %>%
	group_by(sample) %>%
	mutate(n_metabolites = n()) %>%
	filter(n_metabolites > 5) %>%
	select(-n_metabolites) %>%
	arrange(sample) %>%
	group_split() %>%
	setNames(.,unlist(lapply(.,function(x)x$sample[[1]]))) %>%
	lapply(.,function(df){
		mtb_data = df$binarised_score
		names(mtb_data) = df$cosmosID
		return(mtb_data)
	})

stopifnot(all(order(names(metabolomics_data)) == sort(order(names(metabolomics_data)))))
write_rds(metabolomics_data,"./data/carnival_inputs_cluster/metabolomics_data.rds")


# diff gene expressed --------------------------------
# merge the 2 organ,
# format sampleID
# use the stat as statistics
diff_expr_genes <- bind_rows(rna_seq_colon %>% mutate(organ = "colon"),
		  rna_seq_si %>% mutate(organ = "si")) %>%
	mutate(sample = paste(organ,fileID,sep = "_")) %>%
	mutate(ID = convert_ensembl_to_entrezid(ensembl_gene_id)) %>%
	mutate(ID = paste0("X",ID)) %>%
	select(sample,ID,stat) %>%
	# remove duplicated:
	group_by(sample,ID) %>%
	summarise(stat = mean(stat)) %>%
	ungroup() %>%
	# format to list
	group_by(sample) %>%
	arrange(sample) %>%
	group_split() %>%
	setNames(.,unlist(lapply(.,function(x)x$sample[[1]]))) %>%
	lapply(.,function(df){
		mtb_data = df$stat
		names(mtb_data) = df$ID
		return(mtb_data)
	})


stopifnot(all(order(names(diff_expr_genes)) == sort(order(names(diff_expr_genes)))))

diff_expr_genes <- diff_expr_genes[names(diff_expr_genes) %in% keep_sample]

write_rds(diff_expr_genes,"./data/carnival_inputs_cluster/diff_expr_genes.rds")

### expressed genes ----------------------
# expressed genes will be the same across treatments of the organ

expressed_colon_gene_names <- expressed_genes_colon %>%
	filter(keep == TRUE) %>%
	mutate(entrezID = paste0("X",entrezID)) %>%
	pull(entrezID)

expressed_si_gene_names <- expressed_genes_SI %>%
	filter(keep == TRUE) %>%
	mutate(entrezID = paste0("X",entrezID)) %>%
	pull(entrezID)

sort(keep_sample)
expressed_genes_list <- vector("list",length(keep_sample))
names(expressed_genes_list) <- keep_sample

expressed_genes_list[1:4] = rep(list(expressed_colon_gene_names),4)
expressed_genes_list[5:7] =  rep(list(expressed_si_gene_names),3)



write_rds(expressed_genes_list,"./data/carnival_inputs_cluster/expressed_genes.rds")






