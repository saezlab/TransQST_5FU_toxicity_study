library(readr)
library(biomaRt)
library(piano)

setwd("~/Dropbox/COSMOS/")
source("scripts/revision_COSMOS_functions.R")
source("scripts/revision_piano_function.R")

#Load inputs
meta_network <- as.data.frame(
  read_csv("support/meta_network_carnival_ready_exch_solved_fullomni_metfiltered_expfiltered.csv"))

prots <- unique(c(meta_network$source,meta_network$target))
prots <- prots[!grepl("XMetab",prots)]
prots <- gsub("^X","",prots)
prots <- gsub("Gene[0-9]+__","",prots)
prots <- gsub("_reverse","",prots)
prots <- gsub("EXCHANGE.*","",prots)
prots <- unique(prots)
prots <- unlist(sapply(prots, function(x){unlist(strsplit(x,split = "_"))}))

gene_mapping <- "ensembl"

if(gene_mapping == "ensembl")
{
  ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  
  G_list <- getBM(filters = "entrezgene_id",
                  attributes = c('hgnc_symbol','entrezgene_id', "description"),
                  values = prots, mart = ensembl)
  
} else
{
  G_list <- as.data.frame(read_csv(gene_mapping))
}

gene_mapping <- G_list[,1]
names(gene_mapping) <- G_list[,2]

background <- sapply(prots, function(x, translation_dictionary){
  return(translation_dictionary[x])
},translation_dictionary = gene_mapping)
names(background) <- 1:length(background)

background <- unique(background)

hallmarks <- gmt_to_df("support/h.all.v7.1.symbols.gmt_20200709.gmt")

cosmos_att <- as.data.frame(
  read_csv("results_revisions/COSMOS_cluster_runs/run_fullomni/met_exp_filtered_continuous_input_sigfilter/subnet_combined_att_full_newDoro_long.csv"))

cosmos_att <- cosmos_att[cosmos_att$Nodes %in% background,]
cosmos_att <- unique(cosmos_att$Nodes)


hallmarks_ora <- runGSAhyper(genes = cosmos_att, universe = background, gsc = loadGSC(hallmarks))
hallmarks_ora <- as.data.frame(hallmarks_ora$resTab)
hallmarks_ora$pathway <- row.names(hallmarks_ora)

cgp <- gmt_to_df("support/c2.cgp.v7.1.symbols.gmt")
cgp <- cgp[cgp$gene %in% background,]

cgp_ora <- runGSAhyper(genes = cosmos_att, universe = background, gsc = loadGSC(cgp))
cgp_ora <- as.data.frame(cgp_ora$resTab)
cgp_ora$pathway <- row.names(cgp_ora)

cp <- gmt_to_df("support/c2.cp.v7.1.symbols.gmt")
cp <- cp[cp$gene %in% background,]

cp_ora <- runGSAhyper(genes = cosmos_att, universe = background, gsc = loadGSC(cp))
cp_ora <- as.data.frame(cp_ora$resTab)
cp_ora$pathway <- row.names(cp_ora)
