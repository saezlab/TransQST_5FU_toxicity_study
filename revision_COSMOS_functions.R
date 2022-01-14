download_omnipath <- function(url){

  return(read.table(url, sep = '\t', header = TRUE))
}

format_dorothea <- function(url = 'http://omnipathdb.org/interactions?datasets=tfregulons&tfregulons_levels=A,B,C&genesymbols=1&fields=sources,tfregulons_level')
{

  dorothea <- download_omnipath(url)
  dorothea <- dorothea[,c(4,3,6,7)]
  dorothea$sign <- dorothea$is_stimulation - dorothea$is_inhibition
  dorothea <- dorothea[dorothea$sign != 0,]
  dorothea <- dorothea[,c(1,2,5)]

  ensembl = biomaRt::useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")

  G_list <- biomaRt::getBM(filters = "hgnc_symbol",
                  attributes = c('hgnc_symbol','entrezgene_id', "description"),
                  values = unique(dorothea$source_genesymbol), mart = ensembl)

  gene_mapping <- as.character(G_list[,2])
  names(gene_mapping) <- G_list[,1]

  for(i in 1:length(dorothea[,1]))
  {
    dorothea[i,2] <- gene_mapping[dorothea[i,2]]
  }

  dorothea$source_genesymbol <- paste0("X",dorothea$source_genesymbol)

  return(dorothea)
}



####################################################################################################
####################################################################################################
####################################################################################################

downstream_neighbours <- function(meta_network, n_steps, input_names)
{
  meta_g <- igraph::graph_from_data_frame(meta_network[,c(1,3,2)])

  input_names = input_names[input_names %in% igraph::V(meta_g)$name]

  dn_nbours <- igraph::ego(graph = meta_g, order = n_steps, nodes = input_names, mode = "out")

  sub_nodes <- unique(names(unlist(dn_nbours)))

  meta_network <- meta_network[meta_network$source %in% sub_nodes & meta_network$target %in% sub_nodes,]

  return(meta_network)
}

filter_inputs <- function(inputs, meta_network)
{
  meta_g <- igraph::graph_from_data_frame(meta_network[,c(1,3,2)])

  inputs <- inputs[,names(inputs) %in% igraph::V(meta_g)$name]

  return(inputs)
}

####################################################################################################
####################################################################################################
####################################################################################################

filter_TF_sign <- function(meta_network, ttop_RNA, inputs, TF_targets, my_stat = "t", my_threshold = 1)
{
  meta_network$target_sign <- sapply(meta_network$target,function(x,ttop_RNA)
  {
    if(x %in% ttop_RNA$ID)
    {
      return(as.numeric(ttop_RNA[ttop_RNA$ID == x,my_stat]))
    } else
    {
      return(NA)
    }
  },ttop_RNA = ttop_RNA)

  inputs <- inputs[,names(inputs) %in% TF_targets[,2]]

  meta_network$TF_sign <- sapply(meta_network$source, function(x, inputs)
  {
    if(x %in% names(inputs))
    {
      return(as.numeric(inputs[,x]))
    } else
    {
      return(NA)
    }
  },inputs = inputs)

  meta_network$target_sign <- ifelse(abs(meta_network$target_sign) > my_threshold, meta_network$target_sign, 0)

  meta_network <- meta_network[!(meta_network$source %in% TF_targets[,2] & meta_network$target_sign == 0),]

  meta_network <- meta_network[is.na(meta_network$TF_sign) | sign(meta_network$target_sign*meta_network$interaction) == sign(meta_network$TF_sign),]
  meta_network <- meta_network[,c(1,2,3)]

  return(meta_network)
}

get_TF_sign_from_CARNI <- function(carni_res, TF_targets)
{
  TF_signs <- as.data.frame(carni_res$nodesAttributes)
  TF_signs <- TF_signs[TF_signs$AvgAct != 0,]

  TF_signs <- TF_signs[TF_signs$Node %in% TF_targets[,2],]
  row.names(TF_signs) <- TF_signs$Node
  TF_signs <- TF_signs[,"AvgAct",drop = F]
  TF_signs <- as.data.frame(t(TF_signs))

  return(TF_signs)
}
