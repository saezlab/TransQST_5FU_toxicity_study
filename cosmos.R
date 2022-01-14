

#' prepare_metabolomics_names
#'
#' takes the metabolites by pubchemID. Assignes compartments to metabolites
#' based on prior knowledge and filters for metabolites that are in the prior
#' knowledge network.
#'
#' @param metabolomics_pubchem pubChemID of metabolites
#' @param meta_pkn COSMOS meta prior knowledge network
#' @return tibble with pubchemID and cosmosID
prepare_metabolomics_names <- function(metabolomics_pubchem,meta_pkn){

	# add prefix
	working_data <- tibble(pubchemID = metabolomics_pubchem) %>%
		filter(!is.na(pubchemID)) %>%
		unique() %>%
		mutate(XMid = paste0("XMetab__",pubchemID))

	# find possible compartments in PKN
	compartment_codes <- unique(c(meta_pkn$source,meta_pkn$target))
	compartment_codes <- compartment_codes[grepl("Metab",compartment_codes)]
	compartment_codes <- unique(str_match(compartment_codes,"___.____"))

	# generate all combinations of metabolite and compartments and filter for
	# existing ones in the PKN
	pubchem_compart <- expand.grid(working_data$XMid,compartment_codes,stringsAsFactors = FALSE) %>%
		as_tibble() %>%
		mutate(cosmosID = paste0(Var1,Var2)) %>%
		rename(XMid = Var1) %>%
		select(XMid,cosmosID) %>%
		filter(cosmosID %in% c(meta_pkn$source,meta_pkn$target))

	# return tibble with pubchemID and possible cosmos id-s
	pubchem_compart %>% left_join(working_data,by = "XMid") %>%
		select(pubchemID,cosmosID)

}

#' convert gene symbols to entrez id
#'
#'
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
gene_symbols_to_entrez <- function(symbols){

	require(org.Hs.eg.db)
	require(AnnotationDbi)

	map_table <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, symbols, 'ENTREZID', 'SYMBOL')

	stopifnot(length(map_table) == length(symbols))
	stopifnot(all(names(map_table) ==symbols))

	map_table
}


#' convert gene symbols to entrez id
#'
#' @param ensembl vector of genes with ensembl id
#'
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
gene_ensembl_to_entrez <- function(ensembl){

	require(org.Hs.eg.db)
	require(AnnotationDbi)

	map_table <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, ensembl, 'ENTREZID', 'ENSEMBL')

	stopifnot(length(map_table) == length(ensembl))
	stopifnot(all(names(map_table) ==ensembl))

	map_table
}

#' find_unexpressed_genes
#'
#' finds the unexpressed
#'
find_unexpressed_genes <- function(expression_matrix, sample_design){
	require(edgeR)

	designMatrix <- model.matrix( ~ group + time + 0,  data = sample_design)
	keep <- edgeR::filterByExpr(as.matrix(expression_matrix), design = designMatrix)

}



#' filter_pkn_expressed_genes
#'
#' filters the non-expressed genes from the prior knowledge network
#'
#' @param expressed_genes_symbols  EntrezID of expressed genes
#' @param meta_pkn  COSMOS prior knowledge

filter_pkn_expressed_genes <- function(expressed_genes_entrez,meta_pkn){

	expressed_genes_entrez <- as.character(as.vector(as.matrix(expressed_genes_entrez)))

	is_expressed <- function(x)
	{
		if(!grepl("Metab",x))
		{
			if(gsub("X","",x) %in% expressed_genes_entrez)
			{
				return(x)
			} else
			{
				if(grepl("XGene[0-9]+__[0-9_]+$",x))
				{
					genes <- gsub("XGene[0-9]+__","",x)
					genes <- strsplit(genes,"_")[[1]]
					if(sum(genes %in% expressed_genes_entrez) != length(genes))
					{
						return(NA)
					} else
					{
						return(x)
					}
				} else
				{
					if(grepl("XGene[0-9]+__[0-9_]+reverse",x))
					{
						genes <- gsub("XGene[0-9]+__","",x)
						genes <- gsub("_reverse","",genes)
						genes <- strsplit(genes,"_")[[1]]
						if(sum(genes %in% expressed_genes_entrez) != length(genes))
						{
							return(NA)
						} else
						{
							return(x)
						}
					} else
					{
						return(NA)
					}
				}
			}
		} else
		{
			return(x)
		}
	}

	# is_expressed("XGene3004__124975_91227")

	meta_pkn$source <- sapply(meta_pkn$source,is_expressed)
	meta_pkn <- meta_pkn[complete.cases(meta_pkn),]
	meta_pkn$target <- sapply(meta_pkn$target,is_expressed)
	meta_pkn <- meta_pkn[complete.cases(meta_pkn),]

	return(meta_pkn)
}



#' format_COSMOS_res
#'
#' formats the network with readable names
#'
#' @param cosmos_res  results of CARNIVAL run
#' @param metab_mapping mapping table between metabolite naming conventions,
#' a two column dataframe with names: c("name","pubchem")
#' @param gene_mapping -- not sure !!!
#' @param measured_nodes vector of nodes that are measured or inputs
#' @param omnipath_ptm ptms database from OmnipathR
format_COSMOS_res <- function(cosmos_res,
							  metab_mapping,
							  gene_mapping = "ensembl",
							  measured_nodes,
							  omnipath_ptm)
{
	require(dorothea)

	sif <- as.data.frame(cosmos_res$weightedSIF)

	sif$Node1 <- gsub("^X","",sif$Node1)
	sif$Node2 <- gsub("^X","",sif$Node2)


	att <- as.data.frame(cosmos_res$nodesAttributes)#[,c(1,2)]
	names(att)[1] <- "Nodes"
	att$measured <- ifelse(att$Nodes %in% measured_nodes, 1, 0)
	att$Nodes <- gsub("^X","",att$Nodes)


	att$type <- ifelse(grepl("Metab",att$Nodes), "metabolite","protein")

	att <- att[abs(as.numeric(att$AvgAct)) > 0,]

	########################
	prots <- unique(att$Nodes)
	prots <- prots[!(grepl("Metab",prots))]
	prots <- gsub("Gene[0-9]+__","",prots)
	prots <- gsub("_reverse","",prots)
	prots <- gsub("EXCHANGE.*","",prots)
	prots <- unique(prots)
	prots <- unlist(sapply(prots, function(x){unlist(strsplit(x,split = "_"))}))


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


	sif$Node1 <- gsub("_reverse","",sif$Node1)
	sif$Node2 <- gsub("_reverse","",sif$Node2)

	att$Nodes <- gsub("_reverse","",att$Nodes)


	sif$Node1 <- gsub("EXCHANGE.*","",sif$Node1)
	sif$Node2 <- gsub("EXCHANGE.*","",sif$Node2)

	att$Nodes <- gsub("EXCHANGE.*","",att$Nodes)


	metabs <- unique(c(att$Nodes))
	metabs <- metabs[(grepl("Metab",metabs))]

	metab_to_pubchem <- as.data.frame(metab_mapping)
	metab_to_pubchem_vec <- metab_to_pubchem$name
	names(metab_to_pubchem_vec) <- metab_to_pubchem$pubchem

	for(i in 1:length(sif$Node1))
	{
		if(gsub("Gene[0-9]+__","",sif[i,1]) %in% names(gene_mapping))
		{
			if(grepl("Gene",sif[i,1]))
			{
				prefix <- gsub("__.*","",sif[i,1])
				sif[i,1] <- paste(prefix,gene_mapping[gsub("Gene[0-9]+__","",sif[i,1])], sep = "__")
			}
			else
			{
				sif[i,1] <- gene_mapping[gsub("Gene[0-9]+__","",sif[i,1])]
			}
		}
		if(gsub("Gene[0-9]+__","",sif[i,3]) %in% names(gene_mapping))
		{
			if(grepl("Gene",sif[i,3]))
			{
				prefix <- gsub("__.*","",sif[i,3])
				sif[i,3] <- paste(prefix,gene_mapping[gsub("Gene[0-9]+__","",sif[i,3])], sep = "__")
			}
			else
			{
				sif[i,3] <- gene_mapping[gsub("Gene[0-9]+__","",sif[i,3])]
			}
		}
		if(gsub("Metab__","",gsub("___[a-z]____","",sif[i,1])) %in% names(metab_to_pubchem_vec))
		{
			suffix <- str_match(sif[i,1],"___.____")
			metab <- metab_to_pubchem_vec[gsub("Metab__","",gsub("___[a-z]____","",sif[i,1]))]
			sif[i,1] <- paste(metab,suffix,sep = "")
		}
		if(gsub("Metab__","",gsub("___[a-z]____","",sif[i,3])) %in% names(metab_to_pubchem_vec))
		{
			suffix <- str_match(sif[i,3],"___.____")
			metab <- metab_to_pubchem_vec[gsub("Metab__","",gsub("___[a-z]____","",sif[i,3]))]
			sif[i,3] <- paste(metab,suffix,sep = "")
		}
		if(length(intersect(unlist(strsplit(gsub("Gene[0-9]+__","",sif[i,1]),split = "_")),names(gene_mapping))) > 0)
		{
			genes <- unlist(strsplit(gsub("Gene[0-9]+__","",sif[i,1]),split = "_"))
			genes <- unlist(sapply(genes, function(x){gene_mapping[x]}))
			genes <- paste0(genes, collapse = "_")
			sif[i,1] <- genes
		}
		if(length(intersect(unlist(strsplit(gsub("Gene[0-9]+__","",sif[i,3]),split = "_")),names(gene_mapping))) > 0)
		{
			genes <- unlist(strsplit(gsub("Gene[0-9]+__","",sif[i,3]),split = "_"))
			genes <- unlist(sapply(genes, function(x){gene_mapping[x]}))
			genes <- paste0(genes, collapse = "_")
			sif[i,3] <- genes
		}
	}


	for(i in 1:length(att$Nodes))
	{
		if(gsub("Gene[0-9]+__","",att[i,1]) %in% names(gene_mapping))
		{
			if(grepl("Gene",att[i,1]))
			{
				prefix <- gsub("__.*","",att[i,1])
				att[i,1] <- paste(prefix,gene_mapping[gsub("Gene[0-9]+__","",att[i,1])], sep = "__")
			}
			else
			{
				att[i,1] <- gene_mapping[gsub("Gene[0-9]+__","",att[i,1])]
			}
		}
		if(gsub("Metab__","",gsub("___[a-z]____","",att[i,1])) %in% names(metab_to_pubchem_vec))
		{
			suffix <- str_match(att[i,1],"___.____")
			metab <- metab_to_pubchem_vec[gsub("Metab__","",gsub("___[a-z]____","",att[i,1]))]
			att[i,1] <- paste(metab,suffix,sep = "")
		}
		if(length(intersect(unlist(strsplit(gsub("Gene[0-9]+__","",att[i,1]),split = "_")),names(gene_mapping))) > 0)
		{
			genes <- unlist(strsplit(gsub("Gene[0-9]+__","",att[i,1]),split = "_"))
			genes <- unlist(sapply(genes, function(x){gene_mapping[x]}))
			genes <- paste0(genes, collapse = "_")
			att[i,1] <- genes
		}
	}



	########################


	omnipath_ptm <- omnipath_ptm[omnipath_ptm$modification %in% c("dephosphorylation","phosphorylation"),]
	KSN <- omnipath_ptm[,c(4,3)]
	KSN$substrate_genesymbol <- paste(KSN$substrate_genesymbol,omnipath_ptm$residue_type, sep ="_")
	KSN$substrate_genesymbol <- paste(KSN$substrate_genesymbol,omnipath_ptm$residue_offset, sep = "")
	KSN$sign <- ifelse(omnipath_ptm$modification == "phosphorylation", 1, -1)

	att$type <- ifelse(att$Nodes %in% KSN$enzyme_genesymbol, "Kinase",att$type)

	dorothea <- as.data.frame(dorothea::dorothea_hs[dorothea::dorothea_hs$confidence %in% c("A","B","C","D"),c(3,1,4)])
	names(dorothea) <- c("target_genesymbol","source_genesymbol","sign")

	att$type <- ifelse(att$Nodes %in% dorothea$source_genesymbol, "TF",att$type)

	i <- 1
	for(node in att$Nodes) #check if a node name is basename appears more than once in nodes. If so, it is a metabolic enzyme
	{
		if(sum(gsub("Gene.*_","",node) == gsub("Gene.*_","",att$Nodes)) > 1)
		{
			att[i,"type"] <- "metab_enzyme"
		}
		i <- i+1
	}

	att$type <- ifelse(grepl("Gene.*_",att$Nodes), "metab_enzyme",att$type)


	att$Activity <- sign(as.numeric(as.character(att$AvgAct)))

	sif$Node1 <- gsub("Gene[0-9]+__","",sif$Node1)
	sif$Node2 <- gsub("Gene[0-9]+__","",sif$Node2)
	att$Nodes <- gsub("Gene[0-9]+__","",att$Nodes)


	sif <- sif[sif$Node1 != sif$Node2,]

	sif$Node1 <- gsub("[_]{3,5}","_",sif$Node1)
	sif$Node2 <- gsub("[_]{3,5}","_",sif$Node2)
	att$Nodes <- gsub("[_]{3,5}","_",att$Nodes)

	sif <- sif[sif$Node1 %in% att$Nodes & sif$Node2 %in% att$Nodes,]

	return(list(sif,att))
}



display_node_neighboorhood <- function(central_node,sif, att, n = 100)
{
	full_sif <- sif
	full_att <- att

	ig_net <- graph_from_data_frame(full_sif)

	ig_net <- make_ego_graph(ig_net, nodes = central_node, order = n, mode = "all")[[1]]

	to_keep <- V(ig_net)$name

	full_sif <- full_sif[full_sif$Node1 %in% to_keep & full_sif$Node2 %in% to_keep,]
	full_att <- full_att[full_att$Nodes %in% to_keep,]

	center_node <- shortest_paths(ig_net, from = central_node, to = full_att[full_att$measured == 1,1])
	center_node_out <- full_sif[full_sif$Node1 %in% names(unlist(center_node$vpath)) & full_sif$Node2 %in% names(unlist(center_node$vpath)),]

	# write_csv(center_node_net,"center_node_sif_newDoro.csv")

	center_node <- shortest_paths(ig_net, from = central_node, to = full_att[full_att$measured == 1,1], mode = "in")
	center_node_in <- full_sif[full_sif$Node1 %in% names(unlist(center_node$vpath)) & full_sif$Node2 %in% names(unlist(center_node$vpath)),]

	center_node_net <- as.data.frame(rbind(center_node_in,center_node_out))

	nodes <- full_att[full_att$Nodes %in% center_node_net$Node1 | full_att$Nodes %in% center_node_net$Node2,]
	edges <- center_node_net

	names(edges) <- c("from","to","sign","weigth")
	edges$color <- ifelse(edges$sign == 1, "green","red")
	edges$arrows <- "to"
	edges <- unique(edges)

	names(nodes)[1] <- "id"
	nodes$label <- nodes$id
	nodes$color <- ifelse(nodes$Activity > 0, "green","red")
	nodes <- nodes[!duplicated(nodes$id),]
	nodes$shape <- "dot"
	nodes[nodes$type == "metab_enzyme","shape"] <- "square"
	nodes[nodes$type == "protein","shape"] <- "square"
	nodes[nodes$type == "Kinase","shape"] <- "triangle"
	nodes[nodes$type == "TF","shape"] <- "diamond"
	nodes <- nodes[order(nodes$id),]
	nodes$shadow <- ifelse(nodes$measured == 1, T, F)

	return(visNetwork(nodes, edges, width = "100%") %>%
		   	visOptions(highlightNearest = TRUE,
		   			   nodesIdSelection = list(enabled = TRUE,
		   			   						style = 'width: 200px; height: 26px;
                                              background: #f8f8f8;
                                              color: darkblue;
                                              border:none;
                                              outline:none;')))

}
