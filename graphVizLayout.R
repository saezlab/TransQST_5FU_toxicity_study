#' Create graphViz layout from graph
#'
#' Creates a matrix containing the graphViz layout for an igraph
#' object.
#'
#' @param g The graph to process.
#' @param mode The rank direction passed to the \link[Rgraphviz]{layoutGraph} function.
#'
#' @return A matrix with two columns containing the graph layout.
#'
#' @author Martin Garrido-Rodriguez, \email{mgrcprof@gmail.com}
#'
#' @importFrom graph graphNEL
#' @importFrom Rgraphviz layoutGraph
#' @importFrom igraph get.data.frame
#'
graphVizLayout <- function(g, mode = "LR") {

	# get simple graph format
	simple <- igraph::get.data.frame(g)[,1:2]
	# create graphNEL
	rEG <- new("graphNEL", nodes = V(g)$name, edgemode = "directed")
	for(i in 1:nrow(simple)){
		rEG <- graph::addEdge(simple[i,1], simple[i,2], rEG, 1)
	}
	# set graph atts
	att <- list(graph = list(rankdir = mode))
	# render layout
	rEG <- Rgraphviz::layoutGraph(rEG, attrs = att)
	# get matrix with X and Y
	l <- cbind(rEG@renderInfo@nodes$nodeX, y = rEG@renderInfo@nodes$nodeY)
	colnames(l) <- c("x","y")
	return(l)

}
