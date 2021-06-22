# This script contains helper functions for summarizing the ancestral state
# reconstruction across multiple trees onto a single consensus tree. 


# ctreeAnc takes a multiphylo object (trees), a corresponding set of 
# marginal asr's from applying ancRECON across trees (clist), using rates 
# estimated from corHMM, and a consensus
# tree (ctree). It makes a datatable of probability of being in each state
# for nodes that are bifurcating in the ctree. These are indexed by their 
# corresponding ctree node. 
ctreeAnc <- function (trees, ctree, clist){
  # Return nodes in trees that match with ctree nodes
  matched.nodes <- lapply(trees, function (x) matchNodes(ctree, x)[,2])
  # Make vector of nodes in ctree 
  ctree.nodes <- unique(ctree$edge[,1])
  # Make vectors of nodes in trees. This will be used along with matched.nodes
  # to create an index of state probs to select. State probs are listed in same
  # order as nodes on the tree. 
  trees.nodes <- lapply(trees, function(x) unique(x$edge[,1]))
  # Get index of nodes that will be used to select state probs
  state.index <- mapply(match, matched.nodes, trees.nodes, SIMPLIFY = FALSE)
  # Select state probs using state.index. These are in order of ctree nodes
  state.probs <- mapply(function(x,y) x$lik.anc.states[y,], clist, state.index, 
                        SIMPLIFY = FALSE)
  # Add column that gives corresponding ctree node for each state prob
  state.probs <- lapply(state.probs, cbind, ctree.nodes)
  # Convert list of node probs to datatable 
  nodeframe <- rbindlist(lapply(state.probs, data.table))
  return(nodeframe)
}
