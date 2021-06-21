maxPaint <- function(tree, anc, rate.cat)
{
  
  ### This function takes a tree with branch lengths, and a set of ancestral
  ### state reconstructions made either with ancRECON or from summarizing
  ### across all trees in subsample. ASR should be a datatable with cols
  ### labeled with rate categories AS, ZS, etc. 
  ### The function returns a tree with branches painted according to the most 
  ### likely state at their ancestor node. 
  
  ### tree is a tree of class phylo with branch lengths
  ### anc is a datatable that lists the probability at each internal node of 
  ### each state. 
  ### rate.cat specifies the number of rate categories that were fit to the tree
  ### to estimate the rates used for the ancestral state reconstruction
  
  anc <- copy(anc)
  anc$lik.anc.states[,max :=  names(.SD)[max.col(.SD)]]
  if (rate.cat == 1)
  {
    
    as.edges <- unique(tree$edge[,1])[which(anc$lik.anc.states[,max] == "A")]
    as.paint <- unlist(Descendants(tree, node = as.edges, type = "children"))
    
    zs.edges <- unique(tree$edge[,1])[which(anc$lik.anc.states[,max] == "Z")]
    zs.paint <- unlist(Descendants(tree, node = zs.edges, type = "children"))
    
    plist <- list(as.paint, zs.paint)
    nvec <- c("A", "Z")
    
    for (i in 1:length(plist))
    {
      tryCatch(
        {
          tree <- paintBranches(tree, edge = plist[[i]], state = nvec[i])
        }, error = function(e)
        {
          message("Some regimes could not be printed on tree.")
        })
    }
    
  }
  else if (rate.cat == 2)
  {
    as.edges <- unique(tree$edge[,1])[which(anc$lik.anc.states[,max] == "AS")]
    as.paint <- unlist(Descendants(tree, node = as.edges, type = "children"))
    
    zs.edges <- unique(tree$edge[,1])[which(anc$lik.anc.states[,max] == "ZS")]
    zs.paint <- unlist(Descendants(tree, node = zs.edges, type = "children"))
    
    af.edges <- unique(tree$edge[,1])[which(anc$lik.anc.states[,max] == "AF")]
    af.paint <- unlist(Descendants(tree, node = af.edges, type = "children"))
    
    zf.edges <- unique(tree$edge[,1])[which(anc$lik.anc.states[,max] == "ZF")]
    zf.paint <- unlist(Descendants(tree, node = zf.edges, type = "children"))
    
    plist <- list(as.paint, zs.paint, af.paint, zf.paint)
    nvec <- c("AS", "ZS", "AF", "ZF")
    
    for (i in 1:length(plist))
    {
      tryCatch(
        {
          tree <- paintBranches(tree, edge = plist[[i]], state = nvec[i])
        }, error = function(e)
        {
          message("Some regimes could not be printed on tree.")
        })
    }
  } else
  {
    # Designate edges that should be painted in each state
    as.edges <- unique(tree$edge[,1])[which(anc$lik.anc.states[,max] == "AS")]
    as.paint <- unlist(Descendants(tree, node = as.edges, type = "children"))
    
    zs.edges <- unique(tree$edge[,1])[which(anc$lik.anc.states[,max] == "ZS")]
    zs.paint <- unlist(Descendants(tree, node = zs.edges, type = "children"))
    
    am.edges <- unique(tree$edge[,1])[which(anc$lik.anc.states[,max] == "AM")]
    am.paint <- unlist(Descendants(tree, node = am.edges, type = "children"))
    
    zm.edges <- unique(tree$edge[,1])[which(anc$lik.anc.states[,max] == "ZM")]
    zm.paint <- unlist(Descendants(tree, node = zm.edges, type = "children"))
    
    af.edges <- unique(tree$edge[,1])[which(anc$lik.anc.states[,max] == "AF")]
    af.paint <- unlist(Descendants(tree, node = af.edges, type = "children"))
    
    zf.edges <- unique(tree$edge[,1])[which(anc$lik.anc.states[,max] == "ZF")]
    zf.paint <- unlist(Descendants(tree, node = zf.edges, type = "children"))
    
    plist <- list(as.paint, zs.paint, am.paint, zm.paint, af.paint, zf.paint)
    nvec <- c("AS", "ZS", "AM", "ZM", "AF", "ZF")
    
    for (i in 1:length(plist))
    {
      tryCatch(
        {
          tree <- paintBranches(tree, edge = plist[[i]], state = nvec[i])
        }, error = function(e)
        {
          message("Some regimes could not be printed on tree.")
        })
    }
  }
  return(tree)
}

treePaint <- function(tree, anc, rate.cat)
{
  ### This function takes a tree with branch lengths, and a set of ancestral
  ### state reconstructions made either with ancRECON or from summarizing
  ### across all trees in subsample. ASR should be a datatable with cols
  ### labeled with rate categories AS, ZS, etc. 
  ### The function returns a tree with branches painted according to the most 
  ### likely state at their ancestor node. Unlike the maxPaint function, this 
  ### uses a cutoff. If the probability of any given state is less than 0.75,
  ### the branch is painted grey to denote uncertainty. 
  
  ### tree is a tree of class phylo with branch lengths
  ### anc is a datatable that lists the probability at each internal node of 
  ### each state. 
  ### rate.cat specifies the number of rate categories that were fit to the tree
  ### to estimate the rates used for the ancestral state reconstruction
  
  if (rate.cat == 1)
  {
    as.edges <- unique(tree$edge[,1])[which(anc$lik.anc.states[,1]>0.75)]
    as.paint <- unlist(Descendants(tree, node = as.edges, type = "children"))
    
    zs.edges <- unique(tree$edge[,1])[which(anc$lik.anc.states[,2]>0.75)]
    zs.paint <- unlist(Descendants(tree, node = zs.edges, type = "children"))
    
    plist <- list(as.paint, zs.paint)
    nvec <- c("A", "Z")
    
    for (i in 1:length(plist))
    {
      tryCatch(
        {
          tree <- paintBranches(tree, edge = plist[[i]], state = nvec[i])
        }, error = function(e)
        {
          message("Some regimes could not be printed on tree.")
        })
    }
    
  }
  else if (rate.cat == 2)
  {
    as.edges <- unique(tree$edge[,1])[which(anc$lik.anc.states[,1]>0.75)]
    as.paint <- unlist(Descendants(tree, node = as.edges, type = "children"))
    
    zs.edges <- unique(tree$edge[,1])[which(anc$lik.anc.states[,2]>0.75)]
    zs.paint <- unlist(Descendants(tree, node = zs.edges, type = "children"))
    
    af.edges <- unique(tree$edge[,1])[which(anc$lik.anc.states[,3]>0.75)]
    af.paint <- unlist(Descendants(tree, node = af.edges, type = "children"))
    
    zf.edges <- unique(tree$edge[,1])[which(anc$lik.anc.states[,4]>0.75)]
    zf.paint <- unlist(Descendants(tree, node = zf.edges, type = "children"))
    
    plist <- list(as.paint, zs.paint, af.paint, zf.paint)
    nvec <- c("AS", "ZS", "AF", "ZF")
    
    for (i in 1:length(plist))
    {
      tryCatch(
        {
          tree <- paintBranches(tree, edge = plist[[i]], state = nvec[i])
        }, error = function(e)
        {
          message("Some regimes could not be printed on tree.")
        })
    }
  } else
  {
    # Designate edges that should be painted in each state
    as.edges <- unique(tree$edge[,1])[which(anc$lik.anc.states[,1]>0.75)]
    as.paint <- unlist(Descendants(tree, node = as.edges, type = "children"))
    
    zs.edges <- unique(tree$edge[,1])[which(anc$lik.anc.states[,2]>0.75)]
    zs.paint <- unlist(Descendants(tree, node = zs.edges, type = "children"))
    
    am.edges <- unique(tree$edge[,1])[which(anc$lik.anc.states[,3]>0.75)]
    am.paint <- unlist(Descendants(tree, node = am.edges, type = "children"))
    
    zm.edges <- unique(tree$edge[,1])[which(anc$lik.anc.states[,4]>0.75)]
    zm.paint <- unlist(Descendants(tree, node = zm.edges, type = "children"))
    
    af.edges <- unique(tree$edge[,1])[which(anc$lik.anc.states[,5]>0.75)]
    af.paint <- unlist(Descendants(tree, node = af.edges, type = "children"))
    
    zf.edges <- unique(tree$edge[,1])[which(anc$lik.anc.states[,6]>0.75)]
    zf.paint <- unlist(Descendants(tree, node = zf.edges, type = "children"))
    
    plist <- list(as.paint, zs.paint, am.paint, zm.paint, af.paint, zf.paint)
    nvec <- c("AS", "ZS", "AM", "ZM", "AF", "ZF")
    
    for (i in 1:length(plist))
    {
      tryCatch(
        {
          tree <- paintBranches(tree, edge = plist[[i]], state = nvec[i])
        }, error = function(e)
        {
          message("Some regimes could not be printed on tree.")
        })
    }
  }
  return(tree)
}

azPaint <- function(tree, anc, rate.cat)
{
  ### This function takes a tree with branch lengths, and a set of ancestral
  ### state reconstructions made either with ancRECON or from summarizing
  ### across all trees in subsample. ASR should be a datatable with cols
  ### labeled with rate categories AS, ZS, etc. The function sums across 
  ### all zoox states and azoox states to plot an asr of just zoox vs azoox
  ### The function returns a tree with branches painted according to whether
  ### their ancestral node is zoox or azoox
  
  ### tree is a tree of class phylo with branch lengths
  ### anc is a datatable that lists the probability at each internal node of 
  ### each state. 
  ### rate.cat specifies the number of rate categories that were fit to the tree
  ### to estimate the rates used for the ancestral state reconstruction
  
  anc <- copy(anc)
  
  if (rate.cat == 1)
  {
    z.edges <- unique(tree$edge[,1])[which(anc[,Z]>0.75)]
    z.paint <- unlist(Descendants(tree, node = z.edges, type = "children"))
    
    a.edges <- unique(tree$edge[,1])[which(anc[,A]>0.75)]
    a.paint <- unlist(Descendants(tree, node = a.edges, type = "children"))
    
    tree <- paintBranches(tree, edge = z.paint, state = "Z")
    tree <- paintBranches(tree, edge = a.paint, state = "A")
    
    outlist <- list(tree = tree, nodeprobs = anc[, .(A = A, Z = Z)])
    
  } else
  {
    anc[, zsum := rowSums(.SD), .SDcols = grep("Z", names(anc))]
    anc[, asum := rowSums(.SD), .SDcols = grep("A", names(anc))]
    
    z.edges <- unique(tree$edge[,1])[which(anc[,zsum]>0.75)]
    z.paint <- unlist(Descendants(tree, node = z.edges, type = "children"))
    
    a.edges <- unique(tree$edge[,1])[which(anc[,asum]>0.75)]
    a.paint <- unlist(Descendants(tree, node = a.edges, type = "children"))
    
    tree <- paintBranches(tree, edge = z.paint, state = "Z")
    tree <- paintBranches(tree, edge = a.paint, state = "A")
    
    outlist <- list(tree = tree, nodeprobs = anc[, .(A = asum, Z = zsum)])
    
  }
  return(outlist)
}

