# Helper functions for counting number of gains and losses from ancestral state
# recontructions. 

getChanges <- function(tree, states, rate.cat) 
{
  
  # This function takes a tree and associated states at internal nodes and tips
  # and counts the number of transitions between states. Tree must be a single 
  # tree in phylo format. states is a numeric vector with tip states first and
  # node states second. Tip states must be in order of tip labels of tree. Node
  # states must be in order of internal nodes. rate.cat must be an integer 1-3.
  # Returns a dataframe where each column is the number of transitions from or
  # to each state. 
  
  if(rate.cat == 1)
  {
    al <- sum(apply(tree$edge,1,function(x,s) 
      s[x[1]] == 1 & s[x[1]]!=s[x[2]], s=states))
    
    ag <- sum(apply(tree$edge,1,function(x,s) 
      s[x[1]] != 1 & s[x[2]] == 1, s=states))
    
    zl <- sum(apply(tree$edge,1,function(x,s) 
      s[x[1]] == 2 & s[x[1]]!=s[x[2]], s=states))
    
    zg <- sum(apply(tree$edge,1,function(x,s) 
      s[x[1]] != 2 & s[x[2]] == 2,s=states))
    
    changes <- data.table(A.losses = al, A.gains = ag, 
                          Z.losses = zl, Z.gains = zg)
  } 
  else if (rate.cat == 2)
  {
    # AS losses
    asl <- sum(apply(tree$edge,1,function(x,s) 
      s[x[1]] == 1 & s[x[1]]!=s[x[2]], s=states))
    
    # AS gains
    asg <- sum(apply(tree$edge,1,function(x,s) 
      s[x[1]] != 1 & s[x[2]] == 1, s=states))
    
    # ZS losses
    zsl <- sum(apply(tree$edge,1,function(x,s) 
      s[x[1]] == 2 & s[x[1]]!=s[x[2]], s=states))
    
    # ZS gains
    zsg <- sum(apply(tree$edge,1,function(x,s) 
      s[x[1]] != 2 & s[x[2]] == 2, s=states))
    
    # AF losses
    afl <- sum(apply(tree$edge,1,function(x,s) 
      s[x[1]] == 3 & s[x[1]]!=s[x[2]], s=states))
    
    # AF gains
    afg <- sum(apply(tree$edge,1,function(x,s) 
      s[x[1]] != 3 & s[x[2]] == 3, s=states))
    
    # ZF losses
    zfl <- sum(apply(tree$edge,1,function(x,s) 
      s[x[1]] == 4 & s[x[1]]!=s[x[2]], s=states))
    
    # ZF gains 
    zfg <- sum(apply(tree$edge,1,function(x,s) 
      s[x[1]] != 4 & s[x[2]] == 4, s=states))
    
    changes <- data.table(AS.losses = asl, AS.gains = asg, ZS.losses = zsl, 
                          ZS.gains = zsg, AF.losses = afl, AF.gains = afg, 
                          ZF.losses = zfl, ZF.gains = zfg)
  } 
  else if (rate.cat == 3)
  {
    # AS losses
    asl <- sum(apply(tree$edge,1,function(x,s) 
      s[x[1]] == 1 & s[x[1]]!=s[x[2]], s=states))
    
    # AS gains
    asg <- sum(apply(tree$edge,1,function(x,s) 
      s[x[1]] != 1 & s[x[2]] == 1, s=states))
    
    # ZS losses
    zsl <- sum(apply(tree$edge,1,function(x,s) 
      s[x[1]] == 2 & s[x[1]]!=s[x[2]], s=states))
    
    # ZS gains
    zsg <- sum(apply(tree$edge,1,function(x,s) 
      s[x[1]] != 2 & s[x[2]] == 2, s=states))
    
    # AM losses
    aml <- sum(apply(tree$edge,1,function(x,s) 
      s[x[1]] == 3 & s[x[1]]!=s[x[2]], s=states))
    
    # AM gains
    amg <- sum(apply(tree$edge,1,function(x,s) 
      s[x[1]] != 3 & s[x[2]] == 3, s=states))
    
    # ZM losses
    zml <- sum(apply(tree$edge,1,function(x,s) 
      s[x[1]] == 4 & s[x[1]]!=s[x[2]], s=states))
    
    # ZM gains 
    zmg <- sum(apply(tree$edge,1,function(x,s) 
      s[x[1]] != 4 & s[x[2]] == 4, s=states))
    
    # AF losses
    afl <- sum(apply(tree$edge,1,function(x,s) 
      s[x[1]] == 5 & s[x[1]]!=s[x[2]], s=states))
    
    # AF gains
    afg <- sum(apply(tree$edge,1,function(x,s) 
      s[x[1]] != 5 & s[x[2]] == 5, s=states))
    
    # ZF losses
    zfl <- sum(apply(tree$edge,1,function(x,s) 
      s[x[1]] == 6 & s[x[1]]!=s[x[2]], s=states))
    
    # ZF gains 
    zfg <- sum(apply(tree$edge,1,function(x,s) 
      s[x[1]] != 6 & s[x[2]] == 6, s=states))
    
    changes <- data.table(AS.losses = asl, AS.gains = asg, 
                          ZS.losses = zsl, ZS.gains = zsg,
                          AM.losses = aml, AM.gains = amg,
                          ZM.losses = zml, ZM.gains = zmg,
                          AF.losses = afl, AF.gains = afg, 
                          ZF.losses = zfl, ZF.gains = zfg)
    
  }
  
  return(changes)
}