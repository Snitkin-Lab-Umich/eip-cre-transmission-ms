# from hogwash

is_tip <- function(node_num, tr){
  return(node_num <= ape::Ntip(tr))
}

identify_transition_edges <- function(tr, mat, num, node_recon, disc_cont, same_thresh = 1, disc_bin = TRUE){

  # FUNCTION -------------------------------------------------------------------
  transition <- transition_direction <-
    parent_node <- child_node <- integer(ape::Nedge(tr))
  older <- 1 # older node is 1st column in tr$edge
  younger <- 2 # younger node is 2nd column in tr$edge
  parent_0_child_1 <- 1
  parent_1_child_0 <- -1
  parent_equals_child <- 0
  both_parent_and_child_are_one <- 2

  for (i in 1:ape::Nedge(tr)) {
    if (is_tip(tr$edge[i, older], tr)) {
      stop("tree invalid")
    }
    parent_node[i] <- node_recon[tr$edge[i, older] - ape::Ntip(tr)]
    if (is_tip(tr$edge[i, younger], tr)) {
      # child is a tip
      child_node[i]  <- mat[, num][tr$edge[i, younger]]
    } else {
      # child is internal nodes
      child_node[i]  <- node_recon[tr$edge[i, younger] - ape::Ntip(tr)]
    }
    if(!disc_bin){
      nodes_diff = parent_node[i] != child_node[i]
      if(is.na(parent_node[i]) | is.na(child_node[i])){
        transition[i] <- NA
      }else if(!nodes_diff){
        transition[i] <- parent_node[i]
      }else{
        transition[i] <- 'trans'#paste0(parent_node[i],'-',child_node[i])

      }
      transition_direction[i] <- NA
      next
    }
    transition[i] <- sum(parent_node[i] + child_node[i])
    # transition[i] is either 0, 1, or 2 for discrete traits
    # transition[i] is not to be used when the trait is continuous because all,
    # or very nearly all edges are transition edges.

    if (parent_node[i] > child_node[i]) {
      transition_direction[i] <- parent_1_child_0
    } else if (parent_node[i] < child_node[i]) {
      transition_direction[i] <- parent_0_child_1
    }
    if(disc_cont == "continuous"){
      diff_parent_child <- abs(parent_node[i] - child_node[i])
      if (diff_parent_child < same_thresh) {
        transition_direction[i] <- parent_equals_child
      }
    }
    if (disc_cont == "discrete") {
      transition[transition == both_parent_and_child_are_one] <-
        parent_equals_child # parent_node == child_node, then no transition.
    } else {
      transition <- NA # Not used when the trait is continuous.
    }
  }
  # Check and return output ----------------------------------------------------
  results <- list("transition" = transition, "trans_dir" = transition_direction)
  return(results)
}

# new functions
transition_edge_pipeline <- function(df, tree, continuous = TRUE, same_thresh = 1, ar_conf = 0.875){

  type = 'discrete'
  if(continuous) type = 'continuous'

  # can't have nas
  df = df %>% drop_na()
  tree = drop.tip(tree,tree$tip.label[!tree$tip.label %in% df$isolate_no])

  ar = lapply(2:ncol(df), function(x) {
    dat = df %>% select(isolate_no,{{x}})
    dat_vec = unname(unlist(c(dat[,2])))
    names(dat_vec) = dat$isolate_no
    ape::ace(x = dat_vec, phy = tree, type = type)
  })
  names(ar) = names(df)[2:ncol(df)]

  node_recons = lapply(ar, function(p){
    if(continuous){
      nrs = p$ace
    }else{
      lik = p$lik.anc
      nrs = apply(lik,1,function(x) ifelse(x[1]>x[2],0,1))
    }
    return(nrs)
  })

  if(is.null(tree$node.label)){
    warning('No edge confidences in tree.\n')
    edge_confs = rep(100, (Nnode(tree) + Ntip(tree)))
  }else{
    edge_confs = as.numeric(sapply(tree$edge[,1], function(x) tree$node.label[x-Ntip(tree)]))
    # make root edge confidence zero (not sure if this is good?)
    edge_confs[is.na(edge_confs)] = 0
  }

  ar_confs = lapply(ar, function(p){
    if(continuous){
      recon = p$ace
      ci = p$CI95
      ci[ci < 0] = 10e-10
      lb = ci[,1]
      ub = ci[,2]
      lb_ci_diff = recon - lb
      ub_ci_diff = ub - recon
      confs = lb_ci_diff + ub_ci_diff
      return(confs)
    }else{
      lik = p$lik.anc
      nr_inds = apply(lik,1,function(x) ifelse(x[1]>x[2],1,2))
      confs = sapply(1:length(nr_inds),function(x) ifelse(unname(lik[x,nr_inds[x]]) <= ar_conf,0,1))
    }
  })

  ar_confs_edge = lapply(ar_confs, function(conf){
    apply(tree$edge, 1, function(x){
      parent = x[1]-Ntip(tree)
      parent_conf = unlist(conf[parent])
      if(x[2] > Ntip(tree)){
        child = x[2]-Ntip(tree)
        child_conf = unlist(conf[child])
      }else{
        child = x[2]
        child_conf = 1
      }
      if(parent_conf < 1 | child_conf < 1) return(0)
      else return(1)
    })
  })

  transitions = lapply(1:(ncol(df)-1), function(x) {
    if(length(unlist(node_recons[[x]])) != Nnode(tree)) return(NULL)
    trans = identify_transition_edges(tree, data.frame(df[,x+1]), 1, node_recons[[x]], type, same_thresh = same_thresh)
    if(continuous){
      trans$transition[edge_confs <= 70 | ar_confs_edge[[x]] > 1] = NA
      trans$trans_dir[edge_confs <= 70 | ar_confs_edge[[x]] > 1] = NA
    }else{
      trans$transition[edge_confs <= 70 | ar_confs_edge[[x]] <= ar_conf] = NA
      trans$trans_dir[edge_confs <= 70 | ar_confs_edge[[x]] <= ar_conf] = NA
    }

    return(trans)
  })
  names(transitions) = names(df)[2:ncol(df)]

  return(transitions)

}
