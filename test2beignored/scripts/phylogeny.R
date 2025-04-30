sp_list <- data_ordered_seeds
taxon <- "Plant"
scenario <- "at_basal_node"
show_grafted <- TRUE
tree_by_user <- FALSE
mc_cores <- future::availableCores() - 2
.progress <- "text"
fish_tree <- c("timetree", "all-taxon")
mammal_tree <- c("vertlife", 
                 "phylacine")
bee_tree <- c("maximum-likelihood", "bootstrap")
dt <- TRUE

tree <- megatrees::tree_plant_otl
library(fastmatch)
tree_genus <- unique(gsub("^([-A-Za-z]*)_.*$", "\\1", tree$tip.label))
sp_list <- rtrees::sp_list_df(unique(sp_list))
all_genus_in_tree <- all(unique(sp_list$genus) %fin% tree_genus)


  if (!is.null(taxon)) {
    if (!taxon %fin% rtrees::taxa_supported & !all_genus_in_tree) {
      new_cls = unique(dplyr::select(sp_list, genus, family))
      new_cls$taxon = taxon
      classifications <- dplyr::bind_rows(rtrees::classifications, 
                                          new_cls)
    }
  }
  sp_out_tree <- sp_list[!sp_list$species %fin% tree$tip.label,] #species not in the megatree
  subsp_in_tree <- grep("^.*_.*_.*$", x = tree$tip.label, value = T)
  if (length(subsp_in_tree)) {
    sp_out_tree = dplyr::mutate(sp_out_tree, re_matched = NA, 
                                matched_name = NA)
    for (i in 1:length(sp_out_tree$species)) {
      name_in_tree = grep(paste0("^", sp_out_tree$species[i], 
                                 "_"), x = subsp_in_tree, ignore.case = T, value = T)
      if (length(name_in_tree)) {
        sp_out_tree$re_matched[i] = TRUE
        sp_out_tree$matched_name[i] = sample(name_in_tree, 
                                             1)
        tree$tip.label[tree$tip.label == sp_out_tree$matched_name[i]] = sp_out_tree$species[i]
        if (!is.null(tree$genus_family_root)) {
          tree$genus_family_root$only_sp[tree$genus_family_root$only_sp == 
                                           sp_out_tree$matched_name[i]] = sp_out_tree$species[i]
        }
      }
    }
    sp_out_tree = dplyr::distinct(sp_list[!sp_list$species %fin% 
                                            tree$tip.label, ])
  }
  close_sp_specified = close_genus_specified = FALSE
  if ("close_sp" %fin% names(sp_out_tree)) {
    close_sp_specified = TRUE
    sp_out_tree$close_sp = cap_first_letter(gsub(" +", "_", 
                                                 sp_out_tree$close_sp))
  }
  if ("close_genus" %fin% names(sp_out_tree)) 
    close_genus_specified = TRUE
  # if (nrow(sp_out_tree) == 0) {
  #   message("Wow, all species are already in the mega-tree!")
  #   tree_sub = ape::drop.tip(tree, setdiff(tree$tip.label, 
  #                                          sp_list$species))
  #   return(tree_sub)
  # }
  if (is.null(tree$genus_family_root)) 
    stop("Did you use your own phylogeny? If so, please set `tree_by_user = TRUE`.")
  sp_out_tree$status = ""
  tree_df = tidytree::as_tibble(tree)
  tree_df$is_tip = !(tree_df$node %fin% tree_df$parent)
  node_hts = ape::branching.times(tree)
  all_eligible_nodes = unique(c(tree$genus_family_root$basal_node, 
                                tree$genus_family_root$root_node))
  n_spp_to_show_progress = 200
  if (nrow(sp_out_tree) > n_spp_to_show_progress) {
    progress <- create_progress_bar(.progress)
    progress$init(nrow(sp_out_tree))
    on.exit(progress$term())
  }
  for (i in 1:nrow(sp_out_tree)) {
    if (nrow(sp_out_tree) > n_spp_to_show_progress) 
      progress$step()
    where_loc_i = where_loc_i2 = NA
    if (close_sp_specified) {
      if (!is.na(sp_out_tree$close_sp[i]) & sp_out_tree$close_sp[i] %fin% 
          tree$tip.label) {
        where_loc_i = sp_out_tree$close_sp[i]
      }
    }
    if (close_genus_specified) {
      if (!is.na(sp_out_tree$close_genus[i]) & sp_out_tree$close_genus[i] != 
          "" & sp_out_tree$close_genus[i] %fin% tree_genus) {
        sp_out_tree$genus[i] = sp_out_tree$close_genus[i]
        where_loc_i2 = sp_out_tree$close_genus[i]
      }
      else {
        if (!is.na(sp_out_tree$close_genus[i])) 
          warning("The genus specified for ", sp_out_tree$species[i], 
                  " is not in the phylogeny.")
      }
    }
    if (!all_genus_in_tree & is.na(where_loc_i) & is.na(where_loc_i2)) {
      if (is.na(sp_out_tree$family[i]) | !sp_out_tree$family[i] %fin% 
          tree$genus_family_root$family) {
        sp_out_tree$status[i] = "No co-family species in the mega-tree"
        (next)()
      }
    }
    node_label_new = NULL
    add_above_node = FALSE
    fraction = 1/2
    if (sp_out_tree$genus[i] %fin% tree$genus_family_root$genus | 
        !is.na(where_loc_i2) | !is.na(where_loc_i)) {
      sp_out_tree$status[i] = "*"
      idx_row = which(tree$genus_family_root$genus == 
                        sp_out_tree$genus[i])
      root_sub = tree$genus_family_root[idx_row, ]
      if (root_sub$n_spp == 1 | !is.na(where_loc_i)) {
        if (!is.na(where_loc_i)) {
          where_loc = where_loc_i
          new_ht = tree_df$branch.length[tree_df$label == 
                                           where_loc_i] * (1 - fraction)
          node_hts = c(new_ht, node_hts)
          node_label_new = paste0("N", length(node_hts))
          names(node_hts)[1] = node_label_new
          all_eligible_nodes = c(all_eligible_nodes, 
                                 node_label_new)
          add_above_node = TRUE
          if (!sp_out_tree$genus[i] %fin% tree$genus_family_root$genus) {
            tree$genus_family_root = tibble::add_row(tree$genus_family_root, 
                                                     family = sp_out_tree$family[i], genus = sp_out_tree$genus[i], 
                                                     basal_node = node_label_new, basal_time = new_ht, 
                                                     root_node = tree_df$label[tree_df$node == 
                                                                                 tree_df$parent[tree_df$label == where_loc_i]], 
                                                     root_time = tree_df$branch.length[tree_df$node == 
                                                                                         tree_df$parent[tree_df$label == where_loc_i]], 
                                                     n_genus = 1, n_spp = 1, only_sp = sp_out_tree$species[i])
          }
        }
        else {
          where_loc = root_sub$only_sp
          new_ht = root_sub$basal_time * (1 - fraction)
          node_hts = c(new_ht, node_hts)
          node_label_new = paste0("N", length(node_hts))
          names(node_hts)[1] = node_label_new
          all_eligible_nodes = c(all_eligible_nodes, 
                                 node_label_new)
          add_above_node = TRUE
          tree$genus_family_root$only_sp[idx_row] = NA
          tree$genus_family_root$basal_node[idx_row] = node_label_new
          tree$genus_family_root$basal_time[idx_row] = unname(new_ht)
        }
      }
      else {
        where_loc = root_sub$basal_node
        if (scenario == "random_below_basal") {
          tree_df_sub = tidytree::offspring(tree_df, 
                                            where_loc)
          tree_df_sub = tree_df_sub[tree_df_sub$is_tip == 
                                      FALSE, ]
          if (nrow(tree_df_sub) > 0) {
            potential_locs = c(where_loc, tree_df_sub$label)
            bls = tree_df_sub$branch.length
            names(bls) = tree_df_sub$label
            bls = c(root_sub$root_time - root_sub$basal_time, 
                    bls)
            names(bls)[1] = root_sub$basal_node
            prob = bls/sum(bls)
            where_loc = sample(potential_locs, 1, prob = prob)
          }
        }
      }
    }
    else {
      sp_out_tree$status[i] = "**"
      idx_row = which(tree$genus_family_root$family == 
                        sp_out_tree$family[i] & is.na(tree$genus_family_root$genus))
      root_sub = tree$genus_family_root[idx_row, ]
      if (root_sub$n_spp == 1) {
        where_loc = root_sub$only_sp
        new_ht = root_sub$basal_time * (1 - fraction)
        node_hts = c(new_ht, node_hts)
        node_label_new = paste0("N", length(node_hts))
        names(node_hts)[1] = node_label_new
        all_eligible_nodes = c(all_eligible_nodes, node_label_new)
        add_above_node = TRUE
        tree$genus_family_root = tibble::add_row(tree$genus_family_root, 
                                                 family = sp_out_tree$family[i], genus = sp_out_tree$genus[i], 
                                                 basal_node = node_label_new, basal_time = unname(new_ht), 
                                                 root_node = node_label_new, root_time = unname(new_ht), 
                                                 n_genus = 1, n_spp = 1, only_sp = sp_out_tree$species[i])
      }
      else {
        where_loc = root_sub$basal_node
        if (scenario == "random_below_basal") {
          tree_df_sub = tidytree::offspring(tree_df, 
                                            where_loc)
          tree_df_sub = tree_df_sub[tree_df_sub$is_tip == 
                                      FALSE, ]
          if (nrow(tree_df_sub) > 0) {
            potential_locs = intersect(c(where_loc, 
                                         tree_df_sub$label), all_eligible_nodes)
            locs_bl = tree_df_sub[tree_df_sub$label %fin% 
                                    potential_locs, ]
            bls = locs_bl$branch.length
            names(bls) = locs_bl$label
            bls = c(root_sub$root_time - root_sub$basal_time, 
                    bls)
            names(bls)[1] = root_sub$basal_node
            prob = bls/sum(bls)
            where_loc = sample(potential_locs, 1, prob = prob)
          }
        }
      }
      tree$genus_family_root$n_genus[idx_row] = tree$genus_family_root$n_genus[idx_row] + 
        1
    }
    if (root_sub$n_spp > 3) 
      use_castor = TRUE
    else use_castor = FALSE
    if (dt) {
      tree_df = rtrees::bind_tip(tree_tbl = tree_df, node_heights = node_hts, 
                         where = where_loc, new_node_above = add_above_node, 
                         tip_label = sp_out_tree$species[i], frac = fraction, 
                         return_tree = FALSE, node_label = node_label_new, 
                         use_castor = use_castor)
    }
    tree$genus_family_root$n_spp[idx_row] = tree$genus_family_root$n_spp[idx_row] + 1
  }
  tree_df = dplyr::arrange(tree_df, node)
  if (any(sp_out_tree$status == "*")) {
    message("\n", sum(sp_out_tree$status == "*"), " species added at genus level (*) \n")
  }
  if (any(sp_out_tree$status == "**")) {
    message(sum(sp_out_tree$status == "**"), " species added at family level (**) \n")
  }
  # tree_sub = castor::get_subtree_with_tips(tidytree::as.phylo(tree_df), 
  #                                          sp_list$species)$subtree
  tree_sub <- tidytree::as.phylo(tree_df)
  grafted = sp_out_tree[sp_out_tree$status %fin% c("*", "**"), 
  ]
  grafted$sp2 = paste0(grafted$species, grafted$status)
  wid = which(tree_sub$tip.label %fin% grafted$species)
  tree_sub$tip.label[wid] = dplyr::left_join(tibble::tibble(species = tree_sub$tip.label[wid]), 
                                             grafted, by = "species")$sp2
  graft_status = tibble::tibble(tip_label = tree_sub$tip.label)
  graft_status$species = gsub("\\*", "", graft_status$tip_label)
  graft_status$status = ifelse(grepl("\\*{2}$", graft_status$tip_label), 
                               "grafted at family level", ifelse(grepl("[^*]\\*{1}$", 
                                                                       graft_status$tip_label), "grafted at genus level", 
                                                                 "exisiting species in the megatree"))
  if (any(sp_out_tree$status == "No co-family species in the mega-tree")) {
    sp_no_family = sp_out_tree$species[sp_out_tree$status == 
                                         "No co-family species in the mega-tree"]
    message(length(sp_no_family), " species have no co-family species in the mega-tree, skipped\n(if you know their family, prepare and edit species list with `rtrees::sp_list_df()` may help): \n", 
            paste(sp_no_family, collapse = ", "))
    graft_status = dplyr::bind_rows(graft_status, data.frame(species = sp_no_family, 
                                                             status = rep("skipped as no co-family in the megatree", 
                                                                          length(sp_no_family))))
  }
  if (!show_grafted) {
    tree_sub <- rm_stars(tree_sub)
    graft_status$tip_label = gsub("\\*", "", graft_status$tip_label)
  }
  tree_sub$graft_status = graft_status
  tree_sub <- ape::ladderize(tree_sub)

  
  ## get a list of all genera
  tips <- tree_sub$tip.label
  genera <- unique(sapply(strsplit(tips,"_"),function(x) x[1]))

  ## now drop all but one of each
  ii <- sapply(genera,function(x,y) grep(x,y)[1],y=tips)
  tree_oneper_genus <- ape::drop.tip(tree_sub,setdiff(tree_sub$tip.label,tips[ii]))
  ggtree::ggtree(tree_oneper_genus,layout = "fan")
  
  tree_oneper_genus$tip.label<-sapply(strsplit(tree_oneper_genus$tip.label,"_"),function(x) x[1])
  
  tab_tree <- as_tibble(tree_oneper_genus)
  tab_tree%<>%
    mutate(is_in_our_dataset=ifelse(label%in%data_ordered_seeds$genus,"yes","nope"),
           node_depth= ape::node.depth(ape::as.phylo(tab_tree)))%>%
    mutate(deep_node_label= ifelse(node_depth>1&!grepl("N[0-9]*$|mrcaott.*|ae$",label),stringr::str_extract(label,"^[a-zA-Z]*"),NA))
  
  
  test <-
  tab_tree%>%
    filter(is_in_our_dataset=="yes")
  pathtocollapse <- NULL
  for(nodez in test$node){
    path2genus <- ggtree::get.path(ape::as.phylo(tab_tree),10515,nodez)
    pathtocollapse <- c(pathtocollapse,path2genus[-length(path2genus)])
  }
  nodeztocollapse <- unique(pathtocollapse)
  
  # outgroupz <- tab_tree%>%
  #   arrange(desc(branch.length))%>%
  #   slice(1:10)
  # nodeztocollapse <- c(nodeztocollapse,outgroupz$node)
  p <- ape::as.phylo(tab_tree)%>%
    ggtree::ggtree(tab_tree)
  
  library(ggtree)
  p2 <- p  %<+% tab_tree + 
    ggtree::geom_tippoint(aes(color = is_in_our_dataset))+
    ggtree::  geom_nodelab(
      mapping = aes(
        label = deep_node_label),
      fill = "lightgrey",
      geom = "label")
p3 <-  purrr::reduce(
  tab_tree$node[!tab_tree$node%in%nodeztocollapse],
  \(x,y) collapse(x,y),
  .init = p
)
warnings()
p3

new_tree <-
tab_tree%>%
  ape::as.phylo()%>%
  tidytree::tree_subset(node=10515)

tab_new_tree <- new_tree%>%
  as_tibble()%>%
  mutate(is_in_our_dataset=ifelse(label%in%data_ordered_seeds$genus,"yes","nope"),
         node_depth= ape::node.depth(new_tree))%>%
  mutate(deep_node_label= ifelse(node_depth>1&!grepl("N[0-9]*$|mrcaott.*|ae$",label),stringr::str_extract(label,"^[a-zA-Z]*"),NA))

pt1 <- ggtree(new_tree)

pt1 %<+% tab_new_tree + 
  ggtree::geom_tippoint(aes(color = is_in_our_dataset))+
  ggtree::  geom_nodelab(
    mapping = aes(
      label = deep_node_label),
    fill = "lightgrey",
    geom = "label")

# get paths that we want to keep
test <-
  tab_new_tree%>%
  filter(is_in_our_dataset=="yes")

path2keep<- NULL
for(nodez in test$node){
  path2genus <- ggtree::get.path(ape::as.phylo(tab_new_tree),10501,nodez) #path from the tip to magnoliophyta
  path2keep <- c(path2keep,path2genus[-length(path2genus)])
}
nodez2keep <- unique(path2keep)

# get paths that we want to collapse
alltips <- tab_new_tree%>%
  filter(node_depth==1)

path2collapse <- NULL
for(tipz in alltips$node){
  path2magnoliophyta <- ggtree::get.path(ape::as.phylo(tab_new_tree),10501,tipz) #path from the tip to magnoliophyta
  path2magnoliophyta <- rev(path2magnoliophyta) #reverse it
  pathpart2remove <- path2magnoliophyta[!path2magnoliophyta%in%nodez2keep] # remove those that we want to keep
  highestnode2remove <- pathpart2remove[length(pathpart2remove)] # keep only highest node
  
  if(highestnode2remove%in%path2collapse| # if alrdy tagged as to remove then don't add it
     highestnode2remove%in%unique(unlist(tidytree::offspring(ape::as.phylo(tab_new_tree),path2collapse)))| # if offspring of a tagged one dont add
     highestnode2remove%in%alltips$node){ # if a tip don't add
    path2collapse <- path2collapse
  }else{
    path2collapse <- c(path2collapse,highestnode2remove) # add it to path2collapse
  }
  print(tipz)
}
length(path2collapse)==length(unique(path2collapse))
nodeztocollapse <- unique(path2collapse)


pt3 <-  purrr::reduce(
  tab_new_tree$node[!tab_new_tree$node%in%nodeztocollapse],
  \(x,y) collapse(x,y),
  .init = pt1
)
warnings()
pt3
