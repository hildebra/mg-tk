# script to plot MGS phylogenies in MF and automatically color them
#(c) Klara Cerk

#-------
#Import libraries needed:
#-------
suppressPackageStartupMessages({
library("ape")
library("ggplot2")
library("ggtree")
library("tidytree")
library("phangorn")
library("phytools")
library("dplyr")
library("tidyr")
library("reshape2")
})

#-------
#Read in and prepare tree and metadata:
#-------
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least two arguments must be supplied (.treefile, MGS.matL7.txt).n", call.=FALSE)
}

inMeta = args[1] #MGS.matL7.txt
inTree = args[2] #IQtree_allsites.treefile
outPDF = args[3] #phyloTree.pdf

if(0){
##for manual input:
  inMeta = "IQtree_allsites.treefile"
  inTree = "Abundance/MGS.matL7.txt"
  outPDF = "myOut.pdf"
}


tree <- read.tree(inTree)
metadata_raw <- read.table(inMeta, header=TRUE, row.names=1, sep="\t")


#Prepare and filter metadata as needed for later:
metadata <- colsplit(rownames(metadata_raw), ";", names=c("superkingdom", "phylum", "class", "order", "family", "genus", "species", "MGS")) %>% #reshape2
  filter(grepl("MGS*", MGS)) %>% #filter only taxa for MGS labels
  select("MGS", "superkingdom", "phylum", "class", "order", "family", "genus", "species") %>% #put it in correct order
  as.data.frame() %>% #put it in correct format
  mutate(species = ifelse(species == "?",paste(genus, "unclass") ,species)) #if there is ?, change it to genus unclass


#Lets prepare tree and use ggtree to get info needed:
gg_spTree <- ggtree(tree, layout = 'circular') 
gg_spTree.ann <- gg_spTree %<+% metadata #ggtree


#Re-rooting the tree:
#Identify the Archaea edge: if it's present root the tree based on it, if not use midpoint rooting:
if (length((gg_spTree.ann$data %>% filter(superkingdom == "Archaea"))$label)>0 ) {
  edgeArch <- (gg_spTree.ann$data %>% filter(superkingdom == "Archaea"))$label # identify Archaea edge #tidytree
  print(paste0("New root is based on Archea outgroup ", edgeArch))
  spTree.NewRoot <- root(tree, edgeArch, resolve.root = TRUE, edgelabel=TRUE) # reroot tree at Archaea
  print(paste0("Tree was re-rooted"))
}else {
  spTree.NewRoot <- midpoint.root(tree)# or re-root tree at midpoint
  print(paste0("Tree was re-rooted based on midpoint rooting."))
}


#-------
#All the custom functions needed: #phagorn
#-------
#######
###get_unique_ancestors: returns the list of all unique ancestors of the MRCA nodes in the tree:
#@cladedf - dataframe with 2 columns: clade and node (MRCA of the clade)
#@spTree - phylo tree from which the cladedf was derived (in MFF is IQtree_allsites.treefile)

get_unique_ancestors <- function(cladesdf, spTree) {
  # Initialize an empty list to store the unlisted ancestors
  all_ancestors_combined <- list()
  
  # Loop over each node in clades.df$node to get ancestors for the current node
  for (node in cladesdf$node) {
    ancestors <- Ancestors(spTree, node)
    all_ancestors_combined <- c(all_ancestors_combined, unlist(ancestors))
  }
  
  # Convert the combined list into a single list and make it unique
  all_ancestors_unique <- unique(unlist(all_ancestors_combined))
  # return the unique list of all ancestors
  return(all_ancestors_unique)
}

######
###get_overlaping_node: return the overlaping nodes between ancestors of MRCA tree nodes and all tree nodes, 
#and overlaping clade names:
#@cladedf - dataframe with 2 columns: clade and node (MRCA of the clade)
#@spTree - phylo tree from which the cladedf was derived (in MFF is IQtree_allsites.treefile)

get_overlaping_node <- function(cladesdf, spTree) {
  # Initialize an empty list to store the overlaping nodes
  anb_only <- list()
  # Get ancestors for the current node, and check if they overlap with any other node's ancestors
  for (node in cladesdf$node) {
    anb <- Ancestors(spTree, node)[which(Ancestors(spTree, node) %in% cladesdf$node)]
    anb_only <- c(anb_only, unlist(anb))
  }
  
  all_anb_unique <- unique(unlist(anb_only))  
  # Initialize an empty list to store the overlaping node's clade names
  name_only <- list()
  for (i in all_anb_unique) {
    clade_name <- cladesdf$clade[cladesdf$node == i]
    name_only <- c(name_only, unlist(clade_name))
  }
  
  all_clade_name <- unique(unlist(name_only))
  #return the overlaping nodes and their clade names
  return(list(overlaping_nodes=all_anb_unique,  clade_name=all_clade_name))
}

########
###find_and_check_children: Finds  children of each ancestor of overlapping node and of each MRCA tree node, and checks
#if any children are not in the unique ancestors and if any ancestors children are not in the children's tree:
#@cladedf - dataframe with 2 columns: clade and node (MRCA of the clade)
#@spTree - phylo tree from which the cladedf was derived (in MFF is IQtree_allsites.treefile)
#@node - the overlaping nodes identified with get_overlaping_node function
#@ancestors - all unique ancestors of the MRCA tree nodes identified with get_unique_ancestors

find_and_check_children <- function(node, spTree, ancestors, clades_df) {
  children <- list()
  
  # Find children of each ancestor of overlapping node (node)
  for (ancestor in node) { 
    ancestor_children <- Children(spTree, ancestor)
    children[[ancestor]] <- ancestor_children
  }
  
  children <- unlist(children)
  
  # Find children of each MRCA tree node:
  children_tree <- list()
  for (nodet in ancestors) { 
    ancestort_children <- Children(spTree, nodet)
    children_tree[[nodet]] <- ancestort_children
  }
  
  
  children_tree <- unlist(children_tree)
  unmatched_tree1 <- setdiff(children_tree, ancestors)
  unmatched_tree2 <- setdiff(unmatched_tree1, clades_df$node)
  
  # Check if any children are not in the ancestors
  unmatched_children <- setdiff(children, ancestors) #ancestors
  # Check if any ancestors children are not in the children's tree
  unmatched_children_tree <- union(unmatched_tree2, unmatched_children)
  rest_children <- setdiff(children, unmatched_children)
  
  return(list(children = children, unmatched_children = unmatched_children, rest_children = rest_children, unmatched_children_tree = unmatched_children_tree))
}

#######
###unnest_dataframes: since the output can consist of several lists inside of lists, here is function to change those into dataframes: 
unnest_dataframes <- function(x) {
  y <- do.call(data.frame, x)
  if("data.frame" %in% sapply(y, class)) unnest_dataframes(y)
  y
}


#-------
#Analysis:
#-------

#Detect clade nodes and their names, based on metadata (here we are using phylum):
#Make dataframe for clade nodes
clades.df <- data.frame(
  clade=unique(metadata$phylum),
  node=NA
)
#Find the most recent common ancestor for each clade
for (i in 1:length(clades.df$clade)) {
  
  clades.df$node[i] <- MRCA(spTree.NewRoot, metadata$MGS[metadata$phylum == clades.df$clade[i]])
  
}


#Sometimes the clades.df loop can miss clades that are not monophyletic (Paraphyly/Polyphyly); 
#therefore when you assign colour to them, the colour can overlap the whole tree: 
#Main loop in order to identify overlaping nodes, their names and their new children (which don't overlap):

# Iterate over each node
ancestors <- get_unique_ancestors (clades.df, spTree.NewRoot)
anb <- get_overlaping_node(clades.df, spTree.NewRoot)

if(is.null(anb) == "FALSE"){
  
#Initialize data frame to collect old and new info:
clades.df.old <- data.frame(
  clade = NA,
  node = anb$overlaping_nodes,
  new_nodes = NA,
  extra_nodes = NA)

# Initialize the initial set of children
for (i in 1:length(clades.df.old$node)) {
  clades.df.old$clade[i] <- anb$clade_name[anb$overlaping_nodes == clades.df.old$node[i]]
  rest_children <- clades.df.old$node[i] 
  
  # Initialize an empty list to store children lists
  all_new_nodes <- list()
  all_extra_node <- list()
  
  # Loop until both rest_children have exactly two elements
  while (!(length(rest_children) == 2 && all(length(rest_children) == 2))) {
    # Call the function to find and check children
    newnodes_list <- list()
    extranode_list <- list()
    
    result <- find_and_check_children(rest_children, spTree.NewRoot, ancestors, clades.df)
    newnodes_list <- setdiff(c(newnodes_list, list(result$unmatched_children)), clades.df$node)
    extranode_list <- setdiff(c(extranode_list, list(result$unmatched_children_tree)), newnodes_list)
    
    all_new_nodes <- append(all_new_nodes, newnodes_list)
    all_extra_node <- append(all_extra_node, extranode_list)
    # Retrieve the updated rest_children from the result
    rest_children <- result$rest_children
    
    # If there are no more rest_children, break the loop
    if (length(rest_children) == 0) {
      break
    }}
  
  # Add the current set of children to the list 
  clades.df.old$new_nodes[i] <- list(setdiff(unlist(all_new_nodes), clades.df$node)) 
  clades.df.old$extra_nodes[i] <- list(setdiff(unlist(all_extra_node), clades.df$node))
  
}

clades.df.old$new_nodes <- as.data.frame(do.call(rbind,clades.df.old$new_nodes))


#Reorganise the output, so you can combine it with previous clade.df data:
output_df <- clades.df.old %>%
  select(c(clade, node, new_nodes)) %>%
  unnest_dataframes() %>%
  pivot_longer(!c(clade, node), names_to = "nodes_name", values_to = "nodes") %>% #tidyverse
  select(c(clade, nodes)) %>%
  rename(node = nodes)


#Capture extra nodes that are not part of ancestry(children) of one node - polyphyly:
extra_nodes2 <- setdiff(unique(unlist(clades.df.old$extra_nodes)), unique(unlist(clades.df.old$new_nodes)))

if (length(extra_nodes2)>0){
  
  clades.df.extra <- data.frame(
    clade = NA,
    node = extra_nodes2)
  
  for (node in extra_nodes2) {
    ances.extra <- Ancestors(spTree.NewRoot, node)
    for (i in 1:length(ances.extra)) {
      overlap <- anb$overlaping_nodes == ances.extra[i]
      if(any(overlap) == TRUE){
        clades.df.extra$clade <- anb$clade_name[overlap]
      }
    }
  }
  
  output_df.all <- rbind(output_df, clades.df.extra)
  
} else {
  output_df.all <- output_df #if there is no polyphyly in the tree
}


#Combine old clade.df with new info and update clade.df:
clades.df.new <- rbind(clades.df, output_df.all) 

for (i in 1:length(anb$overlaping_nodes)){
  clades.df.new <- filter(clades.df.new, !node == anb$overlaping_nodes[i])
}

}else {
  clades.df.new <- clades.df
}

#-------
#Figures
#-------
#With all the info, we can color the new tree:
gg_spTree_new<-ggtree(spTree.NewRoot, layout = 'circular')

#Figure0:
gg.tree.new0 <- gg_spTree_new %<+% metadata + 
  geom_tree(aes(color=phylum), size=0.8) +
  geom_tiplab(aes(label=label), size=1.4) +
  xlim(NA, 3.7) +
  theme(legend.position = 'bottom',
        legend.background = element_rect(),
        legend.key = element_blank(), # removes the border
        legend.key.size = unit(0.4, 'cm'), # sets overall area/size of the legend
        legend.text = element_text(size = 6), # text size
        title = element_text(size = 8))

#Figure1:
#Try to add highlights based on previous clade info
gg.tree.new <- gg_spTree_new %<+% metadata +
  geom_highlight(data=clades.df.new, aes(node=node, fill=clade),
                 type = "roundrect",
                 #align="right",
                 #extend=0.1,
                 show.legend=TRUE) +
  xlim(NA, 2.4) +
  geom_tiplab(aes(label=species), size=1.8) +
  theme(legend.position = 'bottom',
        legend.background = element_rect(),
        legend.key = element_blank(), # removes the border
        legend.key.size = unit(0.4, 'cm'), # sets overall area/size of the legend
        legend.text = element_text(size = 5), # text size
        title = element_text(size = 8))


#Figure2:
gg.tree.new2 <- gg_spTree_new %<+% metadata +
  geom_highlight(data=clades.df.new, 
                 aes(node=node, fill=as.factor(clade)),
                 alpha=1,
                 align="right",
                 extend=0.04,
                 show.legend=FALSE) +
  geom_cladelab(data=clades.df.new,
                mapping=aes(node=node, label=clade),
                fontsize=3,
                align="TRUE",
                angle="auto",
                #offset=0.04,
                offset.text=0.28) +
  geom_tree(linewidth=0.3) +
  geom_tippoint() +
  xlim(NA, 5) +
  scale_fill_manual(values=c("#F5F5F5", "#ECECEC", "#C1C1C1", "#FFFFFF", "#F0F0F0", "#CCCCCC", "#c4c4c4", "#FAFAFA","#ebebeb", 
                             "#d4d4d4", "#CECECE", "#EBEBEB" , "#e0e0e0","#c8c8c8",
                             "#e5e5e5", "#DADADA", "#cecece", "#C1C1C1")) 


#-------
#Output
#-------
pdf(outPDF)
plot(gg.tree.new0)
plot(gg.tree.new)
plot(gg.tree.new2)
dev.off()

