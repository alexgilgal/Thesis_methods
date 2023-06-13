

## We need to create a table that has the TF in a way that TOBIAs likes it
## We neeed that the ensemxg genes match the output of TOBIAS binddetect TF names

TF_ast <- read.delim("data/Astyanax_mexicanus_TF.txt", stringsAsFactors = F) # data taken from here: http://bioinfo.life.hust.edu.cn/static/AnimalTFDB3/download/Astyanax_mexicanus_TF

TF_ast <- TF_ast[,c("Symbol","Ensembl")] # we select only those columns that interest us

orthologs <- read.delim("data/Astyanax_mexicanus_TF_human_orthologs.txt", stringsAsFactors = F) # we load the human orthologs of these TF

colnames(TF_ast)
colnames(orthologs)[1] <- "Ensembl" # we tweak a little the colnames so we can do a proper merge of the 2 tables

TF_ast_with_orth <- merge(TF_ast, orthologs)

TF_list_jaspar <- list.dirs(path = "TOBIAS_snakemake/results_tobias/TFBS/", full.names = F, recursive = F) # this is the total list of jaspar TFBS

## We will use this list of jaspar motifs to generate the table that we need.

TF_final_table <- TF_ast_with_orth

for (i in TF_list_jaspar) {
  
  symbol <- gsub("_.*$", "", i)
  
  finding <- grep(pattern = paste0("^",symbol,"$"), x = TF_final_table$Human.gene.name, ignore.case = T )
  
  if (!is.null(finding)) { 
    
    TF_final_table[finding,"Symbol"] <- i
    
  }
  
}

# finally we will keep only those entries that have the proper format:

TF_final_table_filt <- TF_final_table[grepl("^.*_MA.*$", TF_final_table$Symbol), c(1,2)] 

# TF_final_table_filt$Symbol <- gsub("_MA.*$", "", TF_final_table_filt$Symbol)

write.table(TF_final_table_filt, file = "data/Astyanax_mexicanus_TF_for_TOBIAS.txt", sep = "\t",
            quote = F, col.names = T, row.names = F)

#### we now want to check what Tf have a mutation inside

genes_with_mut <- unique(read.delim("TOBIAS_snakemake/data/genes_with_freebayes_mut_Ens_ids.txt", header = F, stringsAsFactors = F)[,1])

mut_TF <- TF_final_table_filt[TF_final_table_filt$Ensembl %in% genes_with_mut, ]

### Now we will take these genes, and see the downstream TF that they affect

### see:
# https://www.jessesadler.com/post/network-analysis-with-r/
# https://igraph.org/r/
# https://kateto.net/network-visualization

library(igraph)
library(tidyverse)
library(ggraph)

# we need to load the Network tha Tobias generated for us

TF_adjacency <- read.delim("TOBIAS_snakemake/results_tobias/createnetwork_output/All_TF_path_edges.txt", header = T)

# we need to pass this to a proper edge list.

TF_edge_list <- TF_adjacency[,1:2]

# #we separate each Target and each Source so we have only 2 columns and unique genes. But we can have repeated interactions
# 
# for (j in 1:length(TF_adjacency$Source)) {
#   
#   # print(paste('Row', j, 'of',nrow(table)))
#   
#   targets <- unlist(strsplit( TF_adjacency[j,"Targets"], split = ', '))
#   
#   TF_edge_list <- rbind(TF_edge_list, data.frame(Source = TF_adjacency$Source[j], Target = targets))
#   
# }

sum(duplicated(paste(TF_edge_list$Source, TF_edge_list$Target, sep = "_"))) # well aparently there are a lot of duplicates ...

# generating the weights column in the edge list

TF_edge_list <- TF_adjacency %>%  
  group_by(Source, Target) %>%
  summarise(weight = n()) %>% 
  ungroup()

#  we need now to create a node/vertice list eg: ID - Gene

TF_ids <- unique(c(TF_edge_list$Source, TF_edge_list$Target))

TF_node_list <- data.frame(Id = 1:length(TF_ids), label = TF_ids)

TF_node_list$Mutated <- TF_node_list$label %in% mut_TF$Symbol

# Now we will modify the edge list to follow the ID column (eg 1 --> 38 means FOX ---> SOX)

  edges <- TF_edge_list %>% 
  left_join(TF_node_list, by = c("Source" = "label")) %>% 
  rename(from = Id)

edges <- edges %>% 
  left_join(TF_node_list, by = c("Target" = "label")) %>% 
  rename(to = Id)

TF_edge_list_final <- edges %>% select("from", "to", "weight")

write.table(TF_edge_list_final, "results/TF_network_edge_list.txt", col.names = T, row.names = F, quote = F, sep = "\t")

write.table(TF_node_list, "results/TF_network_node_list.txt",col.names = T, row.names = F, quote = F, sep = "\t" )

write.table(TF_node_list[,1:2],"results/TF_network_node_list_simple.txt",col.names = T, row.names = F, quote = F, sep = "\t")



### now we load the data nicely cleaned

TF_network_ig <- graph_from_data_frame(d = TF_edge_list_final, vertices = TF_node_list, directed = TRUE)


# we will delete those edges that have less than a cutoff value of weight

hist(TF_edge_list$weight, breaks = 100)

mean(TF_edge_list$weight)

median(TF_edge_list$weight)

cutoff <- mean(TF_edge_list$weight)

TF_network_ig <- delete_edges(TF_network_ig, E(TF_network_ig)[weight < cutoff])

## We will generate a force directed layout

l <- layout_with_lgl(TF_network_ig,area = 308^5, maxiter = 300)

l <- norm_coords(l, ymin=-1, ymax=1, xmin=-1.5, xmax=1.5)

plot(TF_network_ig, rescale= F, layout = l*1.5,
     edge.arrow.size = 0.2, vertex.shape="none",
     vertex.label.color="steelblue") 



#### some statistics

# density: how many edges there are in a network divided by the total possible number of edges

ecount(TF_network_ig)/(vcount(TF_network_ig)*(vcount(TF_network_ig)-1))

# reciprocity: measure of the likelihood of vertices in a directed network to be mutually linked

reciprocity(TF_network_ig)

# Degree Number of conections of a node:

deg <- degree(TF_network_ig, mode = "all")

plot(TF_network_ig, vertex.size=deg*0.11)

hist(deg, breaks = 100, main="Histogram of node degree")

# Take the top 10 with more degree out: Master regulators

deg_out <- degree(TF_network_ig, mode = "out")

# saving this into the network 

V(TF_network_ig)$deg.out <- deg_out

TF_network_ig

master_regulators <- names(sort(deg_out,decreasing = T)[1:10]) # well, we need the gene name, not the network ID

master_regulators <- TF_node_list[master_regulators,]

# Take the top 10 with more degree in: Master integrators

deg_in <- degree(TF_network_ig, mode = "in")

V(TF_network_ig)$deg.in <- deg_in

master_integrators <- TF_node_list[names(sort(deg_in,decreasing = T)[1:10]),]

# lets make a network only with those, as an example

example_net <- delete_vertices(TF_network_ig, V(TF_network_ig)[!(V(TF_network_ig)  %in% c(master_integrators$Id, master_regulators$Id))])

graph_attr(example_net)

deg_out_ex <- degree(example_net, mode = "out")

l_example <- layout_with_lgl(example_net)

ggraph(example_net, layout = l_example ) + 
  geom_edge_fan(aes(width = weight),color = "gray50", alpha= 0.5, arrow = arrow(length = unit(4, 'mm')),
                 end_cap = circle(3, 'mm'))+ scale_edge_width(range = c(0.2, 1)) +
  geom_node_text(aes(label=label, color = Mutated), size = 5, repel = T) + scale_color_brewer(palette = "Dark2") +
  geom_node_point(aes(size = deg.out)) +
  theme_void()

example_net

# hubs and authorities 

hs <- hub_score(TF_network_ig, weights=NA)$vector

as <- authority_score(TF_network_ig, weights=NA)$vector

plot(TF_network_ig, vertex.size=hs*50, main="Hubs")

plot(TF_network_ig, vertex.size=as*30, main="Authorities")

# comunity detection based on propagating labels

clp <- cluster_label_prop(TF_network_ig)

plot(clp, TF_network_ig)
