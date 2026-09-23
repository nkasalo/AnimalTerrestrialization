
library(dplyr)
library(tidyr)
library(readr)
library(openxlsx)
library(plyr)
library(ontologyIndex)
library(ape)
library(treeio)
library(tidytree)
library(treetools)
library(ggplot2)
library(ggbreak)
library(ggpubr)


##RAW NUMBER INTERSECTION####
cluster_folder = "Clustering"
og_or_not = "noOG"
ps_mapping_folder = "ps_mapping"
terrestrial_nodes = c("Bdelloidea", "Clitellata", "Epiperipatus broadwayi (Onycophora)", "Stylommatophora", "Nematoda", "Tardigrada", "Arachnida", "Myriapoda", "Armadillidium nasatum", "Hexapoda", "Tetrapoda")
gain_or_loss = "gain"


# Load go.obo file with explanations for GO ids
ontology = get_ontology("go.obo", extract_tags = "everything")
GO_table = cbind(ontology$id, ontology$name, ontology$namespace)
GO_table = as.data.frame(GO_table, row.names = NULL)
colnames(GO_table) <- c("GOs", "name", "namespace")
GO_table$Name_ID = paste(GO_table$GOs, GO_table$name)
#Removing technical rows (not connected to biology)
GO_table = GO_table[GO_table$namespace != "external", ]
#Get all non-parent GO terms and remove obsolete
all_parents = unlist(ontology$parents)
terminal_terms = GO_table$GOs[!(GO_table$GOs %in% all_parents)]
#GO_table = subset(GO_table, GOs %in% terminal_terms)
GO_table = GO_table[!grepl("obsolete", GO_table$name), ]


#Load nodes
nodes = read_tsv("nodes.txt", col_names = c("child", "parent"))

#Load names
names = read_tsv("names.txt", col_names = c("taxID", "name"))

#Make a mapping table
name_mapping = read.xlsx("name_mapping.xlsx")
colnames(name_mapping) = c("name","new_label")
name_mapping_full = join(name_mapping, names, by = "name", type = "right")
name_mapping_full$new_label = ifelse(is.na(name_mapping_full$new_label), name_mapping_full$name, name_mapping_full$new_label)


#Get taxIDs of the nodes being compared
nodes_to_compare_taxIDs = names[1:301,]
nodes_to_compare_taxIDs = nodes_to_compare_taxIDs$taxID

#Load the annotations
annotations = read_tsv("nature_all_annotated_default.emapper.annotations")
annotations$query <- sub("\\\\t.*", "", annotations$query)
annotations = annotations %>% dplyr::rename("member" = "query")
annotations$taxID = gsub(".*tx|", "", annotations$member)
annotations$taxID = gsub("|", "", annotations$taxID, fixed = TRUE)


final_result = data.frame()
final_result_terrestrial = data.frame()
for (c in 10:10) {
  print(paste0("c-value: ", c))
  
  #Load the ps mapping
  if (gain_or_loss == "gain") {
    ps_mapping = read_tsv(paste0(ps_mapping_folder, "/gain_c0", c, ".tsv"))
    
    #Subset the ps_mapping file to only the nodes of interest
    ps_mapping = subset(ps_mapping, minTaxID %in% nodes_to_compare_taxIDs)
  }
  
  if(gain_or_loss == "loss") {
    ps_mapping = read_tsv(paste0(ps_mapping_folder, "/loss_c0", c, ".tsv"))
    
    #Subset the ps_mapping file to only the nodes of interest
    ps_mapping = subset(ps_mapping, maxTaxID %in% nodes_to_compare_taxIDs)
  }
  
  
  #Load the clustering file
  current_cluster = read_tsv(paste0(cluster_folder, "/results_0_", c, "/db_clu_all.tsv"), col_names = c("representative", "member"))
  
  #Subset to relevant clusters
  current_cluster = subset(current_cluster, representative %in% ps_mapping$representative)
  
  #Add annotations
  current_cluster = join(current_cluster, annotations, by = "member", type = "left")
  current_cluster = subset(current_cluster, select = c(representative, GOs))
  current_cluster$GOs = current_cluster$GOs %>% replace_na("-")
  
  print("Resolving annotations")
  #Group by cluster representative and combine all annotations of that cluster into a list of unique values
  current_cluster <- current_cluster %>%
    group_by(representative) %>%
    summarise(
      GOs = {
        # Split comma-separated strings into vectors
        vals <- unique(unlist(strsplit(GOs, ",")))
        
        # Remove "-" if any real GO exists
        if (any(vals != "-")) {
          vals <- vals[vals != "-"]
        }
        
        # Collapse back to one string
        paste(vals, collapse = ",")
      },
      .groups = "drop"
    )
  current_cluster = subset(current_cluster, GOs != "-")
  
  #Go through each node and get all of its GO terms
  GO_list = vector("list", length(nodes_to_compare_taxIDs))
  
  print("Creating list of GO terms per node")
  
  position = 1
  for (spec in nodes_to_compare_taxIDs){
    
    #print(paste0(position, " of ", length(nodes_to_compare_taxIDs)))
    
    #Get the name of the node
    spec_name = subset(name_mapping_full, taxID == spec)
    spec_name = spec_name[["new_label"]]
    
    if (gain_or_loss == "gain") {
      subset_mapping = subset(ps_mapping, minTaxID == spec)
    }
    if (gain_or_loss == "loss") {
      subset_mapping = subset(ps_mapping, maxTaxID == spec)
    }
    
    subset_cluster = subset(current_cluster, representative %in% subset_mapping$representative)
    
    current_GOs = unlist(strsplit(unlist(subset_cluster$GOs), ","))
    current_GOs = unique(current_GOs)
    
    #Remove the GOs that are not in the cleaned GO table
    current_GOs = subset(current_GOs, current_GOs %in% c(GO_table$GOs, "-"))
    
    if(is.null(current_GOs)) {
      print("null")
      current_GOs = character(0)
    }
    
    
    #Add to list
    GO_list[[position]] = current_GOs
    names(GO_list)[position] = spec_name
    
    position = position + 1
  }
  
  #Initialize the result vector
  result_vector = c()
  
  #Create a table for saving the results
  result_table = data.frame()
  result_table_terrestrial = data.frame()
  
  print("Calculating intersections for all nodes")
  
  #Do n random samplings of all nodes and calculate their intersection
  n_permutations = 10000
  results = replicate(n_permutations, {
    random_nodes = sample(GO_list, 11)
    
    #intersection within the random set
    intersection = Reduce(intersect, random_nodes)
  })
  
  print("Calculating intersections for terrestrial nodes")
  
  result_vector_terrestrial = c()
  
  GO_list_terrestrial = GO_list[terrestrial_nodes]
  GO_list_terrestrial = GO_list_terrestrial[!is.na(names(GO_list_terrestrial))]
  
  intersection_terrestrial = Reduce(intersect, GO_list_terrestrial)
  
  #Make a result table, with the first row containing the actual value, and other rows permutations
  result_table = as.data.frame(cbind(c, length(intersection_terrestrial)))
  colnames(result_table) = c("cval", "intersection_length")
  
  result_table_all = as.data.frame(cbind(c, results))
  result_table_all$results = sapply(result_table_all$results, length)
  colnames(result_table_all) = c("cval", "intersection_length")
  
  result_table_all = rbind(result_table, result_table_all)

  
  if (nrow(final_result)==0){
    final_result = result_table_all
  } else {
    final_result = rbind(final_result, result_table_all)
  }
  
}

final_result$cval = as.numeric(final_result$cval)

write_tsv(final_result, paste0("permutations_noenrichment_", gain_or_loss, "_", og_or_not, "_", n_permutations, "all11.tsv"))


###Draw permutation results####

data = read_tsv("permutations_noenrichment_gain_noOG_10000all11.tsv")

#data = final_result
n_perm = 10000

obs_raw  <- data$intersection_length[1]
perm_raw <- data$intersection_length[2:(n_perm+1)]

#How many permutations are larger than the observed value
p_value_raw = mean(perm_raw >= obs_raw)
#What is the 5% right bound
average_value = quantile(perm_raw, 0.95)

#Determine the distribution of all the values
breaks = data[2:(n_perm+1),]
breaks = as.data.frame(table(breaks$intersection_length))
breaks$Var2 = round(as.numeric(levels(breaks$Var1))[breaks$Var1], digits = 2)
breaks$probability = breaks$Freq / sum(breaks$Freq)

#Order the breaks from left to right (smallest value of intersection to the highest)
breaks = breaks %>% arrange(Var2)

color_vector = c("black", "red", "black", "black")
plot_raw_permut = ggplot(data = breaks, aes(x = Var2, y=probability*100)) + 
  geom_bar(stat="identity", aes(fill = Var2 >= average_value)) +
  scale_y_break(c(0.3, 7),
                expand = expansion(mult = c(0, 0.1))) +
  scale_fill_manual(values = c("grey", "black")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text=element_text(size=20, color = "black"),
        legend.position = "none",
        axis.title = element_text(size=20),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(color = color_vector)) +
  annotate("rect", xmin = average_value, xmax = Inf, ymin = 0, ymax = Inf, alpha = 0) +
  scale_x_continuous(breaks = c(round(min(breaks$Var2),digits=2), round(obs_raw, digits = 2), round(max(breaks$Var2), digits = 2))) + 
  geom_vline(aes(xintercept = obs_raw), color = "red", linetype = "dashed", size = 1) +
  annotate("text", x = obs_raw - 20, y = 7.2, label = paste0("11 terrestrial\nnodes\n", "p = ", p_value_raw), hjust = 1, size = 7, color = "red") +
  xlab(expression(bold("Number of shared GO terms"))) +
  ylab(expression(bold("Proportion of node groups (%)"))) +
  theme(plot.title = element_text(size = 20)) + 
  scale_y_continuous(expand = expansion(mult = c(0, .1)))


#NODE-WEIGHTED INTERSECTION####

cluster_folder = "Clustering"
og_or_not = "noOG"
ps_mapping_folder = "ps_mapping"
terrestrial_nodes = c("Bdelloidea", "Clitellata", "Epiperipatus broadwayi (Onycophora)", "Stylommatophora", "Nematoda", "Tardigrada", "Arachnida", "Myriapoda", "Armadillidium nasatum", "Hexapoda", "Tetrapoda")
gain_or_loss = "gain"


# Load go.obo file with explanations for GO ids
ontology = get_ontology("go.obo", extract_tags = "everything")
GO_table = cbind(ontology$id, ontology$name, ontology$namespace)
GO_table = as.data.frame(GO_table, row.names = NULL)
colnames(GO_table) <- c("GOs", "name", "namespace")
GO_table$Name_ID = paste(GO_table$GOs, GO_table$name)
#Removing technical rows (not connected to biology)
GO_table = GO_table[GO_table$namespace != "external", ]
#Get all non-parent GO terms and remove obsolete
all_parents = unlist(ontology$parents)
terminal_terms = GO_table$GOs[!(GO_table$GOs %in% all_parents)]
#GO_table = subset(GO_table, GOs %in% terminal_terms)
GO_table = GO_table[!grepl("obsolete", GO_table$name), ]


#Load nodes
nodes = read_tsv("nodes.txt", col_names = c("child", "parent"))

#Load names
names = read_tsv("names.txt", col_names = c("taxID", "name"))

#Make a mapping table
name_mapping = read.xlsx("name_mapping.xlsx")
colnames(name_mapping) = c("name","new_label")
name_mapping_full = join(name_mapping, names, by = "name", type = "right")
name_mapping_full$new_label = ifelse(is.na(name_mapping_full$new_label), name_mapping_full$name, name_mapping_full$new_label)


#Get taxIDs of the nodes being compared
nodes_to_compare_taxIDs = names[1:301,]
nodes_to_compare_taxIDs = nodes_to_compare_taxIDs$taxID

#Load the annotations
annotations = read_tsv("nature_all_annotated_default.emapper.annotations")
annotations$query <- sub("\\\\t.*", "", annotations$query)
annotations = annotations %>% dplyr::rename("member" = "query")
annotations$taxID = gsub(".*tx|", "", annotations$member)
annotations$taxID = gsub("|", "", annotations$taxID, fixed = TRUE)

final_result = data.frame()
final_result_terrestrial = data.frame()
for (c in 10:10) {
  print(paste0("c-value: ", c))
  
  #Load the ps mapping
  if (gain_or_loss == "gain") {
    ps_mapping = read_tsv(paste0(ps_mapping_folder, "/gain_c0", c, ".tsv"))
    
    #Subset the ps_mapping file to only the nodes of interest
    ps_mapping = subset(ps_mapping, minTaxID %in% nodes_to_compare_taxIDs)
  }
  
  if(gain_or_loss == "loss") {
    ps_mapping = read_tsv(paste0(ps_mapping_folder, "/loss_c0", c, ".tsv"))
    
    #Subset the ps_mapping file to only the nodes of interest
    ps_mapping = subset(ps_mapping, maxTaxID %in% nodes_to_compare_taxIDs)
  }
  
  
  #Load the clustering file
  current_cluster = read_tsv(paste0(cluster_folder, "/results_0_", c, "/db_clu_all.tsv"), col_names = c("representative", "member"))
  
  #Subset to relevant clusters
  current_cluster = subset(current_cluster, representative %in% ps_mapping$representative)
  
  #Add annotations
  current_cluster = join(current_cluster, annotations, by = "member", type = "left")
  current_cluster = subset(current_cluster, select = c(representative, GOs))
  current_cluster$GOs = current_cluster$GOs %>% replace_na("-")
  
  print("Resolving annotations")
  #Group by cluster representative and combine all annotations of that cluster into a list of unique values
  current_cluster <- current_cluster %>%
    group_by(representative) %>%
    summarise(
      GOs = {
        # Split comma-separated strings into vectors
        vals <- unique(unlist(strsplit(GOs, ",")))
        
        # Remove "-" if any real GO exists
        if (any(vals != "-")) {
          vals <- vals[vals != "-"]
        }
        
        # Collapse back to one string
        paste(vals, collapse = ",")
      },
      .groups = "drop"
    )
  current_cluster = subset(current_cluster, GOs != "-")
  
  #Go through each node and get all of its GO terms
  GO_list = vector("list", length(nodes_to_compare_taxIDs))
  
  #Create a table in which every node will be assigned a column depicting the relative presence of each GO function
  GO_weighted = t(as.data.frame(unique(GO_table$GOs)))
  colnames(GO_weighted) = c("GOs")
  GO_weighted = as.data.frame(GO_weighted)
  
  print("Creating list of GO terms per node")
  
  position = 1
  for (spec in nodes_to_compare_taxIDs){
    
    print(paste0(position, " of ", length(nodes_to_compare_taxIDs)))
    
    #Get the name of the node
    spec_name = subset(name_mapping_full, taxID == spec)
    spec_name = spec_name[["new_label"]]
    
    
    if (gain_or_loss == "gain") {
      subset_mapping = subset(ps_mapping, minTaxID == spec)
    }
    if (gain_or_loss == "loss") {
      subset_mapping = subset(ps_mapping, maxTaxID == spec)
    }
    
    subset_cluster = subset(current_cluster, representative %in% subset_mapping$representative)
    
    #For weighted counting
    if (nrow(subset_cluster) > 1) {
      current_weighted = unlist(strsplit(unlist(subset_cluster$GOs), ","))
      #length(current_weighted)
      current_weighted = current_weighted[current_weighted %in% GO_table$GOs]
      #length(current_weighted)
      
      
      current_weighted = as.data.frame(table(current_weighted))
      #UNIQUE
      current_weighted$Freq = 1
      #
      
      current_sum = sum(current_weighted$Freq)
      current_weighted$Freq = current_weighted$Freq / current_sum
      
      colnames(current_weighted) = c("GOs", spec_name)
      
      GO_weighted = left_join(GO_weighted, current_weighted, by = "GOs")
    } else {
      current_weighted = cbind(GO_weighted$GOs, NA)
      colnames(current_weighted) = c("GOs", spec_name)
      current_weighted = as.data.frame(current_weighted)
      GO_weighted = left_join(GO_weighted, current_weighted, by = "GOs")
    }
    
    #end weighted counting
    
    
    current_GOs = unlist(strsplit(unlist(subset_cluster$GOs), ","))
    current_GOs = unique(current_GOs)
    
    #Remove the GOs that are not in the cleaned GO table
    current_GOs = subset(current_GOs, current_GOs %in% c(GO_table$GOs, "-"))
    
    if(is.null(current_GOs)) {
      current_GOs = character(0)
    }
    
    
    
    #Add to list
    GO_list[[position]] = current_GOs
    names(GO_list)[position] = spec_name
    
    position = position + 1
  } 
  
  
  
  result_vector = c()
  
  
  #Create a table for saving the results
  result_table = data.frame()
  result_table_terrestrial = data.frame()
  
  print("Calculating intersections for all nodes")
  
  #Do n random samplings of all nodes and calculate their intersection
  n_permutations = 10000
  results = replicate(n_permutations, {
    random_nodes = sample(GO_list, 11)
    
    #Intersection within the random set
    intersection = Reduce(intersect, random_nodes)
    
    #Weighted proportion of functions
    cols = names(random_nodes)
    rows = intersection
    
    current_weighted_intersection = GO_weighted %>%
      select(all_of(c(cols, "GOs"))) %>%
      filter(.data[["GOs"]] %in% rows)
    
    current_weighted_intersection$GOs = NULL
    
    col_sums <- rowMeans(current_weighted_intersection, na.rm = TRUE)
    avg_of_sums <- mean(col_sums)
  })
  
  
  
  print("Calculating intersections for terrestrial nodes")
  
  result_vector_terrestrial = c()
  
  GO_list_terrestrial = GO_list[terrestrial_nodes]
  GO_list_terrestrial = GO_list_terrestrial[!is.na(names(GO_list_terrestrial))]
  
  intersection_terrestrial = Reduce(intersect, GO_list_terrestrial)
  
  #Weighted proportion of functions
  cols = names(GO_list_terrestrial)
  rows = intersection_terrestrial
  
  current_weighted_intersection = GO_weighted %>%
    select(all_of(c(cols, "GOs"))) %>%
    filter(.data[["GOs"]] %in% rows)
  
  current_weighted_intersection$GOs = NULL
  
  col_sums <- rowMeans(current_weighted_intersection, na.rm = TRUE)
  avg_of_sums_terrestrial <- mean(col_sums)
  
  
  
  #Make a result table, with the first row containing the actual value, and other rows permutations
  result_table = as.data.frame(cbind(c, avg_of_sums_terrestrial))
  colnames(result_table) = c("cval", "avg_intersection")
  

  result_table_all = as.data.frame(cbind(c, results))

  colnames(result_table_all) = c("cval", "avg_intersection")
  
  result_table_all = rbind(result_table, result_table_all)
  
  
  if (nrow(final_result)==0){
    final_result = result_table_all
  } else {
    final_result = rbind(final_result, result_table_all)
  }
  
}

final_result$avg_intersection[is.nan(final_result$avg_intersection)] = 0

write_tsv(final_result, paste0("permutations_noenrichment_weighted_", gain_or_loss, "_", og_or_not, "_", n_permutations, "all11.tsv"))


###Draw permutation results####

data = read_tsv("permutations_noenrichment_weighted_gain_noOG_10000all11.tsv")

n_perm = 10000

obs_node  <- data$avg_intersection[1]
perm_node <- data$avg_intersection[2:(n_perm+1)]

#How many permutations are larger than the observed value
p_value_node = mean(perm_node >= obs_node)
#What is the 5% right bound
average_value_node = as.numeric(quantile(perm_node, 0.95))

#Determine the distribution of all the values
breaks_node = data[2:(n_perm+1),]
breaks_node = as.data.frame(table(breaks_node$avg_intersection))
breaks_node$Var2 = round(as.numeric(levels(breaks_node$Var1))[breaks_node$Var1], digits = 2)
breaks_node$probability = breaks_node$Freq / sum(breaks_node$Freq)

#Order the breaks from left to right (smallest value of intersection to the highest)
breaks_node = breaks_node %>% arrange(Var2)

color_vector = c("black", "red", "black", "black")
plot_node = ggplot(data = breaks_node, aes(x = Var2, y=probability*100)) + geom_bar(stat="identity", aes(fill = Var2 >= average_value_node)) +
  scale_fill_manual(values = c("grey", "black")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text=element_text(size=20, color = "black"),
        legend.position = "none",
        axis.title = element_text(size=20),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(color = color_vector)) +
  scale_x_continuous(breaks = c(round(min(breaks_node$Var2), digits = 3), round(obs_node, digits = 3), round(max(breaks_node$Var2), digits = 3))) +
  geom_vline(aes(xintercept = obs_node),
             color = "red",
             linetype = "dashed",
             size = 1) +
  annotate(
    "text",
    x = obs_node - 0.008,
    y = 7,
    label = paste0(
      "11 terrestrial\nnodes\n",
      "p = ", round(p_value_node, digits = 3)
    ),
    hjust = 1,
    size = 7,
    color = "red"
  )+
  xlab(expression(bold("Average frequency of shared GO terms (node-based)"))) +
  ylab(expression(bold("Proportion of node groups (%)"))) +
  theme(plot.title = element_text(size = 20)) + 
  scale_y_continuous(expand = expansion(mult = c(0, .1)))



#CLUSTER-WEIGHTED INTERSECTIONS####

cluster_folder = "Clustering"
og_or_not = "noOG"
ps_mapping_folder = "ps_mapping"
terrestrial_nodes = c("Bdelloidea", "Clitellata", "Epiperipatus broadwayi (Onycophora)", "Stylommatophora", "Nematoda", "Tardigrada", "Arachnida", "Myriapoda", "Armadillidium nasatum", "Hexapoda", "Tetrapoda")
gain_or_loss = "gain"

# Load go.obo file with explanations for GO ids
ontology = get_ontology("go.obo", extract_tags = "everything")
GO_table = cbind(ontology$id, ontology$name, ontology$namespace)
GO_table = as.data.frame(GO_table, row.names = NULL)
colnames(GO_table) <- c("GOs", "name", "namespace")
GO_table$Name_ID = paste(GO_table$GOs, GO_table$name)
#Removing technical rows (not connected to biology)
GO_table = GO_table[GO_table$namespace != "external", ]
#Get all non-parent GO terms and remove obsolete
all_parents = unlist(ontology$parents)
terminal_terms = GO_table$GOs[!(GO_table$GOs %in% all_parents)]
#GO_table = subset(GO_table, GOs %in% terminal_terms)
GO_table = GO_table[!grepl("obsolete", GO_table$name), ]


#Load nodes
nodes = read_tsv("nodes.txt", col_names = c("child", "parent"))

#Load names
names = read_tsv("names.txt", col_names = c("taxID", "name"))

#Make a mapping table
name_mapping = read.xlsx("name_mapping.xlsx")
colnames(name_mapping) = c("name","new_label")
name_mapping_full = join(name_mapping, names, by = "name", type = "right")
name_mapping_full$new_label = ifelse(is.na(name_mapping_full$new_label), name_mapping_full$name, name_mapping_full$new_label)


#Get taxIDs of the nodes being compared
nodes_to_compare_taxIDs = names[1:301,]
nodes_to_compare_taxIDs = nodes_to_compare_taxIDs$taxID

#Load the annotations
annotations = read_tsv("nature_all_annotated_default.emapper.annotations")
annotations$query <- sub("\\\\t.*", "", annotations$query)
annotations = annotations %>% dplyr::rename("member" = "query")
annotations$taxID = gsub(".*tx|", "", annotations$member)
annotations$taxID = gsub("|", "", annotations$taxID, fixed = TRUE)


final_result = data.frame()
final_result_terrestrial = data.frame()
for (c in 10:10) {
  print(paste0("c-value: ", c))
  
  #Load the ps mapping
  if (gain_or_loss == "gain") {
    ps_mapping = read_tsv(paste0(ps_mapping_folder, "/gain_c0", c, ".tsv"))
    
    #Subset the ps_mapping file to only the nodes of interest
    ps_mapping = subset(ps_mapping, minTaxID %in% nodes_to_compare_taxIDs)
  }
  
  if(gain_or_loss == "loss") {
    ps_mapping = read_tsv(paste0(ps_mapping_folder, "/loss_c0", c, ".tsv"))
    
    #Subset the ps_mapping file to only the nodes of interest
    ps_mapping = subset(ps_mapping, maxTaxID %in% nodes_to_compare_taxIDs)
  }
  
  
  #Load the clustering file
  current_cluster = read_tsv(paste0(cluster_folder, "/results_0_", c, "/db_clu_all.tsv"), col_names = c("representative", "member"))
  
  #Subset to relevant clusters
  current_cluster = subset(current_cluster, representative %in% ps_mapping$representative)
  
  #Add annotations
  current_cluster = join(current_cluster, annotations, by = "member", type = "left")
  current_cluster = subset(current_cluster, select = c(representative, GOs))
  current_cluster$GOs = current_cluster$GOs %>% replace_na("-")
  
  print("Resolving annotations")
  #For each cluster, calculate the frequency of each function
  current_cluster_weights = current_cluster %>%
    mutate(GOs = strsplit(GOs, ",")) %>%
    unnest(GOs) %>%
    mutate(GOs = trimws(GOs)) %>%
    filter(GOs %in% GO_table$GOs) %>%
    group_by(representative, GOs) %>%
    summarise(n_GO = 1, .groups = "drop") %>% # or n_GO = n() for within cluster counts
    group_by(representative) %>%
    mutate(
      total_GO = sum(n_GO),
      fraction = n_GO / total_GO
    ) %>%
    ungroup()
  
  
  #Group by cluster representative and combine all annotations of that cluster into a list of unique values
  current_cluster <- current_cluster %>%
    group_by(representative) %>%
    summarise(
      GOs = {
        # Split comma-separated strings into vectors
        vals <- unique(unlist(strsplit(GOs, ",")))
        
        # Remove "-" if any real GO exists
        if (any(vals != "-")) {
          vals <- vals[vals != "-"]
        }
        
        # Collapse back to one string
        paste(vals, collapse = ",")
      },
      .groups = "drop"
    )
  current_cluster = subset(current_cluster, GOs != "-")
  
  #Go through each node and get all of its GO terms
  GO_list = vector("list", length(nodes_to_compare_taxIDs))
  
  #Create a table in which every node will be assigned a column depicting the relative presence of each GO function
  GO_weighted = t(as.data.frame(unique(GO_table$GOs)))
  colnames(GO_weighted) = c("GOs")
  GO_weighted = as.data.frame(GO_weighted)
  
  print("Creating list of GO terms per node")
  
  position = 1
  for (spec in nodes_to_compare_taxIDs){
    
    print(paste0(position, " of ", length(nodes_to_compare_taxIDs)))
    
    #Get the name of the node
    spec_name = subset(name_mapping_full, taxID == spec)
    spec_name = spec_name[["new_label"]]
    
    
    if (gain_or_loss == "gain") {
      subset_mapping = subset(ps_mapping, minTaxID == spec)
    }
    if (gain_or_loss == "loss") {
      subset_mapping = subset(ps_mapping, maxTaxID == spec)
    }
    
    subset_cluster = subset(current_cluster, representative %in% subset_mapping$representative)
    
    #For weighted counting
    if (nrow(subset_cluster) > 1) {
      current_weighted = subset(current_cluster_weights, representative %in% subset_mapping$representative)
      
      current_weighted = current_weighted %>% group_by(GOs) %>% summarise(Freq = mean(fraction))
      
      colnames(current_weighted) = c("GOs", spec_name)
      
      GO_weighted = left_join(GO_weighted, current_weighted, by = "GOs")
    } else {
      current_weighted = cbind(GO_weighted$GOs, NA)
      colnames(current_weighted) = c("GOs", spec_name)
      current_weighted = as.data.frame(current_weighted)
      GO_weighted = left_join(GO_weighted, current_weighted, by = "GOs")
    }
    
    #end weighted counting
    
    
    current_GOs = unlist(strsplit(unlist(subset_cluster$GOs), ","))
    current_GOs = unique(current_GOs)
    
    #Remove the GOs that are not in the cleaned GO table
    current_GOs = subset(current_GOs, current_GOs %in% c(GO_table$GOs, "-"))
    
    if(is.null(current_GOs)) {
      current_GOs = character(0)
    }
    
    
    
    #Add to list
    GO_list[[position]] = current_GOs
    names(GO_list)[position] = spec_name
    
    position = position + 1
  } 
  
  
  
  result_vector = c()
  
  
  #Create a table for saving the results
  result_table = data.frame()
  result_table_terrestrial = data.frame()
  
  print("Calculating intersections for all nodes")
  
  #Do n random samplings of all nodes and calculate their intersection
  n_permutations = 10000
  results = replicate(n_permutations, {
    random_nodes = sample(GO_list, 11)
    
    #Intersection size within the random set
    intersection = Reduce(intersect, random_nodes)
    
    #Weighted proportion of functions
    cols = names(random_nodes)
    rows = intersection
    
    current_weighted_intersection = GO_weighted %>%
      select(all_of(c(cols, "GOs"))) %>%
      filter(.data[["GOs"]] %in% rows)
    
    current_weighted_intersection$GOs = NULL
    
    col_sums <- rowMeans(current_weighted_intersection, na.rm = TRUE)
    avg_of_avg <- mean(col_sums)
  })
  
  
  
  print("Calculating intersections for terrestrial nodes")
  
  result_vector_terrestrial = c()
  
  GO_list_terrestrial = GO_list[terrestrial_nodes]
  GO_list_terrestrial = GO_list_terrestrial[!is.na(names(GO_list_terrestrial))]
  
  intersection_terrestrial = Reduce(intersect, GO_list_terrestrial)
  
  #Weighted proportion of functions
  cols = names(GO_list_terrestrial)
  rows = intersection_terrestrial
  
  current_weighted_intersection = GO_weighted %>%
    select(all_of(c(cols, "GOs"))) %>%
    filter(.data[["GOs"]] %in% rows)
  
  current_weighted_intersection$GOs = NULL
  
  col_sums <- rowMeans(current_weighted_intersection, na.rm = TRUE)
  avg_of_sums_terrestrial <- mean(col_sums)
  
  
  
  #Make a result table, with the first row containing the actual value, and other rows permutations
  result_table = as.data.frame(cbind(c, avg_of_sums_terrestrial))
  colnames(result_table) = c("cval", "avg_intersection")
  
  result_table_all = as.data.frame(cbind(c, results))
  colnames(result_table_all) = c("cval", "avg_intersection")
  
  result_table_all = rbind(result_table, result_table_all)

  
  if (nrow(final_result)==0){
    final_result = result_table_all
  } else {
    final_result = rbind(final_result, result_table_all)
  }
  
}


final_result$avg_intersection[is.nan(final_result$avg_intersection)]<-0


write_tsv(final_result, paste0("permutations_noenrichment_weightedCluster_", gain_or_loss, "_", og_or_not, "_", n_permutations, "all11.tsv"))


###Draw permutation results####

library(ggbreak)

data = read_tsv("permutations_noenrichment_weightedCluster_gain_noOG_10000all11.tsv")

#data = final_result
n_perm = 10000

obs_clu  <- data$avg_intersection[1]
perm_clu <- data$avg_intersection[2:(n_perm+1)]

#How many permutations are larger than the observed value
p_value_clu = mean(perm_clu >= obs_clu)
#What is the 5% right bound
average_value_clu = as.numeric(quantile(perm_clu, 0.95))

#Determine the distribution of all the values
breaks_clu = data[2:(n_perm+1),]
breaks_clu = as.data.frame(table(breaks_clu$avg_intersection))
breaks_clu$Var2 = round(as.numeric(levels(breaks_clu$Var1))[breaks_clu$Var1], digits = 3)
breaks_clu$probability = breaks_clu$Freq / sum(breaks_clu$Freq)

#Order the breaks from left to right (smallest value of intersection to the highest)
breaks_clu = breaks_clu %>% arrange(Var2)


color_vector = c("black", "red", "black", "black")
plot_cluster = ggplot(data = breaks_clu, aes(x = Var2, y=probability*100)) + geom_bar(stat="identity", aes(fill = Var2 >= average_value_clu)) +
  scale_fill_manual(values = c("grey", "black")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text=element_text(size=20, color = "black"),
        legend.position = "none",
        axis.title = element_text(size=20),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(color = color_vector)) +
  scale_x_continuous(breaks = c(round(min(breaks_clu$Var2), digits = 3), round(obs_clu, digits = 4), round(max(breaks_clu$Var2), digits = 4))) + 
  geom_vline(aes(xintercept = obs_clu),
             color = "red",
             linetype = "dashed",
             size = 1) +
  annotate(
    "text",
    x = obs_clu - 0.00005,
    y = 7.2,
    label = paste0(
      "11 terrestrial\nnodes\n",
      "p = ", round(p_value_clu, digits = 3)
    ),
    hjust = 1,
    size = 7,
    color = "red"
  ) +
  xlab(expression(bold("Average frequency of shared GO terms (cluster-based)"))) +
  ylab(expression(bold("Proportion of node groups (%)"))) +
  theme(plot.title = element_text(size = 20)) + 
  scale_y_continuous(expand = expansion(mult = c(0, .1)))





#GAIN RATE ANALYSES####

##no enrichment####
cluster_folder = "Clustering"
og_or_not = "noOG"
ps_mapping_folder = "ps_mapping"
terrestrial_nodes = c("Bdelloidea", "Clitellata", "Epiperipatus broadwayi (Onycophora)", "Stylommatophora", "Nematoda", "Tardigrada", "Arachnida", "Myriapoda", "Armadillidium nasatum", "Hexapoda", "Tetrapoda")
aquatic_nodes = c("node136", "node145", "node81", "node18", "node143", "node149", "node30", "node88", "node52", "node61", "node78")
gain_or_loss = "gain"

#Load the results of gene gain/loss mapping
results = read_tsv("Phylo_v1_allc_output.tsv")



#Load nodes
nodes = read_tsv("nodes.txt", col_names = c("child", "parent"))

#Load names
names = read_tsv("names.txt", col_names = c("taxID", "name"))

#Make a mapping table
name_mapping = read.xlsx("name_mapping.xlsx")
colnames(name_mapping) = c("name","new_label")
name_mapping_full = join(name_mapping, names, by = "name", type = "right")
name_mapping_full$new_label = ifelse(is.na(name_mapping_full$new_label), name_mapping_full$name, name_mapping_full$new_label)


#Get taxIDs of the nodes being compared
nodes_to_compare_taxIDs = names[1:301,]
nodes_to_compare_taxIDs = nodes_to_compare_taxIDs$taxID

#Prepare a list of all terrestrial nodes - descendants of the known terrestrial nodes
terrestrial_taxIDs = subset(name_mapping_full, new_label %in% terrestrial_nodes)

#Load the tree
tree = read.phyloxml("tree_v2.xml")
time_reference = read.tree("tree_calib.txt")
# Convert tree to tabular format
tree_tbl = as_tibble(tree)
#Go through terrestrial nodes and make a list of all descendants of all terrestrial nodes
all_terrestrial = c()
for (cur_id in terrestrial_taxIDs$taxID) {
  # Get all descendants of a node
  current_node = subset(tree_tbl, id == cur_id)
  current_node = current_node$node
  descendants = offspring(tree_tbl, current_node, self_include = TRUE)
  descendants = as.numeric(descendants$id)
  all_terrestrial = c(all_terrestrial, descendants)
  all_terrestrial = unique(all_terrestrial)
}

#Load the trees and get branch lengths
all_labels = as.data.frame(tidytree::as_tibble(time_reference))
all_labels_raw = as.data.frame(tidytree::as_tibble(tree))
all_labels_raw$id = as.numeric(all_labels_raw$id)

all_labels = subset(all_labels, select = c(node, branch.length))
colnames(all_labels) = c("id", "branch.length")

calibrations = left_join(all_labels_raw, all_labels, by = "id")
calibrations = subset(calibrations, select = c(label, branch.length))
colnames(calibrations) = c("ps_name", "branch_length")


final_result = data.frame()
final_result_terrestrial = data.frame()
for (c in 10:10) {
  print(paste0("c-value: ", c))
  
  #Subset to the current c-value
  current_results = subset(results, cval == c)
  
  #Subset to gain or loss
  if (gain_or_loss == "gain") {
    current_results = subset(current_results, select = c(codeGained, ps_name_gained, gained_clusters))
  } else {
    current_results = subset(current_results, select = c(codeLost, ps_name_lost, lost_clusters))
  }
  colnames(current_results) = c("ps_taxID", "ps_name", "cluster_count")
  
  #Remove repeated phylostrata
  current_results = na.omit(current_results)
  current_results = unique(current_results)
  #Remove root and outgroups nodes
  current_results = subset(current_results, ps_name != "outgroups")
  current_results = subset(current_results, ps_name != "root")
  
  #Add calibrations
  current_results = left_join(current_results, calibrations, by = "ps_name")
  current_results = na.omit(current_results)
  
  #Add updated node names
  current_results = left_join(current_results, name_mapping_full, by = c("ps_name" = "name"))
  current_results$taxID = NULL
  
  #Multiply branch length by 100 to make the number represent millions of years
  current_results$branch_length = current_results$branch_length * 100
  
  #Calculate the cluster gain rate for each node
  current_results$individual_gain_rate = current_results$cluster_count / current_results$branch_length
  
  #Initialize the result vector
  result_vector = c()
  
  #Create a table for saving the results
  result_table = data.frame()
  result_table_terrestrial = data.frame()
  
  print("Calculating permutations")
  
  #Calculate the rate for the observed set
  terrestrial_subset = subset(current_results, new_label %in% terrestrial_nodes)
  terrestrial_subset = mean(terrestrial_subset$individual_gain_rate)
  
  terrestrial_subset_cumulative = subset(current_results, new_label %in% terrestrial_nodes)
  terrestrial_subset_cumulative = sum(terrestrial_subset_cumulative$cluster_count) / sum(terrestrial_subset_cumulative$branch_length)
  
  
  #Do n random samplings of all nodes and calculate their intersection
  n_permutations = 10000
  result_permutations = replicate(n_permutations, {
    random_nodes = current_results[sample(nrow(current_results), 11), ]
    
    # Average pairwise intersection size within the random set
    mean_rate = mean(random_nodes$individual_gain_rate)
  })
  
  result_permutations_cumulative = replicate(n_permutations, {
    random_nodes = current_results[sample(nrow(current_results), 11), ]
    
    # Average pairwise intersection size within the random set
    cumulative_gain = sum(random_nodes$cluster_count)
    cumulative_divergence = sum(random_nodes$branch_length)
    total = cumulative_gain/cumulative_divergence
  })
  
  
  
  #Calculate permutations without replacement against all aquatic nodes only
  current_results = subset(current_results, !(ps_taxID %in% all_terrestrial))
  
  n_permutations = 10000
  result_permutations_all_aquatic = replicate(n_permutations, {
    random_nodes = current_results[sample(nrow(current_results), 11), ]
    
    # Average pairwise intersection size within the random set
    mean_rate = mean(random_nodes$individual_gain_rate)
  })
  
  result_permutations_cumulative_all_aquatic = replicate(n_permutations, {
    random_nodes = current_results[sample(nrow(current_results), 11), ]
    
    # Average pairwise intersection size within the random set
    cumulative_gain = sum(random_nodes$cluster_count)
    cumulative_divergence = sum(random_nodes$branch_length)
    total = cumulative_gain/cumulative_divergence
  })
  
  
  
  #Now calculate permutations with replacement against 11 aquatic nodes only
  current_results = subset(current_results, new_label %in% aquatic_nodes)
  
  n_permutations = 10000
  result_permutations_aquatic = replicate(n_permutations, {
    random_nodes = current_results[sample(nrow(current_results), 11, replace = TRUE), ]
    
    # Average pairwise intersection size within the random set
    mean_rate = mean(random_nodes$individual_gain_rate)
  })
  
  result_permutations_cumulative_aquatic = replicate(n_permutations, {
    random_nodes = current_results[sample(nrow(current_results), 11, replace = TRUE), ]
    
    # Average pairwise intersection size within the random set
    cumulative_gain = sum(random_nodes$cluster_count)
    cumulative_divergence = sum(random_nodes$branch_length)
    total = cumulative_gain/cumulative_divergence
  })
  
  
  #Make a result table, with the first row containing the actual value, and other rows permutations
  result_table = as.data.frame(cbind(c, terrestrial_subset, terrestrial_subset_cumulative, terrestrial_subset,
                                     terrestrial_subset_cumulative, terrestrial_subset, terrestrial_subset_cumulative))
  colnames(result_table) = c("cval", "gain_rate", "cumulative_gain_rate", "gain_rate_aquatic", "cumulative_gain_rate_aquatic", "gain_rate_all_aquatic", "cumulative_gain_rate_all_aquatic")
  
  result_table_all = as.data.frame(cbind(c, result_permutations, result_permutations_cumulative, result_permutations_aquatic,
                                         result_permutations_cumulative_aquatic, result_permutations_all_aquatic, result_permutations_cumulative_all_aquatic))
  colnames(result_table_all) = c("cval", "gain_rate", "cumulative_gain_rate", "gain_rate_aquatic", "cumulative_gain_rate_aquatic", "gain_rate_all_aquatic", "cumulative_gain_rate_all_aquatic")
  
  result_table_all = rbind(result_table, result_table_all)
  
  
  if (nrow(final_result)==0){
    final_result = result_table_all
  } else {
    final_result = rbind(final_result, result_table_all)
  }
  
}

write_tsv(final_result, paste0("permutations_gainrate_noenrichment_", gain_or_loss, "_", og_or_not, "_", n_permutations, ".tsv"))


###Draw permutation results####

data = read_tsv("permutations_gainrate_noenrichment_gain_noOG_10000.tsv")
n_perm = 10000

#Go through all c-values
pdf("permutations_gain_rate_v3.pdf", width = 17, height = 8)
for (c in 10:10) {
  print = c
  
  current_data = subset(data, cval == c)
  
  plot1 = local({
    
    obs  <- current_data$gain_rate[1]
    perm <- current_data$gain_rate[2:(n_perm+1)]
    
    #How many permutations are larger than the observed value
    p_value = mean(perm >= obs)
    #What is the 5% right bound
    average_value = as.numeric(quantile(perm, 0.95))
    
    #Determine the distribution of all the values
    breaks1 = current_data[2:(n_perm+1),]
    breaks1 = as.data.frame(table(breaks1$gain_rate))
    breaks1$Var2 = round(as.numeric(levels(breaks1$Var1))[breaks1$Var1], digits = 1)
    breaks1$probability = breaks1$Freq / sum(breaks1$Freq)
    
    
    #Order the breaks from left to right (smallest value of intersection to the highest)
    breaks1 = breaks1 %>% arrange(Var2)
    
    
    
    #pdf("Figures_rev/11_paths_PDF.pdf", width = 9, height = 7)
    color_vector = c("black", "red", "black", "black")
    ggplot(data = breaks1, aes(x = Var2, y=probability*100)) + geom_bar(stat="identity", aes(fill = Var2 >= round(average_value, digits = 2))) +
      scale_fill_manual(values = c("grey", "black")) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            axis.text=element_text(size=20, color = "black"),
            legend.position = "none",
            axis.title = element_text(size=20),
            plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(color = color_vector)) +
      scale_x_continuous(
        breaks = c(round(min(breaks1$Var2), digits = 0),
                   round(obs, digits = 0),
                   round(max(breaks1$Var2), digits = 0))
      ) +
      annotate(
        "text",
        x = obs + 3,
        y = 1.3,
        hjust = 0,
        size = 7,
        color = "red",
        label = paste0(
          "11 terrestrial\nnodes\n",
          "p = ", round(p_value, digits = 3)
        )
      ) +
      geom_vline(aes(xintercept = obs),
                 color = "red",
                 linetype = "dashed",
                 size = 1) +
      xlab(expression(bold("Average number of gene clusters gained per MY"))) +
      ylab(expression(bold("Proportion of permutations (%)"))) +
      ggtitle(paste0("Permutations of all 301 nodes")) +
      theme(plot.title = element_text(size = 20)) + 
      scale_y_continuous(expand = expansion(mult = c(0, .1)))

  })
  
  
  #GAIN RATE ALL AQUATIC
  
  plot5 = local({
    
    obs  <- current_data$gain_rate_all_aquatic[1]
    perm <- current_data$gain_rate_all_aquatic[2:(n_perm+1)]
    
    #How many permutations are larger than the observed value
    p_value = mean(perm >= obs)
    #What is the 5% right bound
    average_value = as.numeric(quantile(perm, 0.95))
    
    #Determine the distribution of all the values
    breaks3 = current_data[2:(n_perm+1),]
    breaks3 = as.data.frame(table(breaks3$gain_rate_all_aquatic))
    breaks3$Var2 = round(as.numeric(levels(breaks3$Var1))[breaks3$Var1], digits = 1)
    breaks3$probability = breaks3$Freq / sum(breaks3$Freq)
    
    #Order the breaks from left to right (smallest value of intersection to the highest)
    breaks3 = breaks3 %>% arrange(Var2)
    

    
    
    color_vector = c("black", "red", "black", "black")
    ggplot(data = breaks3, aes(x = Var2, y=probability*100)) + geom_bar(stat="identity", aes(fill = Var2 >= round(average_value, digits = 2))) +
      scale_fill_manual(values = c("grey", "black")) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            axis.text=element_text(size=20, color = "black"),
            legend.position = "none",
            axis.title = element_text(size=20),
            plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(color = color_vector)) +
      scale_x_continuous(
        breaks = c(round(min(breaks3$Var2), digits = 0),
                   round(obs, digits = 0),
                   round(max(breaks3$Var2), digits = 0))
      ) +
      annotate(
        "text",
        x = obs + 2,
        y = 0.8,
        hjust = 0,
        size = 7,
        color = "red",
        label = paste0(
          "11 terrestrial\nnodes\n",
          "p = ", round(p_value, digits = 3)
        )
      ) +
      geom_vline(aes(xintercept = obs),
                 color = "red",
                 linetype = "dashed",
                 size = 1) +
      xlab(expression(bold("Average number of gene clusters gained per MY"))) +
      ylab(expression(bold("Proportion of permutations (%)"))) +
      ggtitle(paste0("Permutations of all 162 aquatic nodes")) +
      theme(plot.title = element_text(size = 20)) + 
      scale_y_continuous(expand = expansion(mult = c(0, .1)))
    
    
  })
  
  
  plot = ggarrange(plot1, plot5, nrow = 1, ncol = 2, labels = paste0("",letters[1:2],""), font.label = list(size = 24))
  
  print(plot)
  
}

dev.off()

