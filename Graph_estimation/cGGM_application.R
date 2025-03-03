setwd("/group/diangelantonio/users/alessia_mapelli/Prot_graphs/UKB_data/APP_82779/4.Graph_estimation")
rm(list=ls(all=TRUE))

library(pak)
library(huge)
library(glasso)
library(glmnet)
library(igraph)
library(Matrix)


source("/group/diangelantonio/users/alessia_mapelli/Prot_graphs/UKB_data/APP_82779/3.Methodology_implementation/Network_estimation_functions_parallel.R")

# PREPARE THE DATA
load("/group/diangelantonio/users/alessia_mapelli/Prot_graphs/UKB_data/APP_82779/1.Data_preprocessing/3.Final_version_of_dataset_for_the_analysis_on_diabetes/complete_DIAB_analysis_82779_20240716.RData")
prot_names <- colnames(prot_data)[-1]

prot_conv_path <- "/exchange/healthds/ukb-ppp/uniprot_panel.rtf"
prot_coversion <- read.csv(prot_conv_path, sep = "\t", header=T)
prot_coversion <- prot_coversion[prot_coversion$Assay %in% prot_names,]
table(prot_coversion$Panel)

Cardiometabolic_prot <- prot_coversion[prot_coversion$Panel == "Cardiometabolic",]

# Cohort selection
table(data$outcome_T2D_diagnosed.HbA1c.drug, useNA="ifany")
table(data$E10)

study_cohort <- data[,c("id","sex","Age_at_recruitment","outcome_T2D_diagnosed.HbA1c.drug", Cardiometabolic_prot$Assay)]
colnames(study_cohort)[2:4] <- c( "Sex", "Age", "Diab")
study_cohort$Diab <- factor(study_cohort$Diab, levels = c("Non diabetic","Incident"))
str(study_cohort)
colSums(is.na(study_cohort))

library(caret)
library(dplyr)
library(tidyr)

set.seed(123)

# Step 1: Split the dataset into 80% and 20% maintaining the proportion of Diab and Sex
index_1 <- createDataPartition(study_cohort$Diab, p = 0.8, list = FALSE, times = 1)
selection_data <- study_cohort[index_1, ]
remaining_data <- study_cohort[-index_1, ]

table(selection_data$Diab, selection_data$Sex)
table(study_cohort$Diab, study_cohort$Sex)


# Step 2: From the remaining 20%, divide into two parts matching the training set on sex and age
incident_participants <- remaining_data %>% filter(Diab == "Incident")
non_diabetic_participants <- remaining_data %>% filter(Diab == "Non diabetic")

set.seed(123)
index_2 <- createDataPartition(incident_participants$Sex, p = 0.5, list = FALSE, times = 1)
incident_training_data <- incident_participants[index_2, ]
incident_testing_data <- incident_participants[-index_2, ]
train_data <- data.frame()
for (i in 1:nrow(incident_training_data)) {
  incident_row <- incident_training_data[i, ]
  set.seed(123)
  matching_non_diabetic <- non_diabetic_participants %>%
    filter(Age == incident_row$Age, Sex == incident_row$Sex) %>%
    sample_n(3, replace = FALSE)
  train_data <- rbind(train_data, incident_row, matching_non_diabetic)
  non_diabetic_participants <- anti_join(non_diabetic_participants, matching_non_diabetic)
}

test_data <- anti_join(remaining_data, train_data)

prop.table(table(train_data$Diab))
table(test_data$Diab)

head(selection_data)


# CHECK
write.table(selection_data$id, "Selection_data_id.txt", col.names = F)
write.table(train_data$id, "Train_data_id.txt",col.names = F)
write.table(test_data$id, "Test_data_id.txt",col.names = F)

# Load the prior

STRING_prior <- readMM("/group/diangelantonio/users/alessia_mapelli/Prot_graphs/UKB_data/APP_82779/1.Data_preprocessing/String_preprocessing/STRING_W")
colnames_W <- read.table("/group/diangelantonio/users/alessia_mapelli/Prot_graphs/UKB_data/APP_82779/1.Data_preprocessing/String_preprocessing/STRING_W_colnames.txt")$V1
rownames_W <- read.table("/group/diangelantonio/users/alessia_mapelli/Prot_graphs/UKB_data/APP_82779/1.Data_preprocessing/String_preprocessing/STRING_W_rownames.txt")$V1

dimnames(STRING_prior)<- list(rownames_W,colnames_W)

STRING_prior <- STRING_prior[Cardiometabolic_prot$Assay,Cardiometabolic_prot$Assay] 

sum(colnames(selection_data)[-(1:4)] == colnames(STRING_prior))


# 
# conditional GGMReg - mean regressed and covariance over the covriate
#
res_4 <- GGReg_full_estimation (x=selection_data[,-(1:4)], # expression data, nrow: subjects; ncol: features.
                                known_ppi = STRING_prior, # previously known PPI 
                                covariates = selection_data[,(2:4)],       #covariates to regress on
                                scr = TRUE, # optional screening to speed up
                                gamma = NULL, #Person correlation threshold for the screening
                                lambda_mean = NULL, #Penalization term in the Lasso for estimation of the mean
                                lambda_mean_type = "1se", # or "min"
                                lambda_prec = NULL, #Penalization term in Lasso for the estimation of the precision matrix
                                lambda_prec_type = "1se", # or "min"
                                weight = 1.1, #Multiplicative factor for the penalization term when prior knowledgw is not available
                                asparse = 0.75, # The relative weight to put on the `1-norm in sparse group lasso
                                verbose = TRUE,
                                eps = 1e-08)

save(res_4, file = "cGGMLReg_Cardiometabolic_17072024.RData")

# load("Graph estimation/cGGMReg_Cardiometabolic_06052024.RData")
# 
# # Describe the results
# # https://kateto.net/netscix2016.html
# # https://r.igraph.org/articles/igraph.html
# D <- res_4$Dic_Delta_hat$DiabPrevalent
# 
# outcome <- as.matrix(D)
# str(outcome)
# g1 <- graph_from_adjacency_matrix(outcome, mode = "max", diag = F, weighted=T,add.colnames = NA, add.rownames = NULL )
# g1
# # layout=layout_with_lgl(g1)
# # layout=layout_with_drl(g1)
# layout = layout_nicely(g1)
# dev.off()
# plot(g1,layout= layout, edge.arrow.size=1, vertex.color="lightskyblue1", vertex.size=10, 
#      
#      vertex.frame.color="lightskyblue1", vertex.label.color="black", 
#      
#      vertex.label.cex=0.5, vertex.label.dist=0, edge.curved=0.2, edge.color= "gray", edge.width =0.1, margin=c(0,0,0,0)) 
# V(g1)
# # 366 nodes
# E(g1)
# # 57693
# 57693/(366*366)
# # 0.4306862
# g1[]
# E(g1)$weight
# edge_attr(g1)
# plot(density(E(g1)$weight))
# vertex_attr(g1)
# 
# diab_eff <- res_4$Cov_effect[,"DiabPrevalent"]
# names(diab_eff) <- colnames(res_4$z)
# head(diab_eff)
# 
# diab_eff <- diab_eff[diab_eff !=0]
# 
# prot_diab_eff_mean <- names(diab_eff)
# 
# # Color the nodes differemtly based on if their mean is influenced by diab
# vertex_colors <- rep("lightskyblue1", vcount(g1))  
# vertex_colors[V(g1)$name %in% prot_diab_eff_mean] <- "darkturquoise"
# 
# plot(g1,
#      layout = layout, 
#      edge.arrow.size = 1, 
#      vertex.color = vertex_colors,  # Apply the custom vertex colors
#      vertex.size = 10,
#      vertex.frame.color = vertex_colors, 
#      vertex.label.color = "black",
#      vertex.label.cex = 0.5,
#      vertex.label.dist = 0,
#      edge.curved = 0.2,
#      edge.color = "gray",
#      edge.width = 0.1,
#      margin = c(0, 0, 0, 0))
# 
# g_abs <- g1
# E(g_abs)$weight <- abs(E(g1)$weight)
# eigen_centrality_result <- eigen_centrality(g_abs)
# eigen_centrality_result
# eigen_centrality_scores <- eigen_centrality_result$vector
# summary(eigen_centrality_result$vector)
# 
# threshold <- quantile(eigen_centrality_scores, 0.90)
# top_25_percent_nodes <- which(eigen_centrality_scores > threshold)
# subgraph <- induced_subgraph(g_abs, top_25_percent_nodes)
# print(subgraph)
# 
# plot(subgraph,
#      layout = layout,  # You can customize the layout as desired
#      edge.arrow.size = 1,
#      vertex.color = "lightskyblue1", 
#      vertex.size = 10,  # You can customize the node size, e.g., based on degree
#      vertex.frame.color = "lightskyblue1",
#      vertex.label.color = "black",
#      vertex.label.cex = 0.5,
#      vertex.label.dist = 0,
#      edge.curved = 0.2,
#      edge.color = "gray",
#      edge.width = 0.1,
#      margin = c(0, 0, 0, 0))
# V(subgraph)
