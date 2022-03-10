###
##
# Cluster ensemble to create bundles of ecosystem services
#   i. hierarchical clustering using euclidean distance and absolute correlation
#   ii. k-means
#   iii. k-medoids (pam)
#   iv. random forest using k-means and k-medoids to cluster proximity matrix
# Required input: The file 'DES.csv'
# Last update: 10/03/2022
# Questions: martin.reader@geo.uzh.ch
##
###



# Load packages
library(dendextend)
library(cluster)      # clustering algorithms (pam)
library(factoextra)   # clustering visualization (fviz)
library(randomForest)
library(tidyverse)

# Load data
DES = read.csv("DES.csv", header = TRUE, sep = ",")
DESes = DES[,c(15:63)]  # select only ES
DESes_scaled = scale(DESes)
DESes_transposed_scaled = t(DESes_scaled)

# Create tibble and insert new ES names
DES_tibble = as_tibble (DESes_transposed_scaled)
DES_names = colnames(DESes)
row.names(DES_tibble) <- DES_names

# Set number of clusters
k = 4



###
##
# i. Hierarchical clustering

##
# Euclidean distance

es_hc_eucl = hclust(dist(DESes_transposed_scaled))
hc_eucl_clust <- cutree(es_hc_eucl, k = k)

# Plot
dend <- DESes_transposed_scaled %>%
	dist %>% hclust %>% as.dendrogram
dend %>% set("labels", DES_names) %>%
	set("labels_cex", 1.45) %>% 
	set("branches_k_color", k = k) %>% 
	plot(horiz = TRUE)

##
# Absolute correlation

es_hc_corr = hclust(as.dist(1-abs(cor(DESes_scaled))))
hc_corr_clust <- cutree(es_hc_corr, k = k)

# Plot
dend <- es_hc_corr %>% as.dendrogram
dend %>% set("labels", DES_names) %>%
	set("labels_cex", 1.45) %>% 
	set("branches_k_color", k = k) %>% 
	plot(horiz = TRUE)



###
##
# ii. k-means

k_clust <- kmeans(DESes_transposed_scaled, centers = k, nstart = 5000)
names(k_clust$cluster) <- DES_names

# Plot
fviz_cluster(k_clust, data = DES_tibble, ggtheme = theme_minimal(base_size = 42.5), 
			 main = FALSE,
			 palette = c("darkviolet", "blue3", "cyan4", "chartreuse4", "goldenrod3", "orange3", "deeppink4"),
			 pointsize = 6, shape = 20, labelsize = 35, repel = TRUE) + 
	         theme(legend.position = "none")



###
##
# iii. k-medoids (pam)

pam_clust <- pam(DES_tibble, k)
names(pam_clust$cluster) <- DES_names

# Plot
fviz_cluster(pam_clust, data = DES_tibble, ggtheme = theme_minimal(base_size = 42.5), 
			 main = FALSE,
			 palette = c("darkviolet", "blue3", "cyan4", "chartreuse4", "goldenrod3", "orange3", "deeppink4"),
			 pointsize = 6, shape = 20, labelsize = 35, repel = TRUE) + 
	         theme(legend.position = "none")



###
##
# iv. Random forest

DESes_rf <- randomForest(DESes_transposed_scaled, ntree = 5000)
prox <- DESes_rf$proximity

##
# k-means of random forest

k_rf_clust <- kmeans(prox, centers = k, nstart = 5000)

# Plot
fviz_cluster(k_rf_clust, data = DES_tibble, ggtheme = theme_minimal(base_size = 42.5), 
			 main = FALSE,
			 palette = c("darkviolet", "blue3", "cyan4", "chartreuse4", "goldenrod3", "orange3", "deeppink4"),
			 pointsize = 6, shape = 20, labelsize = 35, repel = TRUE) + 
	         theme(legend.position = "none")

##
# k-medoids (pam) of random forest

prox_tibble = as_tibble(prox)
row.names(prox_tibble) <- DES_names
colnames(prox_tibble) <- DES_names

pam_rf_clust <- pam(prox_tibble, k)

# Plot
fviz_cluster(pam_rf_clust, data = DES_tibble, ggtheme = theme_minimal(base_size = 42.5), 
			 main = FALSE,
			 palette = c("darkviolet", "blue3", "cyan4", "chartreuse4", "goldenrod3", "orange3", "deeppink4"),
			 pointsize = 6, shape = 20, labelsize = 35, repel = TRUE) + 
	         theme(legend.position = "none")



###
##
# Combine and export clustering results (use Cluster_relabel.m to relabel consistently

ESclusters = cbind(hc_corr_clust, hc_eucl_clust, k_clust$cluster, pam_clust$cluster, k_rf_clust$cluster, pam_rf_clust$cluster)
colnames(ESclusters)[3] <- "k"
colnames(ESclusters)[4] <- "pam"
colnames(ESclusters)[5] <- "k_rf"
colnames(ESclusters)[6] <- "pam_rf"
write.csv(ESclusters, "DESclusters.csv")