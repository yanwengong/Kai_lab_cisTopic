## use the lab mouse data on cisTopic
## focus on luminalSec 1 and 2
## 20191112

library(Matrix)
library(cisTopic)
library(tidyverse)
data_path <- "/Users/yanwengong/Documents/kai_lab/cisTopic/data/lab_mouse"
plot_path <- "/Users/yanwengong/Documents/kai_lab/cisTopic/plot/lab_mouse/luminalSec12"
## read in peak vs cell
## NOTE, this only contain cells from luminal based on FACS!
#peak.matrix.trimmed_luminal <- read.csv(paste(data_path, "luminal_peak_file.csv", sep = "/"))

## read in the complete, intersect meta data 
complete_meta_lu12 <- read.csv(paste(data_path, "mouse_epi_complete_int_meta.csv", sep = "/")) %>%
  mutate(barcode = as.character(barcode))

## luminal 12 from luminal FACS 669+2219+1423 = 4311 
complete_meta_lu12 %>% group_by(orig.ident) %>% count() 

## filter peak.matrix.trimmed_luminal based
barcodes_l12 <- complete_meta_lu12 %>% pull(barcode)
idx <- match(barcodes_l12, names(peak.matrix.trimmed_luminal))
length(which(!is.na(idx)))

peak.matrix.trimmed_luminal12 <- peak.matrix.trimmed_luminal[,c(na.omit(idx))] # 127752 peak x 4311 cells

## RE_START from HERE!!! ## CHANGE THE TOPIC NUMBER LARGER 
## create cisTopic object
mouse_luminal12_cisTopicObject <- createcisTopicObject(peak.matrix.trimmed_luminal12, project.name='mouse_luminal12')

## add detailed meta data and clustering information 
complete_meta_lu12_filtered <- complete_meta_lu12[c(which(!is.na(idx))),] %>%
  remove_rownames %>% column_to_rownames(var="barcode")

mouse_luminal12_cisTopicObject <- addCellMetadata(mouse_luminal12_cisTopicObject, complete_meta_lu12_filtered)

### build model
mouse_luminal12_cisTopicObject <- runModels(mouse_luminal12_cisTopicObject, 
                                            topic=c(5, 10, 15, 30), seed=987, 
                                            nCores=3, burnin = 120, iterations = 150, 
                                            addModels=FALSE) ## selected 15 topics

## note it does not saturate at topic 15, may try higher number next time 
mouse_luminal12_cisTopicObject <- selectModel(mouse_luminal12_cisTopicObject)
logLikelihoodByIter(mouse_luminal12_cisTopicObject, select=c(5, 10, 15, 30))

## select 15
mouse_luminal12_cisTopicObject <- selectModel(mouse_luminal12_cisTopicObject, select=15)



## identify cell state 

mouse_luminal12_cisTopicObject <- runtSNE(mouse_luminal12_cisTopicObject, target='cell', 
                                          seed=123, pca=F, method='Probability')

## unsupervise density cluster~~~ need to check how to control resolution 
cellassign <- modelMatSelection(mouse_luminal12_cisTopicObject, 'cell', 'Probability')
set.seed(123)
library(Rtsne)
DR <- Rtsne(t(cellassign), pca=F)
DRdist <- dist(DR$Y)
library(densityClust)
dclust <- densityClust(DRdist,gaussian=T) #Distance cutoff calculated to 4.052328 
## I changed rho from 50 to 90
dclust <- findClusters(dclust, rho = 95, delta = 2.5) #parameter can be altered
# Check thresholds
options(repr.plot.width=6, repr.plot.height=6)
plot(dclust$rho,dclust$delta,pch=20,cex=0.6,xlab='rho', ylab='delta')
points(dclust$rho[dclust$peaks],dclust$delta[dclust$peaks],col="red",pch=20,cex=0.8)
text(dclust$rho[dclust$peaks]-2,dclust$delta[dclust$peaks]+1.5,labels=dclust$clusters[dclust$peaks])
abline(v=50)
abline(h=2.5)


# Add cluster information
densityClust <- dclust$clusters
densityClust <- as.data.frame(densityClust)
rownames(densityClust) <- mouse_luminal12_cisTopicObject@cell.names
colnames(densityClust) <- 'densityClust'
densityClust[,1] <- as.factor(densityClust[,1])
mouse_luminal12_cisTopicObject <- addCellMetadata(mouse_luminal12_cisTopicObject, densityClust)


## tsne plot
mouse_luminal12_cisTopicObject@cell.data$Cluster.ID <- as.factor(mouse_luminal12_cisTopicObject@cell.data$Cluster.ID)
mouse_luminal12_cisTopicObject@cell.data$seurat_clusters <- as.factor(mouse_luminal12_cisTopicObject@cell.data$seurat_clusters)


pdf(paste(plot_path, "tsne_on_cell_nCount_nAcc_clusterID_seuratCluster.pdf", sep = "/"), width = 5, height = 4)


plotFeatures(mouse_luminal12_cisTopicObject, method='tSNE', target='cell', 
             topic_contr=NULL, colorBy=c('nCounts', 'nAcc','densityClust', 
                                         'Cluster.ID', 'seurat_clusters'), cex.legend = 0.8, 
             factor.max=.75, dim=2, legend=TRUE, cex.dot = 0.5,col.low='darkgreen', 
             col.mid='yellow', col.high='brown1', intervals=10)
dev.off()

## heatmap
pdf(paste(plot_path, "heatmap.pdf", sep = "/"), width = 5, height = 4)
cellTopicHeatmap(mouse_luminal12_cisTopicObject, method='Probability', 
                 colorBy=c('seurat_clusters'))
dev.off()

pdf(paste(plot_path, "tsne_on_cell_colorTopic.pdf", sep = "/"), width = 5, height = 4)
plotFeatures(mouse_luminal12_cisTopicObject, method='tSNE', target='cell', 
             topic_contr='Probability', colorBy=NULL, cex.legend = 0.8, 
             factor.max=.75, dim=2, legend=TRUE)
dev.off()


############## Analyze topics #####################
############## main focus on topic 5 and topic 15 ###########
############## addinng topic 3 and 12 ###########

mouse_luminal12_cisTopicObject <- getRegionsScores(mouse_luminal12_cisTopicObject, method='NormTop', scale=TRUE)

## explore the distribution of the probability of each region to topic
d <- density(mouse_luminal12_cisTopicObject@region.data$Scores_Topic5)
plot(d)

## convert the probability to binary - whether the region is enriched in the topics; by GammaFit
### thrP: Probability threshold to use as cutoff on the probability distribution when using GammaFit as method.
par(mfrow=c(2,5))
mouse_luminal12_cisTopicObject <- binarizecisTopics(mouse_luminal12_cisTopicObject, thrP=0.975, plot=TRUE)
### after binarization, there are around 3000 regions are assign to each topics

## output the bedfiles of region to topic
getBedFiles(mouse_luminal12_cisTopicObject, path=paste(data_path, 'cisTopics_asBed', sep = "/"))

## Topic visualization
### the results are saved in the slot @dr$region
mouse_luminal12_cisTopicObject <- runtSNE(mouse_luminal12_cisTopicObject, target='region', perplexity=200, 
                          check_duplicates=FALSE)

par(mfrow=c(1,1))
plotFeatures(mouse_luminal12_cisTopicObject, method='tSNE', target='region', 
             topic_contr=NULL, colorBy=c('nCells'), cex.legend = 0.8, factor.max=.75, 
             dim=2, legend=TRUE, col.low='darkgreen', 
             col.mid='yellow', col.high='brown1', intervals=10)

par(mfrow=c(2,5))
pdf(paste(plot_path, "tsne_region_colored_byTopics.pdf", sep = "/"), width = 5, height = 4)
plotFeatures(mouse_luminal12_cisTopicObject, method='tSNE', target='region', topic_contr='Z-score', 
             colorBy=NULL, cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE, cex.dot = 0.5,
             col.low='darkgreen', col.mid='yellow', col.high='brown1')
dev.off()

## Link the GO term with topics
## Annotate genes to GO terms 
#BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
#BiocManager::install("org.Mm.eg.db")
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
mouse_luminal12_cisTopicObject <- annotateRegions(mouse_luminal12_cisTopicObject, 
                                           txdb=TxDb.Mmusculus.UCSC.mm10.knownGene, 
                                           annoDb='org.Mm.eg.db')

pdf(paste(plot_path, "heatmap_region_colorTopics.pdf", sep = "/"), width = 5, height = 4)
signaturesHeatmap(mouse_luminal12_cisTopicObject, selected.signatures = 'annotation')
dev.off()

## identifying enriched GO terms per topic
mouse_luminal12_cisTopicObject <- GREAT(mouse_luminal12_cisTopicObject, genome='mm10', fold_enrichment=2, geneHits=1, 
                                 sign=0.05, request_interval=10)

pdf(paste(plot_path, "GO_topic5_15_13.pdf", sep = "/"), width = 8, height = 8)
ontologyDotPlot(mouse_luminal12_cisTopicObject, top=5, topics=c(5,15,13), 
                var.y='name', order.by='Binom_Adjp_BH')
dev.off()

## ClosestFeature to pull out the nearby gene of the open regions
regions_vs_topic_bin <- mouse_luminal12_cisTopicObject@binarized.cisTopics


## install the packages Signac and EnsDb.Mmusculus.v79
install.packages("devtools")
devtools::install_github("timoast/signac")
## err0r: Error: Failed to install 'Signac' from GitHub:
## (converted from warning) packages ‘AnnotationFilter’, ‘TFBSTools’, ‘ggbio’, 
## ‘motifmatchr’, ‘BSgenome’ are not available (for R version 3.5.1)
BiocManager::install("AnnotationFilter")
BiocManager::install("TFBSTools")
BiocManager::install("ggbio")
BiocManager::install("motifmatchr")
BiocManager::install("BSgenome")
### afterwards the package can be installed
library(Signac)



## add mouse data
BiocManager::install(c('BSgenome.Mmusculus.UCSC.mm10', 'EnsDb.Mmusculus.v79'))
BiocManager::install('EnsDb.Mmusculus.v79')
library(EnsDb.Mmusculus.v79)
# annotate the nearest genes
topic5_close_genes <- ClosestFeature(regions = (rownames(regions_vs_topic_bin$Topic5)), 
               annotation = EnsDb.Mmusculus.v79, sep = c(':', '-'))

topic15_close_genes <- ClosestFeature(regions = (rownames(regions_vs_topic_bin$Topic15)), 
                                     annotation = EnsDb.Mmusculus.v79, sep = c(':', '-'))

topic3_close_genes <- ClosestFeature(regions = (rownames(regions_vs_topic_bin$Topic3)), 
                                     annotation = EnsDb.Mmusculus.v79, sep = c(':', '-'))

topic12_close_genes <- ClosestFeature(regions = (rownames(regions_vs_topic_bin$Topic12)), 
                                     annotation = EnsDb.Mmusculus.v79, sep = c(':', '-'))

## output the dataset
output_data_path <- "/Users/yanwengong/Documents/kai_lab/cisTopic/data/lab_mouse/luminal_sec12"
write.csv(topic5_close_genes, paste(output_data_path, "topic5_close_genes.csv", sep = "/"), row.names = TRUE,
          col.names = TRUE)

write.csv(topic15_close_genes, paste(output_data_path, "topic15_close_genes.csv", sep = "/"), row.names = TRUE,
          col.names = TRUE)

write.csv(topic3_close_genes, paste(output_data_path, "topic3_close_genes.csv", sep = "/"), row.names = TRUE,
          col.names = TRUE)

write.csv(topic12_close_genes, paste(output_data_path, "topic12_close_genes.csv", sep = "/"), row.names = TRUE,
          col.names = TRUE)

## install cicero

#devtools::install_github("cole-trapnell-lab/cicero-release")
library(cicero)

## read in the enhancer region
hglft_mm10_enhancer_regions <- read.table(paste(data_path, "hglft_mm10_enhancer_regions.bed", sep = "/"),header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
## change form to A list of coordinates to be searched for overlap in the form chr_100_2000.
res = mutate(states.df,
             concated_column = paste(name, region, division, sep = '_'))
hglft_mm10_enhancer_regions_list <- hglft_mm10_enhancer_regions %>% 
  mutate(concated_col = paste(V1, V2, V3, sep = "_")) %>% 
  pull(concated_col)


## organize the form of topic regions
Topic3_regions <- regions_vs_topic_bin$Topic3 %>% rownames() 
Topic3_regions <- gsub(":", "_", gsub("-", "_", Topic3_regions))

Topic5_regions <- regions_vs_topic_bin$Topic5 %>% rownames() 
Topic5_regions <- gsub(":", "_", gsub("-", "_", Topic5_regions))

Topic12_regions <- regions_vs_topic_bin$Topic12 %>% rownames() 
Topic12_regions <- gsub(":", "_", gsub("-", "_", Topic12_regions))

Topic15_regions <- regions_vs_topic_bin$Topic15 %>% rownames() 
Topic15_regions <- gsub(":", "_", gsub("-", "_", Topic15_regions))

## I set the maxgap to -1, so that the overlap is more accurate
topic3_enhancer <- find_overlapping_coordinates(Topic3_regions, hglft_mm10_enhancer_regions_list,  maxgap = -1)
#topic3_enhancer_1 <- find_overlapping_coordinates(Topic3_regions, hglft_mm10_enhancer_regions_list,  maxgap = 0)

topic5_enhancer <- find_overlapping_coordinates(Topic5_regions, hglft_mm10_enhancer_regions_list,  maxgap = -1)
topic12_enhancer <- find_overlapping_coordinates(Topic12_regions, hglft_mm10_enhancer_regions_list,  maxgap = -1)
topic15_enhancer <- find_overlapping_coordinates(Topic15_regions, hglft_mm10_enhancer_regions_list,  maxgap = -1)

enhancer_inTopics3_df <- data.frame("enhancer" = topic3_enhancer, "topics" = "topics3") %>%
  filter(!is.na(enhancer)) %>% distinct(enhancer, .keep_all = TRUE)
enhancer_inTopics5_df <- data.frame("enhancer" = topic5_enhancer, "topics" = "topics5") %>%
  filter(!is.na(enhancer)) %>% distinct(enhancer, .keep_all = TRUE)
enhancer_inTopics12_df <- data.frame("enhancer" = topic12_enhancer, "topics" = "topics12") %>%
  filter(!is.na(enhancer)) %>% distinct(enhancer, .keep_all = TRUE)
enhancer_inTopics15_df <- data.frame("enhancer" = topic15_enhancer, "topics" = "topics15") %>% 
  filter(!is.na(enhancer)) %>% distinct(enhancer, .keep_all = TRUE)

enhancer_inTopics_df <- rbind(enhancer_inTopics3_df, rbind(enhancer_inTopics5_df, rbind(enhancer_inTopics12_df, enhancer_inTopics15_df)))
write.csv(enhancer_inTopics_df, paste(output_data_path, "region_overlap_enhancer_inTopics_df.csv", sep = "/"), col.names = TRUE)


### TEST
Topic3_regions_sorted <-Topic3_regions[order(Topic3_regions)]




## Gene accessibility scores for marker genes

region2gene <- mouse_luminal12_cisTopicObject@region.data[,'SYMBOL', drop=FALSE]
region2gene <- split(region2gene, region2gene[,'SYMBOL']) 
region2gene <- lapply(region2gene, rownames) 

selectedGenes <-c('Rspo1', 'Lalba', 'Folr1')
region2gene_subset <- region2gene[which(names(region2gene) %in% selectedGenes)]
# Compute AUC rankings based on the predictive distribution
pred.matrix <- predictiveDistribution(mouse_luminal12_cisTopicObject)
predMatSumByGene <- sapply(region2gene_subset, 
                           function(x) apply(pred.matrix[x,, drop=F], 2, sum))

rownames(predMatSumByGene) <- mouse_luminal12_cisTopicObject@cell.names
# Add to cell data
mouse_luminal12_cisTopicObject <- addCellMetadata(mouse_luminal12_cisTopicObject, predMatSumByGene)

## plot
par(mfrow=c(1,2))
pdf(paste(plot_path, "tsne_on_cell_colorGeneAccessibility.pdf", sep = "/"), width = 5, height = 4)
plotFeatures(mouse_luminal12_cisTopicObject, method='tSNE', target='cell', 
             topic_contr=NULL, colorBy=selectedGenes, 
             cex.legend = 0.8, cex.dot = 0.5, factor.max=.75, 
             dim=2, legend=TRUE, intervals=10)
dev.off()
## save data
#saveRDS(mouse_luminal12_cisTopicObject, file=paste(data_path, "mouse_luminal12_cisTopicObject.Rds", sep = "/"))
#mouse_luminal12_cisTopicObject <- readRDS(file=paste(data_path, "mouse_luminal12_cisTopicObject.Rds", sep = "/"))


## export the tsne loc 
server_ouput_path <- "/Volumes/shared/yanwen/mouse_epi_ATACseq/cisTopic/output_data"
tsne_df <- mouse_luminal12_cisTopicObject@dr$cell$tSNE %>% as.data.frame()
write.csv(tsne_df,paste(server_ouput_path, "lu12_tsne_loc.csv", sep = "/"))

## export the topic probability for each cell
cellassign <- modelMatSelection(mouse_luminal12_cisTopicObject, 'cell', 'Probability') %>% as.data.frame()
write.csv(cellassign,paste(server_ouput_path, "lu12_cellassign_probability.csv", sep = "/"))
