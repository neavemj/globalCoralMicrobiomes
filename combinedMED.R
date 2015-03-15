## Phyloseq and r analysis of the Minimum Entropy Decomposition results ##
# 15.1.2015

# load required libraries

library("phyloseq")
library("ggplot2")
library("plyr")
library("vegan")
library('ape')
library('RColorBrewer')
setwd("./data")

# import normal percent matrix

allShared = read.table("all.7974.matrixPercent.txt", header=T)
rownames(allShared) = allShared[,1]
allShared = allShared[,2:length(allShared)]

# import percent matrix modified for phylogenetic tree

allSharedTree = read.table("all.matrixPercentTree.txt", header=T)
rownames(allSharedTree) = allSharedTree[,1]
allSharedTree = allSharedTree[,2:length(allSharedTree)]

# import NoSubsampling count matrix for alpha diversity measures

allSharedDiv = read.table("all.NoSubsampling.matrixCount", header=T)
rownames(allSharedDiv) = allSharedDiv[,1]
allSharedDiv = allSharedDiv[,2:length(allSharedDiv)]

# Import normal taxonomy file from mothur

allTax = read.table('all.7974.nodeReps.nr_v119.knn.taxonomy', header=T, sep='\t')
rownames(allTax) = allTax[,1]
allTax = allTax[,3:9]
allTax = as.matrix(allTax)

# import meta data (and metaData3 for heatmap ordering)

metaFile = read.table('metaData2.MED', header=T, sep='\t')
rownames(metaFile) = metaFile[,1]
metaFile = metaFile[,2:8]

metaFile3 = read.table('metaData3.txt', header=T, sep='\t')
rownames(metaFile3) = metaFile3[,1]
metaFile3 = metaFile3[,2:8]

# import phylogenetic tree of Endozoicomonas types

endoTree = read.tree(file='MEDNJ4.tree')

### Create phyloseq object

OTU = otu_table(allShared, taxa_are_rows = FALSE)
OTUdiv = otu_table(allSharedDiv, taxa_are_rows = FALSE)
OTUtree = otu_table(allSharedTree, taxa_are_rows = FALSE)
TAX = tax_table(allTax) 
META = sample_data(metaFile)
TREE = phy_tree(endoTree)
allPhylo = phyloseq(OTU, TAX, META)
allPhyloDiv = phyloseq(OTUdiv, TAX, META)
endoTree = phyloseq(OTUtree, META, TREE)

# Alpha diversity measures first

theme_set(theme_bw())
allDivPlot <- plot_richness(allPhyloDiv, x = 'species', measures = c('Chao1', 'Shannon', 'observed'), color = 'site')
allDivPlot

# get rid of a few of those low corals

allPhyloDivTmp <- subset_samples(allPhyloDiv, species=="seawater")
allPhyloDivTmp2 <- subset_samples(allPhyloDiv, species=="Stylophora pistillata")
allPhyloDivTmp3 <- subset_samples(allPhyloDiv, species=="Pocillopora verrucosa")
allPhyloDivTmp4 <- subset_samples(allPhyloDiv, species=="Pocillopora damicornis")
allPhyloDiv2 <- merge_phyloseq(allPhyloDivTmp, allPhyloDivTmp2, allPhyloDivTmp3, allPhyloDivTmp4)

# The predicted number of OTUs and diversity indecies look similar across the different coral species. I'll add some whisker plots to make this easier to see. 

allDivPlot2 <- plot_richness(allPhyloDiv2, x = 'species', measures = c('Chao1', 'Shannon', 'observed'), color = 'site')
allDivPlot2

allDivPlot2 + geom_boxplot(data = allDivPlot2$data, aes(x = species, y = value, color = NULL), alpha = 0.1)

# some transformations - top line is to keep sample order in bar graphs

sample_data(allPhylo)$names <- factor(sample_names(allPhylo), levels=unique(sample_names(allPhylo)))
allPhylo <- prune_taxa(taxa_sums(allPhylo) > 0, allPhylo)
allPhyloFilt = filter_taxa(allPhylo, function(x) mean(x) > 0.1, TRUE)

# bar chart with everything included

theme_set(theme_bw())
plot_bar(allPhyloFilt, fill="Phylum") +
  scale_y_continuous(expand = c(0,0), limits = c(0,100)) +
  facet_grid(~site, scales='free', space='free_x')

# plot a tree of all endozoicomonas MED nodes

plot_tree(endoTree, nodelabf = nodeplotboot(), ladderize='left', color='pocType', size='abundance', label.tips = 'taxa_names', base.spacing = 0.005)

# SAVE AS EPS 1500 x 1200

theme_set(theme_bw())
plot_bar(endoTree) +
  scale_y_continuous(expand = c(0,0), limits = c(0,100)) +
  facet_grid(~site, scales='free', space='free_x')

# plot a heat map to show these differences a bit better

plot_heatmap(allPhylo, "RDA", "none", "species", "none")

allPhyloEndo = subset_taxa(allPhylo, Genus=='Endozoicomonas')

coralPhyloEndo = subset_samples(allPhyloEndo, species!='seawater')
spistPhyloEndo <- subset_samples(allPhyloEndo, species=='Stylophora pistillata')
pVerrPhyloEndo <- subset_samples(allPhyloEndo, species=='Pocillopora verrucosa')
spistPverrEndo <- merge_phyloseq(spistPhyloEndo, pVerrPhyloEndo)

spistPverrEndoFilt = filter_taxa(spistPverrEndo, function(x) mean(x) > 0.0, TRUE)
spistPverrEndoFiltPrune = prune_samples(sample_sums(spistPverrEndoFilt) > 0, spistPverrEndoFilt)

spistPverrEndoFiltPrune$names <- factor(spistPverrEndoFiltPrune$Sample, levels=rownames(metaFile), ordered = TRUE)

plot_heatmap(spistPverrEndoFiltPrune, "NMDS", "bray", sample.order=rownames(metaFile3))

plot_heatmap(spistPverrEndoFilt, "RDA", "none", "site", sample.order='species')

# SAVE AS 1500 X 700
# take a quick look at the Archaea samples

archaeaPhylo <- subset_taxa(allPhylo, Domain=="Archaea")
tmp = prune_samples(sample_sums(archaeaPhylo) > 0, archaeaPhylo)
plot_heatmap(archaeaPhylo, "RDA", "none", 'species')

plot_heatmap(tmp)

