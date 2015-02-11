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

allShared = read.table("all.7801.matrixPercent.txt", header=T)
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

allTax = read.table('all.7801.nodeReps.nr_v119.knn.taxonomy', header=T, sep='\t')
rownames(allTax) = allTax[,1]
allTax = allTax[,3:9]
allTax = as.matrix(allTax)

# import meta data

metaFile = read.table('metaData2.MED', header=T, sep='\t')
rownames(metaFile) = metaFile[,1]
metaFile = metaFile[,2:7]

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

plot_tree(endoTree, nodelabf = nodeplotboot(), ladderize='left', color='species', size='abundance', label.tips = 'taxa_names', base.spacing = 0.005)

# SAVE AS EPS 1500 x 1200

theme_set(theme_bw())
plot_bar(endoTree) +
  scale_y_continuous(expand = c(0,0), limits = c(0,100)) +
  facet_grid(~site, scales='free', space='free_x')


# mess around with merging etc.

#endoTreeSite = merge_samples(endoTree, "site", fun=)
#endoTreeSpecies = merge_samples(endoTree, 'species', fun=mean)

# for some reason this just sums the taxa abundances - doesn't calculate the mean! Ahh...

#sample_data(endoTreeSite)$site <- factor(sample_names(endoTreeSite))
#sample_data(endoTreeSpecies)$species <- factor(sample_names(endoTreeSpecies))
#get_taxa(endoTreeSite)
#get_taxa(endoTree)
#endoTreeSiteRel = transform_sample_counts(endoTreeSite, function(x) x / sum(x) )
#endoTreeFilt= filter_taxa(endoTree, function(x) mean(x) > 0.1, TRUE)
