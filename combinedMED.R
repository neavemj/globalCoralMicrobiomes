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

cols <- c("AmericanSamoa" = "#D95F02", "Indonesia" = "#A6761D", "MaggieIs" = "#666666", "Maldives" = "#E6AB02", "Micronesia" = "#66A61E", "Ningaloo" = "#7570B3", "RedSea" = "#E7298A", "other" = "black")

# import normal percent matrix

allShared = read.table("all.7974.matrixPercent.txt", header=T)
rownames(allShared) = allShared[,1]
allShared = allShared[,2:length(allShared)]

# non-subsampled mothur shared file for alpha diversity

all3OTUshared = read.table("all.7974.0.01.shared", header=T)
rownames(all3OTUshared) = all3OTUshared[,2]
all3OTUshared = all3OTUshared[,4:length(all3OTUshared)]

# import percent matrix modified for phylogenetic tree

allSharedTree = read.table("all.7974.matrixPercent.tree.txt", header=T)
rownames(allSharedTree) = allSharedTree[,1]
allSharedTree = allSharedTree[,2:length(allSharedTree)]

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

endoTreeFile = read.tree(file='MEDNJ5.tree')

### Create phyloseq object

OTU = otu_table(allShared, taxa_are_rows = FALSE)
OTUalpha = otu_table(all3OTUshared, taxa_are_rows = FALSE)
OTUtree = otu_table(allSharedTree, taxa_are_rows = FALSE)
TAX = tax_table(allTax) 
META = sample_data(metaFile)
TREE = phy_tree(endoTreeFile)
allPhylo = phyloseq(OTU, TAX, META)
endoTree = phyloseq(OTUtree, META, TREE)
allAlpha = phyloseq(OTUalpha, META)

# Alpha diversity measures first

theme_set(theme_bw())
allDivPlot <- plot_richness(allAlpha, x = 'species', measures = c('Chao1', 'Shannon', 'observed'), color = 'site')
allDivPlot

# get rid of a few of those low corals

allAlphaTmp <- subset_samples(allAlpha, species=="seawater")
allAlphaTmp2 <- subset_samples(allAlpha, species=="Stylophora pistillata")
allAlphaTmp3 <- subset_samples(allAlpha, species=="Pocillopora verrucosa")

allAlpha2 <- merge_phyloseq(allAlphaTmp, allAlphaTmp2, allAlphaTmp3)

#serio <- subset_samples(allPhylo, species=="Seriatopora sp.")
#serio <- filter_taxa(serio, function(x) mean(x) > 0.1, TRUE)
#plot_bar(serio, fill='Genus')

# The predicted number of OTUs and diversity indecies look similar across the different coral species. I'll add some whisker plots to make this easier to see. 

allAlphaPlot2 <- plot_richness(allAlpha2, x = 'species', measures = c('Chao1', 'Shannon', 'observed'), color = 'site')
allAlphaPlot2

allAlphaPlot2 + geom_boxplot(data = allAlphaPlot2$data, aes(x = species, y = value, color = NULL), alpha = 0.1) +
  scale_color_manual(values=cols)

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

# now tree of nodes from the main two corals (plus seawater and the other seqs)

endoTreeSpist <- subset_samples(endoTree, species=='Stylophora pistillata')
endoTreePverr <- subset_samples(endoTree, species=='Pocillopora verrucosa')
endoTreeSea <- subset_samples(endoTree, species=='seawater')
endoTreeOther <- subset_samples(endoTree, species=='other')

endoTreeCorals <- merge_phyloseq(endoTreeSpist, endoTreePverr, endoTreeSea, endoTreeOther)

plot_tree(endoTreeCorals, label.tips = 'taxa_names', color='site', shape='species', size='abundance', nodelabf = nodeplotboot(100, 50, 3), ladderize='left', base.spacing = 0.01) +
  scale_color_manual(values=cols) +
  scale_shape_manual(values = c("other" = 1, "Pocillopora verrucosa" = 17, "Stylophora pistillata" = 16, "seawater" = 15)) 


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

plot_heatmap(spistPverrEndoFiltPrune, "NMDS", "bray", "site", sample.order=rownames(metaFile3))

plot_heatmap(spistPverrEndoFilt, "RDA", "none", "site", sample.order='species')

# SAVE AS 1500 X 700
# take a quick look at the Archaea samples

archaeaPhylo <- subset_taxa(allPhylo, Domain=="Archaea")
archaeaPhyloPrune = prune_samples(sample_sums(archaeaPhylo) > 0, archaeaPhylo)
plot_heatmap(archaeaPhyloPrune, "NMDS", "bray", 'species')


