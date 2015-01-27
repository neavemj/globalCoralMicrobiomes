
## Phyloseq and r analysis of the Minimum Entropy Decomposition results ##
# Pocilloporidae samples extracted from all MED run
# 27.1.2015

# load required libraries

library("phyloseq")
library("ggplot2")
library("plyr")
library("vegan")
library('ape')
setwd("./data")

# import normal percent matrix

allShared = read.table("all.matrixPercent.txt", header=T)
rownames(allShared) = allShared[,1]
allShared = allShared[,2:length(allShared)]

# import percent matrix modified for phylogenetic tree

allSharedTree = read.table("all.matrixPercent.tree.txt", header=T)
rownames(allSharedTree) = allSharedTree[,1]
allSharedTree = allSharedTree[,2:length(allSharedTree)]

# import count matrix for alpha diversity measures

allSharedDiv = read.table("all.MED.matrixCount", header=T)
rownames(allSharedDiv) = allSharedDiv[,1]
allSharedDiv = allSharedDiv[,2:length(allSharedDiv)]

# Import normal taxonomy file from mothur

allTax = read.table('all.nodeReps.taxonomy', header=T, sep='\t')
rownames(allTax) = allTax[,1]
allTax = allTax[,3:9]
allTax = as.matrix(allTax)

# Import taxonomy file modified for the tree

allTaxTree = read.table('all.nodeReps.tree.taxonomy', header=T, sep='\t')
rownames(allTaxTree) = allTaxTree[,1]
allTaxTree = allTaxTree[,3:9]
allTaxTree = as.matrix(allTaxTree)

# import meta data

metaFile = read.table('metaData2.MED', header=T, sep='\t')
rownames(metaFile) = metaFile[,1]
metaFile = metaFile[,2:7]

# import phylogenetic tree of Endozoicomonas types

endoTree = read.tree(file='MEDNJ3.tree')

### Create phyloseq object

OTU = otu_table(allShared, taxa_are_rows = FALSE)
OTUdiv = otu_table(allSharedDiv, taxa_are_rows = FALSE)
OTUtree = otu_table(allSharedTree, taxa_are_rows = FALSE)
TAX = tax_table(allTax) 
Taxtree = tax_table(allTaxTree)
META = sample_data(metaFile)
TREE = phy_tree(endoTree)
allPhylo = phyloseq(OTU, TAX, META)
allPhyloDiv = phyloseq(OTUdiv, TAX, META)
endoTree = phyloseq(OTUtree, META, TREE)

# take a look at the pocilloporid endozoics

pVerr <- subset_samples(allPhylo, species=='Pocillopora verrucosa')
pDami <- subset_samples(allPhylo, species=='Pocillopora damicornis')
poc <- merge_phyloseq(pVerr, pDami)

sample_data(poc)$names <- factor(sample_names(poc), levels=unique(sample_names(poc)))
poc = filter_taxa(poc, function(x) mean(x) > 0.1, TRUE)

theme_set(theme_bw())
plot_bar(poc, fill="Phylum", x="names") +
  scale_y_continuous(expand = c(0,0), limits = c(0,100)) +
  facet_grid(~site, scales='free', space='free_x')





