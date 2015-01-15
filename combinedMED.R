## Phyloseq and r analysis of the Minimum Entropy Decomposition results ##
# 15.1.2015

# load required libraries

library("phyloseq")
library("ggplot2")
library("plyr")
library("vegan")
setwd("./data")

# import data

allShared = read.table("all.matrixPercent.txt", header=T)
rownames(allShared) = allShared[,1]
allShared = allShared[,2:length(allShared)]

# Import taxonomy file from mothur

allTax = read.table('all.MED.nodeReps.nr_v119.knn.taxonomy', header=T, sep='\t')
rownames(allTax) = allTax[,1]
allTax = allTax[,3:9]
allTax = as.matrix(allTax)

# import meta data

metaFile = read.table('metaData2.MED', header=T, sep='\t')
rownames(metaFile) = metaFile[,1]
metaFile = metaFile[,2:7]

### Create phyloseq object

OTU = otu_table(allShared, taxa_are_rows = FALSE)
TAX = tax_table(allTax)
META = sample_data(metaFile)
allPhylo= phyloseq(OTU, TAX, META)

# some transformations

microSub <- prune_taxa(taxa_sums(allPhylo) > 0, allPhylo)
microSubRel = transform_sample_counts(microSub, function(x) x / sum(x) )
microSubRelFilt = filter_taxa(microSubRel, function(x) mean(x) > 1e-5, TRUE)

plot_bar(allPhylo, x='species', fill="Phylum")







