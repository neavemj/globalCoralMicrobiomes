## Phyloseq and r analysis of the Minimum Entropy Decomposition results ##
# 15.1.2015

# load required libraries

library("phyloseq")
library("ggplot2")
library("plyr")
library("vegan")
library('ape')
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

# import phylogenetic tree of Endozoicomonas types

endoTree = read.tree(file='MEDNJ3.tree')

### Create phyloseq object

OTU = otu_table(allShared, taxa_are_rows = FALSE)
TAX = tax_table(allTax)
META = sample_data(metaFile)
TREE = phy_tree(endoTree)
allPhylo = phyloseq(OTU, TAX, META)

endoTree = phyloseq(OTU, META, TREE)

# some transformations

allPhylo <- prune_taxa(taxa_sums(allPhylo) > 0, allPhylo)
allPhyloFilt = filter_taxa(allPhylo, function(x) mean(x) > 0.1, TRUE)

# bar chart with everything included

plot_bar(allPhyloFilt, fill="Phylum") +
  facet_wrap(~species, scales='free')

# let's check what's happening with different Endozoicomonas OTUs

allPhyloEndo = subset_taxa(allPhylo, Genus=='Endozoicomonas')

allBarEndo <- plot_bar(allPhyloEndo, fill="Genus")
allBarEndo + facet_wrap(~site, scales='free')

# add coloring for different Endozoicomonas OTUs

tax_table(allPhyloEndo) <- cbind(tax_table(allPhyloEndo), Strain=taxa_names(allPhyloEndo))
myranks = c("Genus", "Strain")
mylabels = apply(tax_table(allPhyloEndo)[, myranks], 1, paste, sep="", collapse="_")

tax_table(allPhyloEndo) <- cbind(tax_table(allPhyloEndo), catglab=mylabels)

allPhyloEndoFilt = filter_taxa(allPhyloEndo, function(x) mean(x) > 0.1, TRUE)

plot_bar(allPhyloEndoFilt, fill="catglab") +
  facet_wrap(~species+site, scales='free')

# take a look at the pocilloporid endozoics

allPhyloEndoFiltPoc = subset_samples(allPhyloEndoFilt, species=='Pocillopora verrucosa')

plot_bar(allPhyloEndoFiltPoc, fill="catglab") +
  facet_wrap(~species+site, scales='free')

# plot a tree of all endozoicomonas MED nodes

plot_tree(endoTree, nodelabf = nodeplotboot(), ladderize='left', color='site', size='abundance', label.tips = 'taxa_names')


# mess around with merging etc.

endoTreeSite = merge_samples(endoTree, "site", fun=)
endoTreeSpecies = merge_samples(endoTree, 'species', fun=mean)

# for some reason this just sums the taxa abundances - doesn't calculate the mean! Ahh...

sample_data(endoTreeSite)$site <- factor(sample_names(endoTreeSite))
sample_data(endoTreeSpecies)$species <- factor(sample_names(endoTreeSpecies))

get_taxa(endoTreeSite)
get_taxa(endoTree)

endoTreeSiteRel = transform_sample_counts(endoTreeSite, function(x) x / sum(x) )
plot_tree(endoTree, ladderize='left', color='site', size='abundance', label.tips='taxa_names')

#endoTreeFilt= filter_taxa(endoTree, function(x) mean(x) > 0.1, TRUE)
