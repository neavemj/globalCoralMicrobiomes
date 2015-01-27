

## Phyloseq and r analysis of the Minimum Entropy Decomposition results ##
# Stylophora pistillata samples extracted from all MED run
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

# Import normal taxonomy file from mothur

allTax = read.table('all.nodeReps.taxonomy', header=T, sep='\t')
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
allPhylo = phyloseq(OTU, TAX, META)

# bar chart of s.pistillata phyla or class - the names parameter preserves order

spist <- subset_samples(allPhylo, species=='Stylophora pistillata')

sample_data(spist)$names <- factor(sample_names(spist), levels=unique(sample_names(spist)))

spistFilt = filter_taxa(spist, function(x) mean(x) > 0.1, TRUE)

theme_set(theme_bw())
plot_bar(spistFilt, fill="Class", x="names") +
  scale_y_continuous(expand = c(0,0), limits = c(0,100)) +
  facet_grid(~site, scales='free', space='free_x')

# SAVE PLOT: EPS 1500 x 1174

#Top100OTUs = names(sort(taxa_sums(spist), TRUE)[1:100])
#spist100 = prune_taxa(Top50OTUs, spist)

# let's check what's happening with different Endozoicomonas OTUs

allPhyloEndo = subset_taxa(allPhylo, Genus=='Endozoicomonas')
spistEndo = subset_taxa(spist, Genus=='Endozoicomonas')

# add coloring for different Endozoicomonas OTUs

tax_table(allPhyloEndo) <- cbind(tax_table(allPhyloEndo), Strain=taxa_names(allPhyloEndo))
tax_table(spistEndo) <- cbind(tax_table(spistEndo), Strain=taxa_names(spistEndo))

myranks = c("Genus", "Strain")

mylabels = apply(tax_table(allPhyloEndo)[, myranks], 1, paste, sep="", collapse="_")
mylabelsSpist = apply(tax_table(spistEndo)[, myranks], 1, paste, sep="", collapse="_")

tax_table(allPhyloEndo) <- cbind(tax_table(allPhyloEndo), catglab=mylabels)
tax_table(spistEndo) <- cbind(tax_table(spistEndo), catglab=mylabelsSpist)

allPhyloEndoFilt = filter_taxa(allPhyloEndo, function(x) mean(x) > 0.1, TRUE)

sample_data(spistEndo)$names <- factor(sample_names(spistEndo), levels=unique(sample_names(spistEndo)))

spistEndoFilt = filter_taxa(spistEndo, function(x) mean(x) > 0.2, TRUE)

theme_set(theme_bw())
plot_bar(spistEndoFilt, fill="catglab", x="names") +
  scale_y_continuous(expand = c(0,0), limits = c(0,100)) +
  scale_fill_brewer(type='qual', palette = 'Dark2') +
  facet_grid(~site, scales='free', space='free_x')

# SAVE PLOT: EPS 1500 x 1174. Greater than 0.2%







