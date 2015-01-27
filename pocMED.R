
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

# take a look at the pocilloporid endozoics

pVerr <- subset_samples(allPhylo, species=='Pocillopora verrucosa')
pDami <- subset_samples(allPhylo, species=='Pocillopora damicornis')
poc <- merge_phyloseq(pVerr, pDami)

sample_data(poc)$names <- factor(sample_names(poc), levels=unique(sample_names(poc)))
sample_data(pVerr)$names <- factor(sample_names(pVerr), levels=unique(sample_names(pVerr)))
sample_data(pDami)$names <- factor(sample_names(pDami), levels=unique(sample_names(pDami)))

pocFilt = filter_taxa(poc, function(x) mean(x) > 0.1, TRUE)
pVerrFilt = filter_taxa(pVerr, function(x) mean(x) > 0.2, TRUE)

#colourCount = length(unique(mtcars$hp))
#getPalette = colorRampPalette(brewer.pal(9, "Set1"))

pVerrFiltGlom <- tax_glom(pVerrFilt, taxrank="Phylum")

theme_set(theme_bw())
plot_bar(pVerrFiltGlom, fill="Phylum", x="names") +
  scale_y_continuous(expand = c(0,0), limits = c(0,100)) +
  facet_grid(~site, scales='free', space='free_x')

# EXPORT AS EPS. 1500 x 600
# let's check what's happening with different Endozoicomonas OTUs (Pocillopora combined)

pocEndo = subset_taxa(poc, Genus=='Endozoicomonas')

# add coloring for different Endozoicomonas OTUs 

tax_table(pocEndo) <- cbind(tax_table(pocEndo), Strain=taxa_names(pocEndo))
myranks = c("Genus", "Strain")
mylabels = apply(tax_table(pocEndo)[, myranks], 1, paste, sep="", collapse="_")
tax_table(pocEndo) <- cbind(tax_table(pocEndo), catglab=mylabels)

sample_data(pocEndo)$names <- factor(sample_names(pocEndo), levels=unique(sample_names(pocEndo)))

pocEndoFilt = filter_taxa(pocEndo, function(x) mean(x) > 0.2, TRUE)

theme_set(theme_bw())
plot_bar(pocEndoFilt, fill="catglab", x="names") +
  scale_y_continuous(expand = c(0,0), limits = c(0,100)) +
  scale_fill_brewer(type='qual', palette = 'Set1') +
  facet_grid(~site, scales='free', space='free_x')


# let's check what's happening with different Endozoicomonas OTUs (Pocillopora verrucosa)

pVerrEndo = subset_taxa(pVerr, Genus=='Endozoicomonas')

# add coloring for different Endozoicomonas OTUs 

tax_table(pVerrEndo) <- cbind(tax_table(pVerrEndo), Strain=taxa_names(pVerrEndo))
myranks = c("Genus", "Strain")
mylabels = apply(tax_table(pVerrEndo)[, myranks], 1, paste, sep="", collapse="_")
tax_table(pVerrEndo) <- cbind(tax_table(pVerrEndo), catglab=mylabels)

sample_data(pVerrEndo)$names <- factor(sample_names(pVerrEndo), levels=unique(sample_names(pVerrEndo)))

pVerrEndoFilt = filter_taxa(pVerrEndo, function(x) mean(x) > 0.2, TRUE)

theme_set(theme_bw())
plot_bar(pVerrEndoFilt, fill="catglab", x="names") +
  scale_y_continuous(expand = c(0,0), limits = c(0,100)) +
  scale_fill_brewer(type='qual', palette = 'Set1') +
  facet_grid(~site, scales='free', space='free_x')

# SAVE EPS 1500 x 700



