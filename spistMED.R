
## Phyloseq and r analysis of the Minimum Entropy Decomposition results for spist ##
# 19.1.2015

# load required libraries

library("phyloseq")
library("ggplot2")
library("plyr")
library("vegan")
setwd("./data")

# import MED data, 3% OTU and 1% OTU data for comparison

spistShared = read.table("spist.matrixPercent.txt", header=T)
rownames(spistShared) = spistShared[,1]
spistShared = spistShared[,2:length(spistShared)]

spistOTUshared = read.table("spist.subsample.unique.good.filter.precluster.an.0.03.pick.shared", header=T)
rownames(spistOTUshared) = spistOTUshared[,2]
spistOTUshared = spistOTUshared[,4:length(spistOTUshared)]

spist1OTUshared = read.table("spist.subsample.0.01.pick.shared", header=T)
rownames(spist1OTUshared) = spist1OTUshared[,2]
spist1OTUshared = spist1OTUshared[,4:length(spist1OTUshared)]

# Import taxonomy file from mothur

spistTax = read.table('spist.MED.nodeReps.nr_v119.knn.taxonomy', header=T, sep='\t')
rownames(spistTax) = spistTax[,1]
spistTax = spistTax[,3:9]
spistTax = as.matrix(spistTax)

# import meta data

metaFile = read.table('metaData2.MED', header=T, sep='\t')
rownames(metaFile) = metaFile[,1]
metaFile = metaFile[,2:7]

### Create phyloseq object

OTU = otu_table(spistShared, taxa_are_rows = FALSE)
OTUs3 = otu_table(spistOTUshared, taxa_are_rows = FALSE)
OTUs1 = otu_table(spist1OTUshared, taxa_are_rows = FALSE)

TAX = tax_table(spistTax)
META = sample_data(metaFile)

spistPhylo = phyloseq(OTU, TAX, META)
spistOTUphylo = phyloseq(OTUs3, META)
spist1OTUphylo = phyloseq(OTUs1, META)

# ordination comparing MED nodes and 3% OTUs

theme_set(theme_bw())
spistMEDord <- ordinate(spistPhylo, "NMDS", "bray")
plot_ordination(spistPhylo, spistMEDord, type = 'samples', color='site')

spistOTUord <- ordinate(spistOTUphylo, "NMDS", "bray")
plot_ordination(spistOTUphylo, spistOTUord, type = 'samples', color='site')

spist1OTUord <- ordinate(spist1OTUphylo, "NMDS", "bray")
plot_ordination(spist1OTUphylo, spist1OTUord, type = 'samples', color='site')

# stress is fairly high (~0.25). Try to log transform the OTUs

spistSharedLog = decostand(spistShared, method='log', logbase=10)
mds <- metaMDS(spistSharedLog, distance="bray")
ggplot() + 
       geom_point(aes(x=mds$points[,1], y=mds$points[,2])) +
       scale_y_reverse()

# still reasonably high (0.2)
# do some more transforming to get stress down

spistPhyloSqrt = transform_sample_counts(spistPhylo, function(x) sqrt(x))
spistPhyloSqrtord <- ordinate(spistPhyloSqrt, "NMDS", "bray")
plot_ordination(spistPhyloSqrt, spistPhyloSqrtord, type = 'samples', color='site')

spistPhyloFilt = filter_taxa(spistPhylo, function(x) mean(x) > 1, TRUE)
spistPhyloFiltord <- ordinate(spistPhyloFilt, "NMDS", "bray")
plot_ordination(spistPhyloFilt, spistPhyloFiltord, type = 'samples', color='site')

# sqrt doesn't make a difference because this is automatically done in the ordination proceedure anyway
# filtering taxa doesn't make much difference until almost all taxa are gone (only 20 or so left)



