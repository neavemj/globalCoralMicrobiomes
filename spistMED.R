
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

spist3OTUtax = read.table('spist.subsample.unique.good.filter.precluster.an.0.03.cons.taxonomy', header=T, sep='\t')
rownames(spist3OTUtax) = spist3OTUtax[,1]
spist3OTUtax = spist3OTUtax[,3:9]
spist3OTUtax = as.matrix(spist3OTUtax)

spist1OTUtax = read.table('spist.subsample.unique.good.filter.precluster.an.0.01.cons.taxonomy', header=T, sep='\t')
rownames(spist1OTUtax) = spist1OTUtax[,1]
spist1OTUtax = spist1OTUtax[,3:9]
spist1OTUtax = as.matrix(spist1OTUtax)

# import meta data

metaFile = read.table('metaData2.MED', header=T, sep='\t')
rownames(metaFile) = metaFile[,1]
metaFile = metaFile[,2:7]

### Create phyloseq object

OTU = otu_table(spistShared, taxa_are_rows = FALSE)
OTUs3 = otu_table(spistOTUshared, taxa_are_rows = FALSE)
OTUs1 = otu_table(spist1OTUshared, taxa_are_rows = FALSE)

TAX = tax_table(spistTax)
TAX3 = tax_table(spist3OTUtax)
TAX1 = tax_table(spist1OTUtax)

META = sample_data(metaFile)

spistPhylo = phyloseq(OTU, TAX, META)
spistOTUphylo = phyloseq(OTUs3, TAX3, META)
spist1OTUphylo = phyloseq(OTUs1, TAX1, META)

# ordination comparing MED nodes and 3% OTUs

theme_set(theme_bw())
spistMEDord <- ordinate(spistPhylo, "NMDS", "bray")
plot_ordination(spistPhylo, spistMEDord, type = 'samples', color='site', title='MED nodes')

spistOTUord <- ordinate(spistOTUphylo, "NMDS", "bray")
plot_ordination(spistOTUphylo, spistOTUord, type = 'samples', color='site', title='3% OTUs')

spist1OTUord <- ordinate(spist1OTUphylo, "NMDS", "bray")
plot_ordination(spist1OTUphylo, spist1OTUord, type = 'samples', color='site', title='1% OTUs')

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

# do some bar plots to check if MED nodes resolve endozoic types better than OTUs

spistPhyloEndo = subset_taxa(spistPhylo, Genus=='Endozoicomonas')
spistOTUphyloEndo = subset_taxa(spistOTUphylo, Genus=='Endozoicomonas(100)')
spist1OTUphyloEndo = subset_taxa(spist1OTUphylo, Genus=='Endozoicomonas(100)')

spistPhyloEndoBar <- plot_bar(spistPhyloEndo, fill="Genus")
spistPhyloEndoBar + facet_wrap(~site, scales='free')

spistOTUphyloEndoBar <- plot_bar(spistOTUphyloEndo, fill="Genus")
spistOTUphyloEndoBar + facet_wrap(~site, scales='free')

spist1OTUphyloEndoBar <- plot_bar(spist1OTUphyloEndo, fill="Genus")
spist1OTUphyloEndoBar + facet_wrap(~site, scales='free')

# add coloring for different Endozoicomonas OTUs

tax_table(spistPhyloEndo) <- cbind(tax_table(spistPhyloEndo), Strain=taxa_names(spistPhyloEndo))
tax_table(spistOTUphyloEndo) <- cbind(tax_table(spistOTUphyloEndo), Strain=taxa_names(spistOTUphyloEndo))
tax_table(spist1OTUphyloEndo) <- cbind(tax_table(spist1OTUphyloEndo), Strain=taxa_names(spist1OTUphyloEndo))
myranks = c("Genus", "Strain")

mylabelsMED = apply(tax_table(spistPhyloEndo)[, myranks], 1, paste, sep="", collapse="_")
mylabels3OTU = apply(tax_table(spistOTUphyloEndo)[, myranks], 1, paste, sep="", collapse="_")
mylabels1OTU = apply(tax_table(spist1OTUphyloEndo)[, myranks], 1, paste, sep="", collapse="_")

tax_table(spistPhyloEndo) <- cbind(tax_table(spistPhyloEndo), catglab=mylabelsMED)
tax_table(spistOTUphyloEndo) <- cbind(tax_table(spistOTUphyloEndo), catglab=mylabels3OTU)
tax_table(spist1OTUphyloEndo) <- cbind(tax_table(spist1OTUphyloEndo), catglab=mylabels1OTU)

## THE OTU TABLES ARE NOT STANDARDIZED

allPhyloEndoFilt = filter_taxa(allPhyloEndo, function(x) mean(x) > 0.1, TRUE)


plot_bar(spistPhyloEndo, fill="catglab", title='MED nodes') +
  facet_wrap(~species+site, scales='free')

plot_bar(spistOTUphyloEndo, fill="catglab", title='3% OTUs') +
  facet_wrap(~species+site, scales='free')

plot_bar(spist1OTUphyloEndo, fill="catglab", title='1% OTUs') +
  facet_wrap(~species+site, scales='free')


# take a look at the pocilloporid endozoics

allPhyloEndoFiltPoc = subset_samples(allPhyloEndoFilt, species=='Pocillopora verrucosa')

plot_bar(allPhyloEndoFiltPoc, fill="catglab") +
  facet_wrap(~species+site, scales='free')


