
## Phyloseq and r analysis of the Minimum Entropy Decomposition results ##
# seawater samples extracted from #7 all MED run
# 28.1.2015

# load required libraries

library("phyloseq")
library("ggplot2")
library("plyr")
library("vegan")
library('ape')
library("labdsv")
library("grid")
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

sea <- subset_samples(allPhylo, species=='seawater')
seaDiv <- subset_samples(allPhyloDiv, species=='seawater')

sample_data(sea)$names <- factor(sample_names(sea), levels=unique(sample_names(sea)))

seaFilt = filter_taxa(sea, function(x) mean(x) > 0.1, TRUE)

seaFiltGlom <- tax_glom(seaFilt, taxrank="Phylum")

theme_set(theme_bw())
plot_bar(seaFiltGlom, fill="Phylum", x="names") +
  scale_y_continuous(expand = c(0,0), limits = c(0,100)) +
  facet_grid(~site, scales='free', space='free_x')

# SAVE EPS 1500 x 700

# see if SIMPER works

metaFileSIMP <- subset(metaFile, site != 'NA')
seaDistSIMP <- merge(otu_table(sea), metaFileSIMP, by="row.names")

rownames(seaDistSIMP) <- seaDistSIMP[,1]
OTUsSIMP <- seaDistSIMP[,2:743]
factorsSIMP <- seaDistSIMP[,744:749]

seaDistSIMP[,750]

spistSIMP <- simper(OTUsSIMP, factorsSIMP$site)

allTax["MED000008042", "Genus"]

# SIMPROF

seaSIMPROF <- simprof(OTUsSIMP, num.expected=10, num.simulated=9, method.cluster='average', method.distance='braycurtis', method.transform='squareroot', alpha=0.05, sample.orientation='row', silent=FALSE)

simprof.plot(seaSIMPROF, leafcolors=NA, plot=TRUE, fill=TRUE, leaflab="perpendicular", siglinetype=1)





