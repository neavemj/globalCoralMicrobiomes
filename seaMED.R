
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

# import count matrix for alpha diversity measures

allSharedDiv = read.table("all.MED.matrixCount", header=T)
rownames(allSharedDiv) = allSharedDiv[,1]
allSharedDiv = allSharedDiv[,2:length(allSharedDiv)]

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
OTUdiv = otu_table(allSharedDiv, taxa_are_rows = FALSE)
TAX = tax_table(allTax) 
META = sample_data(metaFile)
allPhylo = phyloseq(OTU, TAX, META)
allPhyloDiv = phyloseq(OTUdiv, TAX, META)

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

# ordination for the seawater samples

theme_set(theme_bw())
seaOrd <- ordinate(sea, "NMDS", "bray")
plot_ordination(sea, seaOrd, type = 'samples', color='site', title='seawater')

# overlay the OTUs onto this

seaFilt = filter_taxa(sea, function(x) mean(x) > 0.1, TRUE)

#Top10OTUs = names(sort(taxa_sums(seaFilt), TRUE)[1:10])
#seaFilt10 = prune_taxa(Top10OTUs, seaFilt)

seaOrd <- ordinate(seaFilt, "NMDS", "bray")
plot_ordination(seaFilt, seaOrd, type = 'split', color='site', title='seawater', label="Genus")

# calculate indicator species

seaDist <- vegdist(otu_table(seaFilt), method="bray")
seaInd <- indspc(otu_table(seaFilt), seaDist, numitr=1000)

# sort by most important indicator species

seaIndVals <- seaInd$vals
seaIndValsSorted <- seaIndVals[order(-seaIndVals$indval),,drop=FALSE]

indval numocc  pval
MED000002578 0.7439444     22 0.001
MED000006214 0.6228868     34 0.001
MED000002576 0.6095937     28 0.001
MED000004147 0.5992344     30 0.001
MED000005444 0.5828933     34 0.001
MED000005550 0.5814497     39 0.001
MED000007083 0.5718403     32 0.001
MED000006657 0.5716694     41 0.001
MED000006596 0.5716694     41 0.001
MED000005833 0.5716694     41 0.001

# overlay these onto ordination? 


seaOrdTop10 <- seaOrd$species[c("MED000002578","MED000006214","MED000002576","MED000004147","MED000005444","MED000005550","MED000007083","MED000006657","MED000006596","MED000005833"),]

arrowmatrix = seaOrdTop10
arrowdf <- data.frame(labels = rownames(arrowmatrix), arrowmatrix)

# get taxonomic information from the original tax file

arrowdf <- data.frame(labels = allTax[rownames(arrowmatrix),"Family"], arrowmatrix)

arrowmap <- aes(xend = MDS1, yend = MDS2, x = 0, y = 0, alpha=0.5, shape = NULL, color = NULL, label = labels)
labelmap <- aes(x = MDS1, y = MDS2 + 0.04, shape = NULL, color = NULL, label = labels, size=1.5)
arrowhead = arrow(length = unit(0.02, "npc"))


seaFiltPlot <- plot_ordination(seaFilt, seaOrd, type = 'samples', color='site', title='seawater')

seaFiltPlotArrow <- seaFiltPlot + geom_segment(arrowmap, size = 0.5, data = arrowdf, color = "black",  arrow = arrowhead, show_guide = FALSE) + geom_text(labelmap, size = 3, data = arrowdf)
seaFiltPlotArrow

# check if sites significantly different

df = as(sample_data(sea), "data.frame")
d = distance(sea, "bray")
seaAdonis = adonis(d ~ site + reef, df)
seaAdonis




