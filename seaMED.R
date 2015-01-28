
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

# ordination for the seawater samples

theme_set(theme_bw())
seaOrd <- ordinate(sea, "NMDS", "bray")
plot_ordination(sea, seaOrd, type = 'samples', color='site', title='seawater')

# overlay the OTUs onto this

seaFilt = filter_taxa(sea, function(x) mean(x) > 0.1, TRUE)

#Top10OTUs = names(sort(taxa_sums(seaFilt), TRUE)[1:10])
#seaFilt10 = prune_taxa(Top10OTUs, seaFilt)

seaOrd <- ordinate(sea, "NMDS", "bray")
plot_ordination(sea, seaOrd, type = 'split', color='site', title='seawater', label="Genus")

# calculate indicator species

seaDist <- vegdist(otu_table(sea), method="bray")
seaInd <- indspc(otu_table(sea), seaDist, numitr=1000)

# sort by most important indicator species

seaIndVals <- seaInd$vals
seaIndValsSorted <- seaIndVals[order(-seaIndVals$indval),,drop=FALSE]

indval numocc  pval
MED000006898 0.82280849      2 0.027
MED000001740 0.79870842      2 0.046
MED000004938 0.78373661      2 0.054
MED000007874 0.75106607     14 0.001
MED000005832 0.74325044      2 0.080
MED000003752 0.73354346      9 0.001
MED000006835 0.73130608      6 0.001
MED000002578 0.70127820     22 0.001
MED000007846 0.70127820     22 0.001
MED000000600 0.69992476     23 0.001

# overlay these onto ordination

seaOrdFiltTop10 <- seaOrd$species[c("MED000002578","MED000006214","MED000002576","MED000004147","MED000005444","MED000005550","MED000007083","MED000006657","MED000006596","MED000005833"),]

seaOrdTop10 <- seaOrd$species[c("MED000006898","MED000001740","MED000004938","MED000007874","MED000005832","MED000003752","MED000006835","MED000002578","MED000007846","MED000000600"),]

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
seaAdonis = adonis(d ~ site, df, permutations=999)
seaAdonis




