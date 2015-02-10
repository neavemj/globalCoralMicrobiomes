
## Phyloseq and r analysis of the Minimum Entropy Decomposition results ##
# Stylophora pistillata samples extracted from all MED run
# 27.1.2015

# load required libraries

library("phyloseq")
library("ggplot2")
library("plyr")
library("vegan")
library('ape')
library('grid')
library('RColorBrewer')
setwd("./data")

# generate some colors to be consistent

#display.brewer.all()
#display.brewer.pal(n = 8, name = 'Dark2')
#brewer.pal(n = 8, name = "Dark2")

cols <- c("AmericanSamoa" = "#D95F02", "Indonesia" = "#A6761D", "MaggieIs" = "#666666", "Maldives" = "#E6AB02", "Micronesia" = "#66A61E", "Ningaloo" = "#7570B3", "RedSea" = "#E7298A")

# import normal percent matrix

allShared = read.table("all.7801.matrixPercent.txt", header=T)
rownames(allShared) = allShared[,1]
allShared = allShared[,2:length(allShared)]

# Import normal taxonomy file from mothur

allTax = read.table('all.7801.nodeReps.nr_v119.knn.taxonomy', header=T, sep='\t')
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
# transforming to normal matix and as.factor keeps taxa stacked consistently

spist <- subset_samples(allPhylo, species=='Stylophora pistillata')

sample_data(spist)$names <- factor(sample_names(spist), levels=unique(sample_names(spist)))

spistFilt = filter_taxa(spist, function(x) mean(x) > 0.1, TRUE)

spistFiltGlom <- tax_glom(spistFilt, taxrank="Phylum")
physeqdf <- psmelt(spistFiltGlom)

# get total abundance so can make a 'other' column
# had to add ^ and $ characters to make sure grep matches whole word

taxLevel = "Phylum"

for (j in unique(physeqdf$Sample)) {
  jFirst = paste('^', j, sep='')
  jBoth = paste(jFirst, '$', sep='')
  rowNumbers = grep(jBoth, physeqdf$Sample)
  otherValue = 100 - sum(physeqdf[rowNumbers,"Abundance"])
  newRow = (physeqdf[rowNumbers,])[1,]
  newRow[,taxLevel] = "other"
  newRow[,"Abundance"] = otherValue
  physeqdfOther <- rbind(physeqdfOther, newRow)
}

# need to create my own ggplot colors then replace the last one with gray
# this will ensure that the 'other' category is gray

gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}

ggCols <- gg_color_hue(length(unique(physeqdfOther$Phylum)))
ggCols <- head(ggCols, n=-1)

theme_set(theme_bw())
ggplot(physeqdfOther, aes(x=Sample, y=Abundance, fill=Phylum, order = as.factor(Phylum))) +
  geom_bar(stat="identity", colour="black") +
  scale_fill_manual(values=c(ggCols, "gray")) +
  scale_y_continuous(expand = c(0,0), limits = c(0,100)) +
  facet_grid(~site, scales='free', space='free_x') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


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


# SAVE PLOT: EPS 1500 x 600. Greater than 0.2%

# ordination for the spist samples

theme_set(theme_bw())
spistOrd <- ordinate(spist, "NMDS", "bray")
plot_ordination(spist, spistOrd, type = 'samples', color='site', title='spist')

# overlay the OTUs onto this

spistFilt = filter_taxa(spist, function(x) mean(x) > 0.1, TRUE)

#Top10OTUs = names(sort(taxa_sums(spistFilt), TRUE)[1:10])
#spistFilt10 = prune_taxa(Top10OTUs, spistFilt)

spistOrd <- ordinate(spist, "NMDS", "bray")
plot_ordination(spist, spistOrd, type = 'split', color='site', title='spist', label="Genus")

# calculate indicator species

spistDist <- vegdist(otu_table(spist), method="bray")
spistInd <- indspc(otu_table(spist), spistDist, numitr=1000)

scores.spist <- scores(spistSIMP, display='species')

spist.cor <- cor(otu_table(spist), method='pearson')

# sort by most important indicator species

spistIndVals <- spistInd$vals
spistIndValsSorted <- spistIndVals[order(-spistIndVals$indval),,drop=FALSE]

indval numocc  pval
MED000007667 0.5318256470      7 0.001
MED000007525 0.4988038081      4 0.001
MED000002495 0.4936241515      2 0.006
MED000002426 0.4936241515      2 0.006
MED000005726 0.4634710015      2 0.009
MED000005721 0.4160763361      5 0.001
MED000007914 0.4160763361      5 0.001
MED000005203 0.4133567153      2 0.015
MED000005495 0.4133567153      2 0.015
MED000006562 0.3998108740      5 0.001

# overlay these onto ordination

spistOrdTop10 <- spistOrd$species[c(rownames(spistIndValsSorted)[1:10]),]


spistMothurSpearman <- c("MED000006776",  "MED000009228",  "MED000005973",  "MED000009485",  "MED000005839",  "MED000007679",  "MED000009504",  "MED000009644", "MED000010328",  "MED000009696")


arrowmatrix <- data.frame(labels = spistMothurSpearman)
rownames(arrowmatrix) <- arrowmatrix[,1]
arrowdf <- data.frame(arrowmatrix)

# get taxonomic information from the original tax file

arrowdf <- data.frame(labels = allTax[rownames(arrowmatrix),"Genus"], arrowmatrix)
arrowdf <- cbind(arrowdf, spistOrd$species[spistMothurSpearman, ])

# also combine MED node with genus taxonomy for labelling

mylabels = apply(arrowdf[, c("labels", "labels.1")], 1, paste, sep="", collapse="_")
arrowdf <- cbind(arrowdf, catglab=mylabels)

# now create labels, arrows, etc.

arrowmap <- aes(xend = MDS1, yend = MDS2, x = 0, y = 0, shape = NULL, color = NULL, label = labels)
labelmap <- aes(x = MDS1, y = MDS2 + 0.04, shape = NULL, color = NULL, label = catglab, size=1.5)
arrowhead = arrow(length = unit(0.02, "npc"))

spistFiltPlot <- plot_ordination(spist, spistOrd, type = 'samples', color='site', title='spist') +
  scale_color_manual(values=cols)

spistFiltPlotArrow <- spistFiltPlot + 
  geom_segment(arrowmap, size = 1, data = arrowdf, color = "grey",  arrow = arrowhead, show_guide = FALSE) + 
  scale_color_manual(values=cols) +
  geom_text(labelmap, size = 2, data = arrowdf)

spistFiltPlotArrow








