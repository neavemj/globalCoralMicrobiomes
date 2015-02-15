
## Phyloseq and r analysis of the Minimum Entropy Decomposition results ##
# Pocilloporidae samples extracted from all MED run
# 27.1.2015

# load required libraries

library("phyloseq")
library("ggplot2")
library("plyr")
library("vegan")
library('ape')
library('grid')
library(RColorBrewer)
setwd("./data")

# generate some colors to be consistent

#display.brewer.all()
#display.brewer.pal(n = 8, name = 'Dark2')
#brewer.pal(n = 8, name = "Dark2")

cols <- c("AmericanSamoa" = "#D95F02", "Indonesia" = "#A6761D", "MaggieIs" = "#666666", "Maldives" = "#E6AB02", "Micronesia" = "#66A61E", "Ningaloo" = "#7570B3", "RedSea" = "#E7298A")

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

metaFile = read.table('metaData2.csv', header=T, sep=',')
rownames(metaFile) = metaFile[,1]
metaFile = metaFile[,2:8]

### Create phyloseq object

OTU = otu_table(allShared, taxa_are_rows = FALSE)
TAX = tax_table(allTax) 
META = sample_data(metaFile)
allPhylo = phyloseq(OTU, TAX, META)

# take a look at the pocilloporid samples

pVerr <- subset_samples(allPhylo, species=='Pocillopora verrucosa')
pDami <- subset_samples(allPhylo, species=='Pocillopora damicornis')
poc <- merge_phyloseq(pVerr, pDami)

sample_data(poc)$names <- factor(sample_names(poc), levels=unique(sample_names(poc)))
sample_data(pVerr)$names <- factor(sample_names(pVerr), levels=unique(sample_names(pVerr)))
sample_data(pDami)$names <- factor(sample_names(pDami), levels=unique(sample_names(pDami)))

pocFilt = filter_taxa(poc, function(x) mean(x) > 0.2, TRUE)
pVerrFilt = filter_taxa(pVerr, function(x) mean(x) > 0.2, TRUE)
pDamiFilt = filter_taxa(pDami, function(x) mean(x) > 1, TRUE)

# phylum / class bars
# transforming to normal matix and as.factor keeps taxa stacked consistently

taxLevel <- "Class"

pocFiltGlom <- tax_glom(pocFilt, taxrank=taxLevel)
pocdf <- psmelt(pocFiltGlom)

# get total abundance so can make a 'other' column
# had to add ^ and $ characters to make sure grep matches whole word

pocdfOther <- pocdf

for (j in unique(pocdf$Sample)) {
  jFirst = paste('^', j, sep='')
  jBoth = paste(jFirst, '$', sep='')
  rowNumbers = grep(jBoth, pocdf$Sample)
  otherValue = 100 - sum(pocdf[rowNumbers,"Abundance"])
  newRow = (pocdf[rowNumbers,])[1,]
  newRow[,taxLevel] = "other"
  newRow[,"Abundance"] = otherValue
  pocdfOther <- rbind(pocdfOther, newRow)
}

# need to create my own ggplot colors then replace the last one with gray
# this will ensure that the 'other' category is gray
# need to manually change the tax level after ggplot

gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}

ggCols <- gg_color_hue(length(unique(pocdfOther[,taxLevel])))
ggCols <- head(ggCols, n=-1)

pocdfOther$names <- factor(pocdfOther$Sample, levels=rownames(metaFile), ordered = TRUE)

theme_set(theme_bw())
ggplot(pocdfOther, aes(x=names, y=Abundance, fill=Class, order = as.factor(Class))) +
  geom_bar(stat="identity", colour="black") +
  scale_fill_manual(values=c(ggCols, "gray")) +
  scale_y_continuous(expand = c(0,0), limits = c(0,100)) +
  facet_grid(~site, scales='free', space='free_x') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


# SAVE PLOT: EPS 1500 x 800


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

# ordination for the pocillopora samples

theme_set(theme_bw())
pocOrd <- ordinate(poc, "NMDS", "bray")
plot_ordination(poc, pocOrd, type = 'samples', color='pocType', title='poc')

# not many patterns apparent
# subset based on poc type and re-draw

pocType3 <- subset_samples(poc, pocType=='type3')
pocType3Ord <- ordinate(pocType3, "NMDS", "bray")
plot_ordination(pocType3, pocType3Ord, type = 'samples', color='site', title='pocType3', label='names')

pocTypeUnknown <- subset_samples(poc, pocType=='type?')
pocTypeUnknownOrd <- ordinate(pocTypeUnknown, "NMDS", "bray")
plot_ordination(pocTypeUnknown, pocTypeUnknownOrd, type = 'samples', color='site', title='pocTypeUnknown', label='reef')

pocType5 <- subset_samples(poc, pocType=='type5')
pocType5Ord <- ordinate(pocType5, "NMDS", "bray")
plot_ordination(pocType5, pocType5Ord, type = 'samples', color='site', title='pocType5', label='reef')

pocType1 <- subset_samples(poc, pocType=='type1')

pocType5andUnknown <- merge_phyloseq(pocTypeUnknown, pocType3)
pocType5andUnknownOrd <- ordinate(pocType5andUnknown, "NMDS", "bray")
plot_ordination(pocType5andUnknown, pocType5andUnknownOrd, type = 'samples', color='site', title='pocType5andUnknown', label='names')


# check what happens with Endos across poc Types

pocEndoFilt = filter_taxa(pocEndo, function(x) mean(x) > 0.2, TRUE)

theme_set(theme_bw())
plot_bar(pocEndoFilt, fill="catglab", x="names") +
  scale_y_continuous(expand = c(0,0), limits = c(0,100)) +
  scale_fill_brewer(type='qual', palette = 'Set1') +
  facet_grid(~pocType+site, scales='free', space='free_x')


# ordination for the pocillopora verrucosa

theme_set(theme_bw())
pVerrOrd <- ordinate(pVerr, "NMDS", "bray")
plot_ordination(pVerr, pVerrOrd, type = 'samples', color='site', title='pVerr')

# overlay the best Spearman correlations from mothur

pVerrMothurSpearman <- c("MED000007996","MED000007267","MED000007994","MED000007995","MED000007271","MED000007997","MED000007227","MED000006335","MED000003910","MED000000074")


arrowmatrix <- data.frame(labels = pVerrMothurSpearman)
rownames(arrowmatrix) <- arrowmatrix[,1]
arrowdf <- data.frame(arrowmatrix)

# get taxonomic information from the original tax file

arrowdf <- data.frame(labels = allTax[rownames(arrowmatrix),"Order"], arrowmatrix)
arrowdf <- cbind(arrowdf, pVerrOrd$species[pVerrMothurSpearman, ])

# also combine MED node with genus taxonomy for labelling

mylabels = apply(arrowdf[, c("labels", "labels.1")], 1, paste, sep="", collapse="_")
arrowdf <- cbind(arrowdf, catglab=mylabels)

# now create labels, arrows, etc.

arrowmap <- aes(xend = MDS1, yend = MDS2, x = 0, y = 0, shape = NULL, color = NULL, label = labels)
labelmap <- aes(x = MDS1, y = MDS2 + 0.04, shape = NULL, color = NULL, label = catglab, size=1.5)
arrowhead = arrow(length = unit(0.02, "npc"))

pVerrFiltPlot <- plot_ordination(pVerr, pVerrOrd, type = 'samples', color='site', title='pVerr')

pVerrFiltPlotArrow <- pVerrFiltPlot + 
  geom_segment(arrowmap, size = 1, data = arrowdf, color = "grey",  arrow = arrowhead, show_guide = FALSE) +
  scale_color_manual(values=cols) +
  geom_text(labelmap, size = 2, data = arrowdf)

pVerrFiltPlotArrow






