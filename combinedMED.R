## Phyloseq and r analysis of the Minimum Entropy Decomposition results ##
# 15.1.2015

# load required libraries

library("phyloseq")
library("ggplot2")
library("plyr")
library("vegan")
library('ape')
library('RColorBrewer')
library("dunn.test")
setwd("./data")

cols <- c("AmericanSamoa" = "#D95F02", "Indonesia" = "#A6761D", "MaggieIs" = "#666666", "Maldives" = "#E6AB02", "Micronesia" = "#66A61E", "Ningaloo" = "#7570B3", "RedSea" = "#E7298A", "other" = "black")

# import normal percent matrix

allShared = read.table("all.7974.matrixPercent.txt", header=T, row.names=1)

# non-subsampled mothur shared file for alpha diversity

#all3OTUshared = read.table("all.7974.0.03.shared", header=T)
rownames(all3OTUshared) = all3OTUshared[,2]
all3OTUshared = all3OTUshared[,4:length(all3OTUshared)]

# import percent matrix modified for phylogenetic tree

allSharedTree = read.table("all.7974.matrixPercent.tree.txt", header=T)
rownames(allSharedTree) = allSharedTree[,1]
allSharedTree = allSharedTree[,2:length(allSharedTree)]

# Import normal taxonomy file from mothur

allTax = read.table('all.7974.nodeReps.nr_v119.knn.taxonomy', header=T, sep='\t')
rownames(allTax) = allTax[,1]
allTax = allTax[,3:9]
allTax = as.matrix(allTax)

# import meta data (and metaData3 for heatmap ordering)

metaFile = read.table('metaData2.MED', header=T, sep='\t')
rownames(metaFile) = metaFile[,1]
metaFile = metaFile[,2:8]

metaFile3 = read.table('metaData3.txt', header=T, sep='\t')
rownames(metaFile3) = metaFile3[,1]
metaFile3 = metaFile3[,2:8]

# import phylogenetic tree of Endozoicomonas types

endoTreeFile = read.tree(file='MEDNJ5.tree')

### Create phyloseq object

OTU = otu_table(allShared, taxa_are_rows = FALSE)
OTUalpha = otu_table(all3OTUshared, taxa_are_rows = FALSE)
OTUtree = otu_table(allSharedTree, taxa_are_rows = FALSE)
TAX = tax_table(allTax) 
META = sample_data(metaFile)
TREE = phy_tree(endoTreeFile)
allPhylo = phyloseq(OTU, TAX, META)
endoTree = phyloseq(OTUtree, META, TREE)
allAlpha = phyloseq(OTUalpha, META)


# check at how many reefs we sampled corals

allPhyloTmp2 <- subset_samples(allPhylo, species=="Stylophora pistillata")
allPhyloTmp3 <- subset_samples(allPhylo, species=="Pocillopora verrucosa")
allPhyloCorals <- merge_phyloseq(allPhyloTmp2, allPhyloTmp3)

unique(sample_data(allPhyloCorals)$reef)

rownames(otu_table(allPhyloCorals))

# check which samples we included for each coral
#write.table(rownames(otu_table(allPhyloCorals)), "usedSamples.txt", quote=FALSE, row.names=FALSE)
#write.table(sample_data(allPhyloTmp2), "spistMicrobiomeSamples.txt", quote=FALSE, sep="\t")
#write.table(sample_data(allPhyloTmp3), "pverrMicrobiomeSamples.txt", quote=FALSE, sep="\t")

# Alpha diversity measures first

theme_set(theme_bw())
allDivPlot <- plot_richness(allAlpha, x = 'species', measures = c('Chao1', 'Shannon', 'observed'), color = 'site')
allDivPlot

# get rid of a few of those low corals

allAlphaTmp <- subset_samples(allAlpha, species=="seawater")
allAlphaTmp2 <- subset_samples(allAlpha, species=="Stylophora pistillata")
allAlphaTmp3 <- subset_samples(allAlpha, species=="Pocillopora verrucosa")

allAlpha2 <- merge_phyloseq(allAlphaTmp, allAlphaTmp2, allAlphaTmp3)

# The predicted number of OTUs and diversity indecies look similar across the different coral species. I'll add some whisker plots to make this easier to see. 

allAlphaPlot2 <- plot_richness(allAlpha2, x = 'species', measures = c('Chao1', 'Simpson', 'observed'), color = 'site', sortby = 'Chao1')
allAlphaPlot2

ggplot(data = allAlphaPlot2$data) +
  geom_point(aes(x = species, y = value, color = site), position=position_jitter(width=0.1, height=0)) +
  geom_boxplot(aes(x = species, y = value, color = NULL), alpha = 0.1, outlier.shape = NA) +
  scale_color_manual(values=cols) +
  theme(axis.text.x = element_text(angle = 90)) +
  facet_wrap(~variable, scales='free_y') +
  scale_x_discrete(limits=c("Stylophora pistillata", "Pocillopora verrucosa", "seawater"))

## SAVE AS 900 x 600

## check for significant differences using a kruskal-wallis test
# first get mean and se of the diversity esimates in new dataframe

alphaObserved = (estimate_richness(allAlpha2, measures="Observed"))
alphaSimpson = (estimate_richness(allAlpha2, measures="Simpson"))
alphaChao = (estimate_richness(allAlpha2, measures="Chao1"))

alpha.stats <- cbind(alphaObserved, sample_data(allAlpha2))
alpha.stats2 <- cbind(alpha.stats, alphaSimpson)
alpha.stats3 <- cbind(alpha.stats2, alphaChao)

# krustal-wallis test plus a dunn post hoc test to check which groups were different, including bonferroni multiple testing adjustment

kruskal.test(Observed~species, data = alpha.stats3)
dunn.test(alpha.stats3$Observed, alpha.stats3$species, method="bonferroni")

Kruskal-Wallis chi-squared = 83.5999, df = 2, p-value = 0
  Col Mean-|
  Row Mean |   Pocillop   seawater
---------+----------------------
  seawater |   8.453137
           |     0.0000
           |
  Stylopho |   1.307381  -7.797705
           |     0.2866     0.0000

kruskal.test(Simpson~species, data = alpha.stats3)
dunn.test(alpha.stats3$Simpson, alpha.stats3$species, method="bonferroni")

Kruskal-Wallis chi-squared = 18.8184, df = 2, p-value = 8.197e-05
  Col Mean-|
  Row Mean |   Pocillop   seawater
---------+----------------------
  seawater |   4.071168
           |     0.0001
           |
  Stylopho |   0.784873  -3.609790
           |     0.6488     0.0005

kruskal.test(Chao1~species, data = alpha.stats3)
dunn.test(alpha.stats3$Chao1, alpha.stats3$species, method="bonferroni")

Kruskal-Wallis chi-squared = 87.4746, df = 2, p-value < 2.2e-16
  Col Mean-|
  Row Mean |   Pocillop   seawater
---------+----------------------
  seawater |   8.542250
           |     0.0000
           |
  Stylopho |   1.074309  -8.111644
           |     0.4240     0.0000


## quick check on the seriatopora - they have 10-20% endos..
#serio <- subset_samples(allPhylo, species=="Seriatopora sp.")
#serio <- filter_taxa(serio, function(x) mean(x) > 0.1, TRUE)
#plot_bar(serio, fill='Genus')

# some transformations - top line is to keep sample order in bar graphs

sample_data(allPhylo)$names <- factor(sample_names(allPhylo), levels=unique(sample_names(allPhylo)))
allPhylo <- prune_taxa(taxa_sums(allPhylo) > 0, allPhylo)
allPhyloFilt = filter_taxa(allPhylo, function(x) mean(x) > 0.1, TRUE)

# bar chart with everything included

theme_set(theme_bw())
plot_bar(allPhyloFilt, fill="Phylum") +
  scale_y_continuous(expand = c(0,0), limits = c(0,100)) +
  facet_grid(~site, scales='free', space='free_x')

# plot a tree of all endozoicomonas MED nodes

plot_tree(endoTree, nodelabf = nodeplotboot(), ladderize='left', color='species', size='abundance', label.tips = 'taxa_names', base.spacing = 0.005)

# now tree of nodes from the main two corals (plus seawater and the other seqs)

endoTreeSpist <- subset_samples(endoTree, species=='Stylophora pistillata')
endoTreePverr <- subset_samples(endoTree, species=='Pocillopora verrucosa')
endoTreeSea <- subset_samples(endoTree, species=='seawater')
endoTreeOther <- subset_samples(endoTree, species=='other')

endoTreeCorals <- merge_phyloseq(endoTreeSpist, endoTreePverr, endoTreeSea, endoTreeOther)

plot_tree(endoTreeCorals, label.tips = 'taxa_names', color='site', shape='species', size='abundance', nodelabf = nodeplotboot(100, 50, 3), ladderize='left', base.spacing = 0.01) +
  scale_color_manual(values=cols) +
  scale_shape_manual(values = c("other" = 1, "Pocillopora verrucosa" = 17, "Stylophora pistillata" = 16, "seawater" = 15)) 


# SAVE AS EPS 1500 x 1200

theme_set(theme_bw())
plot_bar(endoTree) +
  scale_y_continuous(expand = c(0,0), limits = c(0,100)) +
  facet_grid(~site, scales='free', space='free_x')

# plot a heat map to show these differences a bit better

plot_heatmap(allPhylo, "RDA", "none", "species", "none")

allPhyloEndo = subset_taxa(allPhylo, Genus=='Endozoicomonas')

coralPhyloEndo = subset_samples(allPhyloEndo, species!='seawater')
spistPhyloEndo <- subset_samples(allPhyloEndo, species=='Stylophora pistillata')
pVerrPhyloEndo <- subset_samples(allPhyloEndo, species=='Pocillopora verrucosa')
spistPverrEndo <- merge_phyloseq(spistPhyloEndo, pVerrPhyloEndo)

spistPverrEndoFilt = filter_taxa(spistPverrEndo, function(x) mean(x) > 0.0, TRUE)
spistPverrEndoFiltPrune = prune_samples(sample_sums(spistPverrEndoFilt) > 0, spistPverrEndoFilt)

plot_heatmap(spistPverrEndoFiltPrune, "NMDS", "bray", "site", low="#000033", high="#FF3300", sample.order=rownames(metaFile3))

#plot_heatmap(spistPverrEndoFilt, "RDA", "none", "site", sample.order='species')

# SAVE AS 1500 X 700
# take a quick look at the Archaea samples

archaeaPhylo <- subset_taxa(allPhylo, Domain=="Archaea")
archaeaPhyloPrune = prune_samples(sample_sums(archaeaPhylo) > 0, archaeaPhylo)
plot_heatmap(archaeaPhyloPrune, "NMDS", "bray", "species", sample.order=rownames(metaFile3))

theme_set(theme_bw())
plot_bar(archaeaPhylo, fill="Phylum") +
  scale_y_continuous(expand = c(0,0)) +
  facet_grid(~species, scales='free', space='free_x')

###########################################################################
# use DESeq2 plugin to calculate sig endos 7.5.2015
###########################################################################

library("DESeq2"); packageVersion("DESeq2")

# DESeq requires count data not percentages, so need to import this

allCounts = read.table("all.7974.matrixCount.txt", header=T, row.names=1)
OTUcounts = otu_table(allCounts, taxa_are_rows = FALSE)
countPhylo = phyloseq(OTUcounts, TAX, META)

# subset and filter endozoicomonas OTUs

countEndos = subset_taxa(countPhylo, Genus=='Endozoicomonas')

spistCountEndo <- subset_samples(countEndos, species=='Stylophora pistillata')
pVerrCountEndo <- subset_samples(countEndos, species=='Pocillopora verrucosa')
coralCountEndo <- merge_phyloseq(spistCountEndo, pVerrCountEndo)

testSpecies <- spistCountEndo

coralCountEndoFilt = filter_taxa(testSpecies, function(x) mean(x) > 0.0, TRUE)
coralCountEndoFiltPrune = prune_samples(sample_sums(coralCountEndoFilt) > 0, coralCountEndoFilt)

# first convert phyloseq object to DESeq object, then run Deseq

endoDeseq <- phyloseq_to_deseq2(coralCountEndoFiltPrune, ~ site)

# need to calculate geometric means separately because there are zeros in the data

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

geomMeans <- apply(counts(endoDeseq), 1, gm_mean)
endoDeseq <- estimateSizeFactors(endoDeseq, geoMeans = geomMeans)
 
# now can run the DESeq tests

endoDeseq <- DESeq(endoDeseq, fitType="local")

# check the results for spistillata

resultsNames(endoDeseq)

res <- results(endoDeseq, contrast = c("site", "AmericanSamoa", "RedSea"))
res = res[order(res$padj, na.last=NA), ]
sigtab = as.data.frame(res[(res$padj < 0.05), ])
sigtab
res <- results(endoDeseq, contrast = c("site", "AmericanSamoa", "Ningaloo"))
res = res[order(res$padj, na.last=NA), ]
sigtab = as.data.frame(res[(res$padj < 0.05), ])
sigtab
res <- results(endoDeseq, contrast = c("site", "AmericanSamoa", "Micronesia"))
res = res[order(res$padj, na.last=NA), ]
sigtab = as.data.frame(res[(res$padj < 0.05), ])
sigtab 
res <- results(endoDeseq, contrast = c("site", "AmericanSamoa", "Indonesia"))
res = res[order(res$padj, na.last=NA), ]
sigtab = as.data.frame(res[(res$padj < 0.05), ])
sigtab 
res <- results(endoDeseq, contrast = c("site", "RedSea", "Ningaloo"))
res = res[order(res$padj, na.last=NA), ]
sigtab = as.data.frame(res[(res$padj < 0.05), ])
sigtab 
res <- results(endoDeseq, contrast = c("site", "RedSea", "Micronesia"))
res = res[order(res$padj, na.last=NA), ]
sigtab = as.data.frame(res[(res$padj < 0.05), ])
sigtab 
res <- results(endoDeseq, contrast = c("site", "RedSea", "Indonesia"))
res = res[order(res$padj, na.last=NA), ]
sigtab = as.data.frame(res[(res$padj < 0.05), ])
sigtab 
res <- results(endoDeseq, contrast = c("site", "Ningaloo", "Micronesia"))
res = res[order(res$padj, na.last=NA), ]
sigtab = as.data.frame(res[(res$padj < 0.05), ])
sigtab 
res <- results(endoDeseq, contrast = c("site", "Ningaloo", "Indonesia"))
res = res[order(res$padj, na.last=NA), ]
sigtab = as.data.frame(res[(res$padj < 0.05), ])
sigtab 
res <- results(endoDeseq, contrast = c("site", "Micronesia", "Indonesia"))
res = res[order(res$padj, na.last=NA), ]
sigtab = as.data.frame(res[(res$padj < 0.05), ])
sigtab 


## OK, now run these test for p.verrucosa

testSpecies <- pVerrCountEndo

coralCountEndoFilt = filter_taxa(testSpecies, function(x) mean(x) > 0.0, TRUE)
coralCountEndoFiltPrune = prune_samples(sample_sums(coralCountEndoFilt) > 0, coralCountEndoFilt)

# first convert phyloseq object to DESeq object, then run Deseq

endoDeseq <- phyloseq_to_deseq2(coralCountEndoFiltPrune, ~ site)
geomMeans <- apply(counts(endoDeseq), 1, gm_mean)
endoDeseq <- estimateSizeFactors(endoDeseq, geoMeans = geomMeans)

# now can run the DESeq tests

endoDeseq <- DESeq(endoDeseq, fitType="local")

res <- results(endoDeseq, contrast = c("site", "Maldives", "RedSea"))
res = res[order(res$padj, na.last=NA), ]
sigtab = res[(res$padj < 0.05), ]
sigtab 
res <- results(endoDeseq, contrast = c("site", "Maldives", "Micronesia"))
res = res[order(res$padj, na.last=NA), ]
sigtab = res[(res$padj < 0.05), ]
sigtab 
res <- results(endoDeseq, contrast = c("site", "Maldives", "Indonesia"))
res = res[order(res$padj, na.last=NA), ]
sigtab = res[(res$padj < 0.05), ]
sigtab 
res <- results(endoDeseq, contrast = c("site", "RedSea", "Micronesia"))
res = res[order(res$padj, na.last=NA), ]
sigtab = res[(res$padj < 0.05), ]
sigtab 
res <- results(endoDeseq, contrast = c("site", "RedSea", "Indonesia"))
res = res[order(res$padj, na.last=NA), ]
sigtab = res[(res$padj < 0.05), ]
sigtab 
res <- results(endoDeseq, contrast = c("site", "Micronesia", "Indonesia"))
res = res[order(res$padj, na.last=NA), ]
sigtab = res[(res$padj < 0.05), ]
sigtab 


###########################################################################
# see if SIMPER works 5.5.2015
###########################################################################

spistEndoSIMP <- subset_samples(spistPverrEndoFiltPrune, species=="Stylophora pistillata")
spistEndoSIMPFilt = filter_taxa(spistEndoSIMP, function(x) mean(x) > 0.0, TRUE)
spistEndoSIMPFiltPrune = prune_samples(sample_sums(spistEndoSIMPFilt) > 0, spistEndoSIMPFilt)

pverrEndoSIMP <- subset_samples(spistPverrEndoFiltPrune, species=="Pocillopora verrucosa")
pverrEndoSIMPFilt = filter_taxa(pverrEndoSIMP, function(x) mean(x) > 0.0, TRUE)
pverrEndoSIMPFiltPrune = prune_samples(sample_sums(pverrEndoSIMPFilt) > 0, pverrEndoSIMPFilt)

metaFileSIMP <- subset(metaFile, site != 'NA')

SimpSpistEndo <- merge(otu_table(spistEndoSIMPFiltPrune), metaFileSIMP, by="row.names")
rownames(SimpSpistEndo) <- SimpSpistEndo[,1]
SimpSpistEndoSIMP <- SimpSpistEndo[c(1:6,8:71),2:46]
factorsSimpSpistEndo <- SimpSpistEndo[c(1:6,8:71),47:53]

spistSIMP <- simper(SimpSpistEndoSIMP, factorsSimpSpistEndo$site)


SimppverrEndo <- merge(otu_table(pverrEndoSIMPFiltPrune), metaFileSIMP, by="row.names")
rownames(SimppverrEndo) <- SimppverrEndo[,1]
SimppverrEndoSIMP <- SimppverrEndo[c(1:18,20:50),2:38]
factorsSimppverrEndo <- SimppverrEndo[c(1:18,20:50),39:45]

pverrSIMP <- simper(SimppverrEndoSIMP, factorsSimppverrEndo$site)



