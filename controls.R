
## controls from MED analysis ##

# 21.4.15


library("phyloseq")
library("ggplot2")
library('RColorBrewer')

setwd("./data")


# import percent matrix

controlsShared = read.table("controls.matrixPercent.txt", header=T)
rownames(controlsShared) = controlsShared[,1]
controlsShared = controlsShared[,2:length(controlsShared)]

# Import taxonomy file from mothur

controlsTax = read.table('controls.NodeReps.nr_v119.knn.taxonomy', header=T, sep='\t')
rownames(controlsTax) = controlsTax[,1]
controlsTax = controlsTax[,3:9]
controlsTax = as.matrix(controlsTax)

# import meta data 

metaFile = read.table('metaData2.MED', header=T, sep='\t')
rownames(metaFile) = metaFile[,1]
metaFile = metaFile[,2:8]

# Create phyloseq object

OTU = otu_table(controlsShared, taxa_are_rows = FALSE)
TAX = tax_table(controlsTax) 
META = sample_data(metaFile)
controlsPhylo = phyloseq(OTU, TAX, META)


# Compare samples that were replicated in each MiSeq run

replicatedSamples = prune_samples(c('RP7-rep2', 'RP8-rep2', 'RS5-rep2', 'END74-rep2', 'END59-rep2', 'END62-rep2', '55-2-rep2', '56-2-rep2', '103-rep2', '105-rep2', '250-rep2', '253-rep2', '254-rep2', 'MIC-53-rep2', 'MIC-85-rep2', 'RP7', 'RP8', 'RS5', 'END74', 'END59', 'END62', '55-2', '56-2', '103', '105', '250', '253', '254', 'MIC-53', 'MIC-85'), controlsPhylo)

replicatedSamplesFilt = filter_taxa(replicatedSamples, function(x) mean(x) > 0.2, TRUE)
replicatedSamplesFiltGlom =tax_glom(replicatedSamplesFilt, taxrank="Genus")

theme_set(theme_bw())
plot_bar(replicatedSamplesFiltGlom, fill="Genus", title='Replicated Samples') +
  scale_y_continuous(expand = c(0,0), limits = c(0,100))

# Mock communities

# Barplots of the mock community compositions

# First the even mock community from pool1. Red dotted line is the expected abundance of each of the taxa. 

evP1 <- subset_samples(controlsPhylo, reef=='even-mock-pool1')
evP1filt = filter_taxa(evP1, function(x) mean(x) > 0.0001, TRUE)
evP1df <- psmelt(evP1filt)

evP1df[order(rownames(evP1df)),]

evP1df$relAbund = evP1df$Abundance / sum(evP1df$Abundance)
evP1df$otuNumber <- rownames(evP1df)

evP1df$otuNumber <- factor(evP1df$otuNumber, levels=evP1df$otuNumber, ordered=TRUE)

barMockPlot <- ggplot(data=evP1df) + 
  geom_bar(aes(x=otuNumber, y=relAbund, fill=Genus), stat="identity") +
  theme(axis.text.x=element_text(angle=90)) +
  geom_hline(aes(yintercept=.047), colour="#990000", linetype="dashed")
  #scale_x_discrete(limits=evP1df$otuNumber)
barMockPlot

# There are 21 OTUs expected from this mock community, so a couple are missing. It seems as though the large pink bar (Staphylococcus) actually contains 2 Staphylococcus OTUs, which accounts for one of the missing members. The other 2 missing are a Methanobrevibacter and a Propionibacterium - apparently these are not amplified with our primers. 

# The very small bars to the right must be the spurious OTUs generated from errors that were coming up in the rarefaction curves. 







