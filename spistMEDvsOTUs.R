
## Phyloseq and r analysis of the Minimum Entropy Decomposition results for spist ##
# 19.1.2015

# load required libraries

library("phyloseq")
library("ggplot2")
library("plyr")
library("vegan")
library("RColorBrewer")
setwd("./data")

# generate some colors to be consistent

#display.brewer.all()
#display.brewer.pal(n = 8, name = 'Dark2')
#brewer.pal(n = 8, name = "Dark2")

cols <- c("AmericanSamoa" = "#D95F02", "Indonesia" = "#A6761D", "MaggieIs" = "#666666", "Maldives" = "#E6AB02", "Micronesia" = "#66A61E", "Ningaloo" = "#7570B3", "RedSea" = "#E7298A")

# import normal percent matrix

allMED = read.table("all.7974.matrixPercent.txt", header=T)
rownames(allMED) = allMED[,1]
allMED = allMED[,2:length(allMED)]

# Import normal taxonomy file from mothur

allMEDTax = read.table('all.7974.nodeReps.nr_v119.knn.taxonomy', header=T, sep='\t')
rownames(allMEDTax) = allMEDTax[,1]
allMEDTax = allMEDTax[,3:9]
allMEDTax = as.matrix(allMEDTax)

# import meta data

metaFile = read.table('metaData2.MED', header=T, sep='\t')
rownames(metaFile) = metaFile[,1]
metaFile = metaFile[,2:7]

# 3% OTU and 1% OTU data for comparison

all3OTUshared = read.table("all.7974.0.03.pick.shared", header=T)
rownames(all3OTUshared) = all3OTUshared[,2]
all3OTUshared = all3OTUshared[,4:length(all3OTUshared)]

all1OTUshared = read.table("all.7974.0.01.pick.shared", header=T)
rownames(all1OTUshared) = all1OTUshared[,2]
all1OTUshared = all1OTUshared[,4:length(all1OTUshared)]

# Import taxonomy file from mothur

all3OTUtax = read.table('all.7974.0.03.taxonomy', header=T, sep='\t')
rownames(all3OTUtax) = all3OTUtax[,1]
all3OTUtax = all3OTUtax[,3:9]
all3OTUtax = as.matrix(all3OTUtax)

all1OTUtax = read.table('all.7974.0.01.taxonomy', header=T, sep='\t')
rownames(all1OTUtax) = all1OTUtax[,1]
all1OTUtax = all1OTUtax[,3:9]
all1OTUtax = as.matrix(all1OTUtax)

### Create phyloseq object

OTU = otu_table(allMED, taxa_are_rows = FALSE)
TAX = tax_table(allMEDTax) 
META = sample_data(metaFile)
allPhylo = phyloseq(OTU, TAX, META)

OTUs3 = otu_table(all3OTUshared, taxa_are_rows = FALSE)
OTUs1 = otu_table(all1OTUshared, taxa_are_rows = FALSE)

TAX3 = tax_table(all3OTUtax)
TAX1 = tax_table(all1OTUtax)

all3OTUphylo = phyloseq(OTUs3, TAX3, META)
all1OTUphylo = phyloseq(OTUs1, TAX1, META)

spistPhylo <- subset_samples(allPhylo, species=='Stylophora pistillata')
spistPhylo = filter_taxa(spistPhylo, function(x) mean(x) > 0, TRUE)

spist3OTUphylo <- subset_samples(all3OTUphylo, species=='Stylophora pistillata')
spist3OTUphylo = filter_taxa(spist3OTUphylo, function(x) mean(x) > 0, TRUE)

spist1OTUphylo <- subset_samples(all1OTUphylo, species=='Stylophora pistillata')
spist1OTUphylo = filter_taxa(spist1OTUphylo, function(x) mean(x) > 0, TRUE)

# make mothur files into relative abundance

spist3OTUphyloRel = transform_sample_counts(spist3OTUphylo, function(x) x / sum(x) )
spist1OTUphyloRel = transform_sample_counts(spist1OTUphylo, function(x) x / sum(x) )

# ordination comparing MED nodes and 3% OTUs

theme_set(theme_bw())
spistMEDord <- ordinate(spistPhylo, "NMDS", "bray")
plot_ordination(spistPhylo, spistMEDord, type = 'samples', color='site', title='MED nodes') +
  geom_point(size=2) +
  scale_color_manual(values=cols)

# SAVE AS 700 x 532
# Best solution 0.24, 762 nodes

spist3OTUord <- ordinate(spist3OTUphylo, "NMDS", "bray")
plot_ordination(spist3OTUphylo, spist3OTUord, type = 'samples', color='site', title='3% OTUs') +
  geom_point(size=2) +
  scale_color_manual(values=cols)

# Best solution 0.24, 680 OTUs

spist1OTUord <- ordinate(spist1OTUphylo, "NMDS", "bray")
plot_ordination(spist1OTUphylo, spist1OTUord, type = 'samples', color='site', title='1% OTUs') +
  geom_point(size=2) +
  scale_color_manual(values=cols)

# SAVE AS 700 x 532
# Best solution 0.22, 833 OTUs

# do some bar plots to check if MED nodes resolve endozoic types better than OTUs

spistPhyloEndo = subset_taxa(spistPhylo, Genus=='Endozoicomonas')
spist3OTUphyloEndo = subset_taxa(spist3OTUphyloRel, Genus=='Endozoicomonas(100)')
spist1OTUphyloEndo = subset_taxa(spist1OTUphyloRel, Genus=='Endozoicomonas(100)')

spistPhyloEndoBar <- plot_bar(spistPhyloEndo, fill="Genus")
spistPhyloEndoBar + facet_wrap(~site, scales='free')

spistOTUphyloEndoBar <- plot_bar(spist3OTUphyloEndo, fill="Genus")
spistOTUphyloEndoBar + facet_wrap(~site, scales='free')

spist1OTUphyloEndoBar <- plot_bar(spist1OTUphyloEndo, fill="Genus")
spist1OTUphyloEndoBar + facet_wrap(~site, scales='free')

# add coloring for different Endozoicomonas OTUs

tax_table(spistPhyloEndo) <- cbind(tax_table(spistPhyloEndo), Strain=taxa_names(spistPhyloEndo))
tax_table(spist3OTUphyloEndo) <- cbind(tax_table(spist3OTUphyloEndo), Strain=taxa_names(spist3OTUphyloEndo))
tax_table(spist1OTUphyloEndo) <- cbind(tax_table(spist1OTUphyloEndo), Strain=taxa_names(spist1OTUphyloEndo))
myranks = c("Genus", "Strain")

mylabelsMED = apply(tax_table(spistPhyloEndo)[, myranks], 1, paste, sep="", collapse="_")
mylabels3OTU = apply(tax_table(spist3OTUphyloEndo)[, myranks], 1, paste, sep="", collapse="_")
mylabels1OTU = apply(tax_table(spist1OTUphyloEndo)[, myranks], 1, paste, sep="", collapse="_")

tax_table(spistPhyloEndo) <- cbind(tax_table(spistPhyloEndo), catglab=mylabelsMED)
tax_table(spist3OTUphyloEndo) <- cbind(tax_table(spist3OTUphyloEndo), catglab=mylabels3OTU)
tax_table(spist1OTUphyloEndo) <- cbind(tax_table(spist1OTUphyloEndo), catglab=mylabels1OTU)

# standardize the OTU shared files

spist3OTUphyloEndoFilt = filter_taxa(spist3OTUphyloEndo, function(x) mean(x) > 0, TRUE)
plot_bar(spist3OTUphyloEndoFilt, fill="catglab") +
facet_wrap(~site, scales='free')

spist1OTUphyloEndoFilt = filter_taxa(spist1OTUphyloEndo, function(x) mean(x) > 0, TRUE)
plot_bar(spist1OTUphyloEndoFilt, fill="catglab") +
  facet_wrap(~site, scales='free')

spistPhyloEndoFilt = filter_taxa(spistPhyloEndo, function(x) mean(x) > 0.1, TRUE)
plot_bar(spistPhyloEndoFilt, fill="catglab") +
  facet_wrap(~site, scales='free')

# let's have a look at what happends to the most abundant Endo OTU

# 3% OTUs - OTU0001

#spist3OTUphyloEndoRel = transform_sample_counts(spist3OTUphyloEndo, function(x) x / sum(x) )

spist3OTUphyloEndo1 = subset_taxa(spist3OTUphyloEndo, catglab=='Endozoicomonas(100)_Otu00001')
end1bar <- plot_bar(spist3OTUphyloEndo1, title='3% OTUs')

ggplot(end1bar$data, aes(x=site, y=Abundance)) +
  geom_boxplot(aes(fill=catglab), alpha=0.8) +
  geom_point(position=position_dodge(width=0.75), aes(group=catglab)) 

# just Micronesia vs Red Sea

spist3OTUphyloEndo1Mic = subset_samples(spist3OTUphyloEndo1, site=='Micronesia')
spist3OTUphyloEndo1RS = subset_samples(spist3OTUphyloEndo1, site=='RedSea')
spist3OTUphyloEndo1MicRS = merge_phyloseq(spist3OTUphyloEndo1Mic, spist3OTUphyloEndo1RS)

spist3OTUphyloEndo1MicRSBar <- plot_bar(spist3OTUphyloEndo1MicRS, title='3% OTUs')

ggplot(spist3OTUphyloEndo1MicRSBar$data, aes(x=site, y=Abundance, title='OTU1')) +
  geom_boxplot(aes(fill = site)) +
  geom_point(position=position_jitter(width=0.2), aes(group=catglab, fill=site), size=5, pch=21, colour = 'black') +
  scale_fill_manual(values=cols) +
  scale_color_manual(values=cols) 

## SAVE AS 492 x 450

# 1% OTUs first, second and third interesting ones

#spist1OTUphyloEndoRel = transform_sample_counts(spist1OTUphyloEndo, function(x) x / sum(x) )
spistOTUphyloEndo1_1 = subset_taxa(spist1OTUphyloEndo, catglab=="Endozoicomonas(100)_Otu00007")
spistOTUphyloEndo1_2 = subset_taxa(spist1OTUphyloEndo, catglab=="Endozoicomonas(100)_Otu00010")
spistOTUphyloEndo1_3 = subset_taxa(spist1OTUphyloEndo, catglab=="Endozoicomonas(100)_Otu00013")


spistOTUphyloEndo1_all3 = merge_phyloseq(spistOTUphyloEndo1_1, spistOTUphyloEndo1_2, spistOTUphyloEndo1_3)

end1bar_1 <- plot_bar(spistOTUphyloEndo1_all3, title='1% OTUs', fill='catglab')

ggplot(end1bar_1$data, aes(x=site, y=Abundance)) +
  geom_boxplot(aes(fill=catglab)) +
  geom_point(position=position_dodge(width=0.75), aes(group=catglab)) 

# just Micronesia vs Red Sea

spistOTUphyloEndo1_1Mic = subset_samples(spistOTUphyloEndo1_all3, site=='Micronesia')
spistOTUphyloEndo1_1RS = subset_samples(spistOTUphyloEndo1_all3, site=='RedSea')
spistOTUphyloEndo1_1MicRS = merge_phyloseq(spistOTUphyloEndo1_1Mic, spistOTUphyloEndo1_1RS)

spistOTUphyloEndo1_1MicRSBar <- plot_bar(spistOTUphyloEndo1_1MicRS, title='1% OTUs')

ggplot(spistOTUphyloEndo1_1MicRSBar$data, aes(x=site, y=Abundance)) +
  geom_boxplot(aes(fill=catglab)) +
  geom_point(position=position_dodge(width=0.75), aes(group=catglab)) 

## SAVE AS 492 x 450

# 1st OTU

spistOTUphyloEndo1_1Mic = subset_samples(spistOTUphyloEndo1_1, site=='Micronesia')
spistOTUphyloEndo1_1RS = subset_samples(spistOTUphyloEndo1_1, site=='RedSea')
spistOTUphyloEndo1_1MicRS = merge_phyloseq(spistOTUphyloEndo1_1RS, spistOTUphyloEndo1_1Mic)

spistOTUphyloEndo1_1MicRSBar <- plot_bar(spistOTUphyloEndo1_1MicRS, title='OTU1')

ggplot(spistOTUphyloEndo1_1MicRSBar$data, aes(x=site, y=Abundance, title='OTU7')) +
  geom_boxplot(aes(fill = site)) +
  geom_point(position=position_jitter(width=0.2), aes(group=catglab, fill=site), size=5, pch=21, colour = 'black') +
  scale_fill_manual(values=cols) +
  scale_color_manual(values=cols) 

## SAVE AS 492 x 450

# 2nd OTU

spistOTUphyloEndo1_2Mic = subset_samples(spistOTUphyloEndo1_2, site=='Micronesia')
spistOTUphyloEndo1_2RS = subset_samples(spistOTUphyloEndo1_2, site=='RedSea')
spistOTUphyloEndo1_2MicRS = merge_phyloseq(spistOTUphyloEndo1_2RS, spistOTUphyloEndo1_2Mic)

spistOTUphyloEndo1_2MicRSBar <- plot_bar(spistOTUphyloEndo1_2MicRS)

ggplot(spistOTUphyloEndo1_2MicRSBar$data, aes(x=site, y=Abundance, title='OTU10')) +
  geom_boxplot(aes(fill = site)) +
  geom_point(position=position_jitter(width=0.2), aes(group=catglab, fill=site), size=5, pch=21, colour = 'black') +
  scale_fill_manual(values=cols) +
  scale_color_manual(values=cols) 

## SAVE AS 492 x 450

# 3rd OTU

spistOTUphyloEndo1_3Mic = subset_samples(spistOTUphyloEndo1_3, site=='Micronesia')
spistOTUphyloEndo1_3RS = subset_samples(spistOTUphyloEndo1_3, site=='RedSea')
spistOTUphyloEndo1_3MicRS = merge_phyloseq(spistOTUphyloEndo1_3RS, spistOTUphyloEndo1_3Mic)

spistOTUphyloEndo1_3MicRSBar <- plot_bar(spistOTUphyloEndo1_3MicRS)

ggplot(spistOTUphyloEndo1_3MicRSBar$data, aes(x=site, y=Abundance, title='OTU13')) +
  geom_boxplot(aes(fill = site)) +
  geom_point(position=position_jitter(width=0.2), aes(group=catglab, fill=site), size=5, pch=21, colour = 'black') +
  scale_fill_manual(values=cols) +
  scale_color_manual(values=cols) 


# MED nodes split from the original Otu0001 at 3%

#spistMEDEndoRel = transform_sample_counts(spistPhyloEndo, function(x) x / sum(x) )

spistMEDEndo_1 = subset_taxa(spistPhyloEndo, catglab=="Endozoicomonas_MED000005672")
spistMEDEndo_2 = subset_taxa(spistPhyloEndo, catglab=="Endozoicomonas_MED000005693")
spistMEDEndo_3 = subset_taxa(spistPhyloEndo, catglab=="Endozoicomonas_MED000008661")
#spistMEDEndoRel1_4669 = subset_taxa(spistMEDEndoRel, catglab=="Endozoicomonas_MED000004669")
#spistMEDEndoRel1_4794 = subset_taxa(spistMEDEndoRel, catglab=="Endozoicomonas_MED000004794")

spistMEDEndo_all3 = merge_phyloseq(spistMEDEndo_1, spistMEDEndo_2, spistMEDEndo_3)

endMEDbar_1 <- plot_bar(spistMEDEndo_all3, title='MED nodes', fill='catglab')

ggplot(endMEDbar_1$data, aes(x=site, y=Abundance)) +
  geom_boxplot(aes(fill=catglab), alpha=0.8) +
  geom_point(position=position_dodge(width=0.75), aes(group=catglab)) 

# just Micronesia vs Red Sea

spistMEDEndo_all3Mic = subset_samples(spistMEDEndo_all3, site=='Micronesia')
spistMEDEndo_all3RS = subset_samples(spistMEDEndo_all3, site=='RedSea')
spistMEDEndo_all3MicRS = merge_phyloseq(spistMEDEndo_all3Mic, spistMEDEndo_all3RS)

spistMEDEndo_all3MicRSBar <- plot_bar(spistMEDEndo_all3MicRS, title='MED nodes')

ggplot(spistMEDEndo_all3MicRSBar$data, aes(x=site, y=Abundance)) +
  geom_boxplot(aes(fill=catglab), alpha=0.8) +
  geom_point(position=position_dodge(width=0.75), aes(group=catglab)) 

# 1st MED node

spistMEDEndoRel1_1Mic = subset_samples(spistMEDEndo_1, site=='Micronesia')
spistMEDEndoRel1_1RS = subset_samples(spistMEDEndo_1, site=='RedSea')
spistMEDEndoRel1_1MicRS = merge_phyloseq(spistMEDEndoRel1_1RS, spistMEDEndoRel1_1Mic)

spistMEDEndoRel1_1MicRSBar <- plot_bar(spistMEDEndoRel1_1MicRS, title='MED nodes')

ggplot(spistMEDEndoRel1_1MicRSBar$data, aes(x=site, y=Abundance, title='MED5672')) +
  geom_boxplot(aes(fill = site)) +
  geom_point(position=position_jitter(width=0.2), aes(group=catglab, fill=site), size=5, pch=21, colour = 'black') +
  scale_fill_manual(values=cols) +
  scale_color_manual(values=cols)

## SAVE AS 492 x 450

# 2nd MED node

spistMEDEndoRel1_2Mic = subset_samples(spistMEDEndo_2, site=='Micronesia')
spistMEDEndoRel1_2RS = subset_samples(spistMEDEndo_2, site=='RedSea')
spistMEDEndoRel1_2MicRS = merge_phyloseq(spistMEDEndoRel1_2RS, spistMEDEndoRel1_2Mic)

spistMEDEndoRel1_2MicRSBar <- plot_bar(spistMEDEndoRel1_2MicRS, title='MED nodes')

ggplot(spistMEDEndoRel1_2MicRSBar$data, aes(x=site, y=Abundance, title='MED5693')) +
  geom_boxplot(aes(fill = site)) +
  geom_point(position=position_jitter(width=0.2), aes(group=catglab, fill=site), size=5, pch=21, colour = 'black') +
  scale_fill_manual(values=cols) +
  scale_color_manual(values=cols) 

## SAVE AS 492 x 450

# 3rd MED node

spistMEDEndoRel1_3Mic = subset_samples(spistMEDEndo_3, site=='Micronesia')
spistMEDEndoRel1_3RS = subset_samples(spistMEDEndo_3, site=='RedSea')
spistMEDEndoRel1_3MicRS = merge_phyloseq(spistMEDEndoRel1_3RS, spistMEDEndoRel1_3Mic)

spistMEDEndoRel1_3MicRSBar <- plot_bar(spistMEDEndoRel1_3MicRS, title='MED nodes')

ggplot(spistMEDEndoRel1_3MicRSBar$data, aes(x=site, y=Abundance, title='MED8661')) +
  geom_boxplot(aes(fill = site)) +
  geom_point(position=position_jitter(width=0.2), aes(group=catglab, fill=site), size=5, pch=21, colour = 'black') +
  scale_fill_manual(values=cols) +
  scale_color_manual(values=cols) 

