
## Phyloseq and r analysis of the Minimum Entropy Decomposition results for spist ##
# 19.1.2015

# load required libraries

library("phyloseq")
library("ggplot2")
?library("plyr")
library("vegan")
setwd("./data")

# import normal percent matrix

spistMED = read.table("spist.7801.MED.matrixPercent", header=T)
rownames(spistMED) = spistMED[,1]
spistMED = spistMED[,2:length(spistMED)]

# Import normal taxonomy file from mothur

spistMEDTax = read.table('spist.7801.MED.nodeReps.nr_v119.knn.taxonomy', header=T, sep='\t')
rownames(spistMEDTax) = spistMEDTax[,1]
spistMEDTax = spistMEDTax[,3:9]
spistMEDTax = as.matrix(spistMEDTax)

# import meta data

metaFile = read.table('metaData2.MED', header=T, sep='\t')
rownames(metaFile) = metaFile[,1]
metaFile = metaFile[,2:7]

# 3% OTU and 1% OTU data for comparison

spist3OTUshared = read.table("spist.7801.0.03.0.03.pick.shared", header=T)
rownames(spist3OTUshared) = spist3OTUshared[,2]
spist3OTUshared = spist3OTUshared[,4:length(spist3OTUshared)]

spist1OTUshared = read.table("spist.7801.0.01.0.01.pick.shared", header=T)
rownames(spist1OTUshared) = spist1OTUshared[,2]
spist1OTUshared = spist1OTUshared[,4:length(spist1OTUshared)]

# Import taxonomy file from mothur

spist3OTUtax = read.table('spist.7801.unique.good.filter.precluster.an.0.03.cons.taxonomy', header=T, sep='\t')
rownames(spist3OTUtax) = spist3OTUtax[,1]
spist3OTUtax = spist3OTUtax[,3:9]
spist3OTUtax = as.matrix(spist3OTUtax)

spist1OTUtax = read.table('spist.7801.unique.good.filter.precluster.an.0.01.cons.taxonomy', header=T, sep='\t')
rownames(spist1OTUtax) = spist1OTUtax[,1]
spist1OTUtax = spist1OTUtax[,3:9]
spist1OTUtax = as.matrix(spist1OTUtax)

### Create phyloseq object

OTU = otu_table(spistMED, taxa_are_rows = FALSE)
TAX = tax_table(spistMEDTax) 
META = sample_data(metaFile)
spistPhylo = phyloseq(OTU, TAX, META)

OTUs3 = otu_table(spist3OTUshared, taxa_are_rows = FALSE)
OTUs1 = otu_table(spist1OTUshared, taxa_are_rows = FALSE)

TAX3 = tax_table(spist3OTUtax)
TAX1 = tax_table(spist1OTUtax)

spist3OTUphylo = phyloseq(OTUs3, TAX3, META)
spist1OTUphylo = phyloseq(OTUs1, TAX1, META)


# ordination comparing MED nodes and 3% OTUs

theme_set(theme_bw())
spistMEDord <- ordinate(spistPhylo, "NMDS", "bray")
plot_ordination(spistPhylo, spistMEDord, type = 'samples', color='site', title='MED nodes')

spist3OTUord <- ordinate(spist3OTUphylo, "NMDS", "bray")
plot_ordination(spist3OTUphylo, spist3OTUord, type = 'samples', color='site', title='3% OTUs')

spist1OTUord <- ordinate(spist1OTUphylo, "NMDS", "bray")
plot_ordination(spist1OTUphylo, spist1OTUord, type = 'samples', color='site', title='1% OTUs')

# do some bar plots to check if MED nodes resolve endozoic types better than OTUs

spistPhyloEndo = subset_taxa(spistPhylo, Genus=='Endozoicomonas')
spist3OTUphyloEndo = subset_taxa(spist3OTUphylo, Genus=='Endozoicomonas(100)')
spist1OTUphyloEndo = subset_taxa(spist1OTUphylo, Genus=='Endozoicomonas(100)')

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

spist3OTUphyloEndoMerged = merge_samples(spist3OTUphyloEndo, "site")
spist3OTUphyloEndoMergedRel = transform_sample_counts(spist3OTUphyloEndoMerged, function(x) x / sum(x) )
plot_bar(spist3OTUphyloEndoMergedRel, fill='catglab', title='3% OTUs')

spist1OTUphyloEndoMerged = merge_samples(spist1OTUphyloEndo, "site")
spist1OTUphyloEndoMergedRel = transform_sample_counts(spist1OTUphyloEndoMerged, function(x) x / sum(x) )
spist1OTUphyloEndoMergedRelFilt = filter_taxa(spist1OTUphyloEndoMergedRel, function(x) mean(x) > 0.01, TRUE)
plot_bar(spist1OTUphyloEndoMergedRelFilt, fill='catglab', title='1% OTUs')

spistPhyloEndoMerged = merge_samples(spistPhyloEndo, "site")
spistPhyloEndoMergedRel = transform_sample_counts(spistPhyloEndoMerged, function(x) x / sum(x) )
spistPhyloEndoMergedRelFilt = filter_taxa(spistPhyloEndoMergedRel, function(x) mean(x) > 0.01, TRUE)
plot_bar(spistPhyloEndoMergedRelFilt, fill='catglab', title='MED otus')

# let's have a look at what happends to the most abundant Endo OTU
# 3% OTUs - OTU0001

spist3OTUphyloEndoRel = transform_sample_counts(spist3OTUphyloEndo, function(x) x / sum(x) )
spist3OTUphyloEndo1 = subset_taxa(spist3OTUphyloEndoRel, catglab=='Endozoicomonas(100)_Otu0001')
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
  scale_fill_manual(values=c("#00BA38", "#F8766D")) +
  geom_point(position=position_jitter(width=0.2), aes(group=catglab), size=5) 

# 1% OTUs

spist1OTUphyloEndoRel = transform_sample_counts(spist1OTUphyloEndo, function(x) x / sum(x) )
spistOTUphyloEndo1_3 = subset_taxa(spist1OTUphyloEndoRel, catglab=="Endozoicomonas(100)_Otu00003")
spistOTUphyloEndo1_6 = subset_taxa(spist1OTUphyloEndoRel, catglab=="Endozoicomonas(100)_Otu00006")
spistOTUphyloEndo1_7 = subset_taxa(spist1OTUphyloEndoRel, catglab=="Endozoicomonas(100)_Otu00007")


spistOTUphyloEndo1_1 = merge_phyloseq(spistOTUphyloEndo1_3, spistOTUphyloEndo1_6, spistOTUphyloEndo1_7)

end1bar_1 <- plot_bar(spistOTUphyloEndo1_1, title='1% OTUs', fill='catglab')

ggplot(end1bar_1$data, aes(x=site, y=Abundance)) +
  geom_boxplot(aes(fill=catglab)) +
  geom_point(position=position_dodge(width=0.75), aes(group=catglab)) 

# just Micronesia vs Red Sea

spistOTUphyloEndo1_1Mic = subset_samples(spistOTUphyloEndo1_1, site=='Micronesia')
spistOTUphyloEndo1_1RS = subset_samples(spistOTUphyloEndo1_1, site=='RedSea')
spistOTUphyloEndo1_1MicRS = merge_phyloseq(spistOTUphyloEndo1_1Mic, spistOTUphyloEndo1_1RS)

spistOTUphyloEndo1_1MicRSBar <- plot_bar(spistOTUphyloEndo1_1MicRS, title='1% OTUs')

ggplot(spistOTUphyloEndo1_1MicRSBar$data, aes(x=site, y=Abundance)) +
  geom_boxplot(aes(fill=catglab)) +
  geom_point(position=position_dodge(width=0.75), aes(group=catglab)) 


# OTU 3

spistOTUphyloEndo1_3Mic = subset_samples(spistOTUphyloEndo1_3, site=='Micronesia')
spistOTUphyloEndo1_3RS = subset_samples(spistOTUphyloEndo1_3, site=='RedSea')
spistOTUphyloEndo1_3MicRS = merge_phyloseq(spistOTUphyloEndo1_3RS, spistOTUphyloEndo1_3Mic)

spistOTUphyloEndo1_3MicRSBar <- plot_bar(spistOTUphyloEndo1_3MicRS, title='OTU3')

ggplot(spistOTUphyloEndo1_3MicRSBar$data, aes(x=site, y=Abundance, title='OTU3')) +
  geom_boxplot(aes(fill = site)) +
  scale_fill_manual(values=c("#00BA38", "#F8766D")) +
  geom_point(position=position_jitter(width=0.2), aes(group=catglab), size=5) 

# OTU 6

spistOTUphyloEndo1_6Mic = subset_samples(spistOTUphyloEndo1_6, site=='Micronesia')
spistOTUphyloEndo1_6RS = subset_samples(spistOTUphyloEndo1_6, site=='RedSea')
spistOTUphyloEndo1_6MicRS = merge_phyloseq(spistOTUphyloEndo1_6RS, spistOTUphyloEndo1_6Mic)

spistOTUphyloEndo1_6MicRSBar <- plot_bar(spistOTUphyloEndo1_6MicRS, title='OTU6')

ggplot(spistOTUphyloEndo1_6MicRSBar$data, aes(x=site, y=Abundance, title='OTU6')) +
  geom_boxplot(aes(fill = site)) +
  scale_fill_manual(values=c("#00BA38", "#F8766D")) +
  geom_point(position=position_jitter(width=0.2), aes(group=catglab), size=5) 

# OTU 7

spistOTUphyloEndo1_7Mic = subset_samples(spistOTUphyloEndo1_7, site=='Micronesia')
spistOTUphyloEndo1_7RS = subset_samples(spistOTUphyloEndo1_7, site=='RedSea')
spistOTUphyloEndo1_7MicRS = merge_phyloseq(spistOTUphyloEndo1_7RS, spistOTUphyloEndo1_7Mic)

spistOTUphyloEndo1_7MicRSBar <- plot_bar(spistOTUphyloEndo1_7MicRS, title='OTU7')

ggplot(spistOTUphyloEndo1_7MicRSBar$data, aes(x=site, y=Abundance, title='OTU7')) +
  geom_boxplot(aes(fill = site)) +
  scale_fill_manual(values=c("#00BA38", "#F8766D")) +
  geom_point(position=position_jitter(width=0.2), aes(group=catglab), size=5) 


# MED nodes split from the original Otu0001 at 3%

spistMEDEndoRel = transform_sample_counts(spistPhyloEndo, function(x) x / sum(x) )

spistMEDEndoRel1_6247 = subset_taxa(spistPhyloEndo, catglab=="Endozoicomonas_MED000006247")
spistMEDEndoRel1_6549 = subset_taxa(spistPhyloEndo, catglab=="Endozoicomonas_MED000006549")
spistMEDEndoRel1_6809 = subset_taxa(spistPhyloEndo, catglab=="Endozoicomonas_MED000006809")
#spistMEDEndoRel1_4669 = subset_taxa(spistMEDEndoRel, catglab=="Endozoicomonas_MED000004669")
#spistMEDEndoRel1_4794 = subset_taxa(spistMEDEndoRel, catglab=="Endozoicomonas_MED000004794")

spistMEDEndoRel1 = merge_phyloseq(spistMEDEndoRel1_6247, spistMEDEndoRel1_6549, spistMEDEndoRel1_6809)

endMEDbar_1 <- plot_bar(spistMEDEndoRel1, title='MED nodes', fill='catglab')

ggplot(endMEDbar_1$data, aes(x=site, y=Abundance)) +
  geom_boxplot(aes(fill=catglab), alpha=0.8) +
  geom_point(position=position_dodge(width=0.75), aes(group=catglab)) 

endMEDbar_6247 <- plot_bar(spistMEDEndoRel1_6247, title='MED nodes', fill='catglab')

ggplot(endMEDbar_6247$data, aes(x=site, y=Abundance)) +
  geom_boxplot(aes(fill=catglab), alpha=0.8) +
  geom_point(position=position_dodge(width=0.75), aes(group=catglab)) 

# just Micronesia vs Red Sea

spistMEDEndoRel1Mic = subset_samples(spistMEDEndoRel1, site=='Micronesia')
spistMEDEndoRel1RS = subset_samples(spistMEDEndoRel1, site=='RedSea')
spistMEDEndoRel1MicRS = merge_phyloseq(spistMEDEndoRel1Mic, spistMEDEndoRel1RS)

spistMEDEndoRel1MicRSBar <- plot_bar(spistMEDEndoRel1MicRS, title='MED nodes')

ggplot(spistMEDEndoRel1MicRSBar$data, aes(x=site, y=Abundance)) +
  geom_boxplot(aes(fill=catglab), alpha=0.8) +
  geom_point(position=position_dodge(width=0.75), aes(group=catglab)) 

# MED6247

spistMEDEndoRel1_6247Mic = subset_samples(spistMEDEndoRel1_6247, site=='Micronesia')
spistMEDEndoRel1_6247RS = subset_samples(spistMEDEndoRel1_6247, site=='RedSea')
spistMEDEndoRel1_6247MicRS = merge_phyloseq(spistMEDEndoRel1_6247RS, spistMEDEndoRel1_6247Mic)

spistMEDEndoRel1_6247MicRSBar <- plot_bar(spistMEDEndoRel1_6247MicRS, title='MED nodes')

ggplot(spistMEDEndoRel1_6247MicRSBar$data, aes(x=site, y=Abundance, title='MED6247')) +
  geom_boxplot(aes(fill = site)) +
  scale_fill_manual(values=c("#00BA38", "#F8766D")) +
  geom_point(position=position_jitter(width=0.2), aes(group=catglab), size=5) 

# MED6549

spistMEDEndoRel1_6549Mic = subset_samples(spistMEDEndoRel1_6549, site=='Micronesia')
spistMEDEndoRel1_6549RS = subset_samples(spistMEDEndoRel1_6549, site=='RedSea')
spistMEDEndoRel1_6549MicRS = merge_phyloseq(spistMEDEndoRel1_6549RS, spistMEDEndoRel1_6549Mic)

spistMEDEndoRel1_6549MicRSBar <- plot_bar(spistMEDEndoRel1_6549MicRS, title='MED nodes')

ggplot(spistMEDEndoRel1_6549MicRSBar$data, aes(x=site, y=Abundance, title='MED6549')) +
  geom_boxplot(aes(fill = site)) +
  scale_fill_manual(values=c("#00BA38", "#F8766D")) +
  geom_point(position=position_jitter(width=0.2), aes(group=catglab), size=5) 

# MED6809

spistMEDEndoRel1_6809Mic = subset_samples(spistMEDEndoRel1_6809, site=='Micronesia')
spistMEDEndoRel1_6809RS = subset_samples(spistMEDEndoRel1_6809, site=='RedSea')
spistMEDEndoRel1_6809MicRS = merge_phyloseq(spistMEDEndoRel1_6809RS, spistMEDEndoRel1_6809Mic)

spistMEDEndoRel1_6809MicRSBar <- plot_bar(spistMEDEndoRel1_6809MicRS, title='MED nodes')

ggplot(spistMEDEndoRel1_6809MicRSBar$data, aes(x=site, y=Abundance, title='MED6809')) +
  geom_boxplot(aes(fill = site)) +
  scale_fill_manual(values=c("#00BA38", "#F8766D")) +
  geom_point(position=position_jitter(width=0.2), aes(group=catglab), size=5) 

