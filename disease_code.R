#Metagenomic analysis
install.packages("dada2")
library(dada2)

install.packages("DECIPHER")
library(DECIPHER)

install.packages("phangorn")
library(phangorn)

install.packages("ggplot2")
library(ggplot2)

install.packages("Rcpp")
library(Rcpp)

install.packages("Rtools")
library(Rtools)

installed.packages("phyloseq")
library(phyloseq)

installed.packages("Biostrings")
library(Biostrings)

pathD <- "C:/Users/kajal/OneDrive/Desktop/BI Project/disease"
pathH <- "C:/Users/kajal/OneDrive/Desktop/BI Project/healthy"
list.files(pathD)
list.files(pathH)

#sort files to ensure diseased/control reads are in same order
DnFs <- sort(list.files(pathD,pattern = "_R1_001.fastq"))
DnRs <- sort(list.files(pathD,pattern = "_R2_001.fastq"))
HnFs <- sort(list.files(pathH,pattern = "_1.fastq"))
HnRs <- sort(list.files(pathH,pattern = "_2.fastq"))

#filtFs <- sort(list.files(path,pattern = "_Fwd_trim.fastq"))
#filtRs <- sort(list.files(path,pattern = "_Rvs_trim.fastq"))

#Extract sample names,assuming filename have format:SAMPLENAMEXXX.fastq
sample.name_D <- sapply(strsplit(DnFs, "_"),'[',1)
sample.name_H <- sapply(strsplit(HnFs, "_"),'[',1)

#specify the full path 
DnFs <- file.path(pathD, DnFs)
DnRs <- file.path(pathD, DnRs)
HnFs <- file.path(pathH, HnFs)
HnRs <- file.path(pathH, HnRs)

#Plot
plotQualityProfile(DnFs)
plotQualityProfile(DnRs)
plotQualityProfile(HnFs)
plotQualityProfile(HnRs)

# Place filtered files in filtered/ subdirectory
filtFsD <- file.path(pathD,"filtered", paste0(sample.name_D, "_F_filt.fastq"))
filtRsD <- file.path(pathD,"filtered", paste0(sample.name_D, "_R_filt.fastq"))
filtFsH <- file.path(pathH,"filtered", paste0(sample.name_H, "_F_filt.fastq"))
filtRsH <- file.path(pathH,"filtered", paste0(sample.name_H, "_R_filt.fastq"))
names(filtFsD) <- sample.name_D
names(filtRsD) <- sample.name_D
names(filtFsH) <- sample.name_H
names(filtRsH) <- sample.name_H

out <- filterAndTrim(HnFs, filtFsH, HnRs, filtRsH, truncLen=c(260,230),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
head(out)

out <- filterAndTrim(DnFs, filtFsD, DnRs, filtRsD, maxN = 0, maxEE = c(2, 2), 
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = FALSE)  # on windows, set multithread = FALSE
head(out)

#learn error rates
errFD <- learnErrors(filtFsD, multithread=FALSE)
errRD <- learnErrors(filtRsD, multithread=FALSE)
errFH <- learnErrors(filtFsH, multithread=FALSE)
errRH <- learnErrors(filtRsH, multithread=FALSE)
plotErrors(errFD, nominalQ=TRUE)
plotErrors(errRD, nominalQ=TRUE)
plotErrors(errFH, nominalQ=TRUE)
plotErrors(errRH, nominalQ=TRUE)
help("plotErrors")


#sample interfernce
dadaFsD <- dada(filtFsD, err=errFD, multithread=FALSE)
dadaRsD <- dada(filtRsD, err=errRD, multithread=FALSE)
dadaFsH <- dada(filtFsH, err=errFH, multithread=FALSE)
dadaRsH <- dada(filtRsH, err=errRH, multithread=FALSE)
dadaFsD[[1]]
dadaFsH[[1]]

#Merged pairs reads
mergersD <- mergePairs(dadaFsD, filtFsD, dadaRsD, filtRsD, verbose=TRUE)
mergersH <- mergePairs(dadaFsH, filtFsH, dadaRsH, filtRsH, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergersD[[1]])
head(mergersH[[1]])

#Construct sequence table
seqtabD <- makeSequenceTable(mergersD)
seqtabH <- makeSequenceTable(mergersH)
dim(seqtabD)
dim(seqtabH)

help("makeSequenceTable")

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtabD)))
table(nchar(getSequences(seqtabH)))

#Remove Chimeras
seqtabD.nochim <- removeBimeraDenovo(seqtabD, method="consensus", multithread=TRUE, verbose=TRUE)
seqtabH.nochim <- removeBimeraDenovo(seqtabH, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtabD.nochim)
dim(seqtabH.nochim)

sum(seqtabD.nochim)/sum(seqtabD)
sum(seqtabH.nochim)/sum(seqtabH)

#Assign Taxonomy
taxaD <- assignTaxonomy(seqtabD.nochim, "C:/Users/kajal/OneDrive/Desktop/BI Project/Result_disease/silva_nr_v132_train_set.fa.gz", multithread=FALSE)
taxaH <- assignTaxonomy(seqtabH.nochim, "C:/Users/kajal/OneDrive/Desktop/BI Project/Result_disease/silva_nr_v132_train_set.fa.gz", multithread=FALSE)
taxaD.print <- taxaD # Removing sequence row names for display only
taxaH.print <- taxaH # Removing sequence row names for display only
rownames(taxaD.print) <- NULL
head(taxaD.print)
rownames(taxaH.print) <- NULL
head(taxaH.print)
help ("assignTaxonomy")


#align sequence
sequenceD <- getSequences(seqtabD.nochim)
sequenceH <- getSequences(seqtabH.nochim)
names(sequenceD) <- sequenceD
names(sequenceH) <- sequenceH
#run sequence alignment using DECIPHER
alignmentD <- AlignSeqs(DNAStringSet(sequenceD), anchor = NA)
alignmentH <- AlignSeqs(DNAStringSet(sequenceH), anchor = NA)

#Construct Phylogenetic tree
#change sequence alignment output into a phyDat structure
phang.alignD <- phyDat(as(alignmentD, "matrix"), type = "DNA")
phang.alignH <- phyDat(as(alignmentH, "matrix"), type = "DNA")
#create distance matrix
dmD <- dist.ml(phang.alignD)
dmH <- dist.ml(phang.alignH)
#perform neighbour joining 
treeNJD <- NJ(dmD)
treeNJH <- NJ(dmH)
#internal maximum likelihood
fitD = pml(treeNJD, data = phang.alignD)
fitH = pml(treeNJH, data = phang.alignH)
fitDGTR <- update(fitD, k= 4, inv = 0.2)
fitHGTR <- update(fitH, k= 4, inv = 0.2)
fitDGTR <- optim.pml(fitDGTR, model = "GTR", optInv = TRUE, optGamma = TRUE,
                     rearrangement = "stochastic", control = pml.control(trace = 0))
fitHGTR <- optim.pml(fitHGTR, model = "GTR", optInv = TRUE, optGamma = TRUE,
                     rearrangement = "stochastic", control = pml.control(trace = 0))
#Import into phyloseq
samplesD.out <- rownames(seqtabD.nochim)
samplesH.out <- rownames(seqtabH.nochim)

samples.out <- c(samplesH.out, samplesD.out)

subjectD <- sapply(strsplit(samplesD.out, "D"), `[`, 2)
subjectD <- substr(subjectD,2,999)

subjectH <- sapply(strsplit(samplesH.out, "9"), `[`, 2)

healthdf <- data.frame(Subject=subjectH, Status='Healthy')
rownames(healthdf) <- samplesH.out

diseasedf <- data.frame(Subject=subjectD, Status="Diseased")
rownames(diseasedf) <- samplesD.out
  
#Construct phyloseq object directly from dada2 
healthyps <- phyloseq(otu_table(seqtabH.nochim, taxa_are_rows=FALSE), 
               sample_data(healthdf), 
               tax_table(taxaH))
healthyps <- prune_samples(sample_names(healthyps) != "Mock", healthyps) # Remove mock sample

diseasedps <- phyloseq(otu_table(seqtabD.nochim, taxa_are_rows=FALSE), 
                      sample_data(diseasedf), 
                      tax_table(taxaD))
diseasedps <- prune_samples(sample_names(diseasedps) != "Mock", diseasedps)# Remove mock sample

ps = merge_phyloseq(healthyps, diseasedps)
  
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps
  
#visualize alpha_diversity
plot_richness(ps, x="Subject", measures=c("Shannon", "Simpson"), color="Status")

#ordinate
# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="PCoA", distance="bray")
plot_ordination(ps.prop, ord.nmds.bray, color="Subject", title="Bray NMDS")


top70 <- names(sort(taxa_sums(ps), decreasing=TRUE))[30:100]
ps.top70 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top70 <- prune_taxa(top70, ps.top70)
plot_bar(ps.top70, x="Subject", fill="Family") + facet_wrap(~Status, scales="free_x")

help(plotQualityProfile)
