library(Rlsd2)
library(ape)
library(ips)
library(phytools)
library(seqinr)
library(TreeSim)

rm(list=ls())
source("lsd_function.r")

## open data
# open data
a <- read.table("../data/typing.txt", header=TRUE)
tree <- read.tree("../data/core_ecoli_rooted.nwk")
aln <- read.fasta("../data/core_gene_alignment_clean.fasta")

seqL <- length(aln[[1]])

# check names
strains <- tree$tip.label
strains[startsWith(strains, "C0")] <- gsub("C0", "0", strains[startsWith(strains, "C0")])
strains[startsWith(strains, "X")] <- sub(".", "", strains[startsWith(strains, "X")])
strains <- gsub("[.]", "-", strains)
strains <- gsub("-", "-", strains)
strains[match(setdiff(strains, a$Strain), strains)] <- gsub("_", "-", strains[match(setdiff(strains, a$Strain), strains)])
strains <- gsub("-PETANC-allcontigs", "", strains)
tree$tip.label <- strains

setdiff(strains, a$Strain) # all good

# get dates
dates <- a$Date[match(tree$tip.label, a$Strain)]
names(dates) <- tree$tip.label

# replace 0 branch length
tree$edge.length[tree$edge.length==0] <- 5e-09

# run lsd2 with root dating
res <- lsd2(inputTree=tree, inputDate=dates, outFile = "lsd2_dating_core", seqLen= seqL, rootDate = -1.02e+08, nullblen = 4e-09, constraint=TRUE, variance=2, confidenceInterval = 1000)


