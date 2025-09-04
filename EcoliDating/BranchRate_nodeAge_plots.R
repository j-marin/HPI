library(ape)
library(phytools)
library(TreeSim)
library(bnpsd)
library(ggplot2)
library(gridExtra)

rm(list=ls())

# open data
a <- read.table("../data/typing.txt", header=TRUE)
tree <- read.tree("../data/core_ecoli_rooted.nwk")
tree.lsd <- read.nexus("../data/lsd2_dating_core.date.nexus")
tree.d <- read.tree("../data/TDR_tree_v2.nwk")

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
setdiff(tree.d$tip.label, a$Strain) # all good
setdiff(tree.bd$tip.label, a$Strain) # all good

cbind(tree.lsd$tip.label, tree.d$tip.label)

## TDR tree
# reorder phylogenetic tree
tree.r <- tree_reorder(tree, tree.d$tip.label)
write.tree(tree.r, "../data/core_ecoli_ordered.nwk")
tree.r <- read.tree("../data/core_ecoli_ordered.nwk")

stopifnot(tree.r$tip.label == tree.d$tip.label)
stopifnot(tree.r$edge[,1] == tree.d$edge[,1])
stopifnot(tree.r$edge[,2] == tree.d$edge[,2])

anc.nd <- tree.d$edge[,1]
b1 <- tree.r$edge.length # phylogenetic tree
b2 <- tree.d$edge.length # timetree
rate <- (b1/b2)
names(rate) <- anc.nd

ages <- getx(tree.d)

ages2 <- numeric()
for (i in 1:length(anc.nd)) {
	ages2[i] <- ages[which(names(ages) == anc.nd[i])]
}

plot(ages2, rate, log="xy", pch=19)

df1 <- data.frame(rate=rate, age=ages2)

p1 <- ggplot(df1, aes(age, rate)) +
     geom_point() +
     scale_y_continuous(trans = "log10", limits = c(1e-015,1)) +
     scale_x_continuous(trans = "log10", limits = c(0.01,100000000)) + 
     theme_minimal() +
     ggtitle ("TDR tree") +
     xlab ("node age (Ma)") + 
     ylab ("branch rate (subs/Ma)")
                
p1 <- p1 + geom_smooth(method = "loess")
  
      
## LSD tree
# compute branch rate vs time
tree.r <- tree_reorder(tree, tree.lsd$tip.label)
write.tree(tree.r, "../data/core_ecoli_ordered_lsd.nwk")
tree.r <- read.tree("../data/core_ecoli_ordered_lsd.nwk")

stopifnot(tree.r$tip.label == tree.lsd$tip.label)
stopifnot(tree.r$edge[,1] == tree.lsd$edge[,1])
stopifnot(tree.r$edge[,2] == tree.lsd$edge[,2])

anc.nd <- tree.lsd$edge[,1]
b1 <- tree.r$edge.length # phylogenetic tree
b2 <- tree.lsd$edge.length # timetree
rate <- (b1/b2)
names(rate) <- anc.nd

ages <- getx(tree.lsd)

ages2 <- numeric()
for (i in 1:length(anc.nd)) {
	ages2[i] <- ages[which(names(ages) == anc.nd[i])]
}

ages2 <- ages2[-c(which(rate == "Inf"))]
rate <- rate[-c(which(rate == "Inf"))]

plot(ages2, rate, log="xy", pch=19)

df2 <- data.frame(rate=rate, age=ages2)

p2 <- ggplot(df2, aes(age, rate)) +
      geom_point() + 
      scale_y_continuous(trans = "log10", limits = c(1e-015,1)) +
      scale_x_continuous(trans = "log10") + 
      theme_minimal() +
      ggtitle ("LSD tree") +
      xlab ("node age (Ma)") + 
      ylab ("branch rate (subs/Ma)")    

p2 <- p2 + geom_smooth(method = "loess")



pdf(file="RateTimePlots_nodeAge.pdf", width=12, height=4)

grid.arrange(p1, p2, nrow=1)

dev.off()

