library(phytools)
library(TreeSim)
library(zoo)
library(coda)
library(phylo_utils)

rm(list=ls())

# get data
a <- read.csv("../data/recapHPI.csv", header=TRUE, sep=";")
table(a$HPI, a$Collection)

tree <- read.tree("../data/TDR_tree.nwk")
tree <- drop.tip(tree, c("001-014", "003-014", "CP054715-1", "CP028169-1", "CP023468-1", "ROAR354")) # not E. coli s.s.

# drop outgroup and murray seq
strains <- tree$tip.label
murray <- strains[startsWith(strains, "M")]
tree <- drop.tip(tree, murray) # 420

# Exclude 15 strains from the `K. pneumoniae clade
hpi_tree <- read.tree("../data/HPItree_rooted_root_to_tips.nwk") #1284 tips
x <- hpi_tree$tip.label
x[startsWith(x, "C0")] <- sub(".", "", x[startsWith(x, "C0")])
x[grep("PAR", x)] <- gsub("_", "-", x[grep("PAR", x)])
x[grep("_DOM", x)] <- gsub("_", "-", x[grep("_DOM", x)])
setdiff(x, a$Strains)

hpi_tree$tip.label <- x
plot(ladderize(hpi_tree))
nodelabels()
cladeK <- extract.clade(hpi_tree, 1615)
table(a$Species[match(cladeK$tip.label, a$Strains)])
coliK <- cladeK$tip.label[which(a$Species[match(cladeK$tip.label, a$Strains)] == "Escherichia coli")]

tree <- drop.tip(tree, coliK) # 405

# check names
setdiff(tree$tip.label, a$Strains)
tree$edge.length[which(tree$edge.length < 0)] <- 0.0000001
tree$node.label <- NULL

# prep the states file
states <- a$HPI[match(tree$tip.label, a$Strains)]
names(states) <- tree$tip.label

## Split tree in subtrees
times <- getx(tree)
hpi_intro <- 4185.86 # HPI introduction
clades.s <- as.numeric(names(times)[which(times<hpi_intro)])

# remove embedded clades
cl <- integer(0)
for (i in 1:length(clades.s)) {
	ta <- extract.clade(tree, clades.s[i])
	cl <- c(cl, list(ta$tip.label))
	}

r <- integer(0)
for(j in 1:length(cl))
for (i in 1:length(cl)) {
	if (length(intersect(cl[[j]], cl[[i]])) > 0 ) {
		if (length(cl[[j]]) > length(cl[[i]])) { r <- c(r, i) } 
		if (length(cl[[j]]) < length(cl[[i]])) { r <- c(r, j) } 
		}}

cl <- cl [-unique(r)]
id_in_clades <- unlist(cl)
stopifnot(all(!duplicated(id_in_clades)))
clades.s <- clades.s [-unique(r)] # nodes nb
	
# singletons (not in clades)
singletons <- tree$tip.label[!tree$tip.label %in% id_in_clades]
singleton_states <- states[match(singletons, names(states))]
table(singleton_states)
	
# prep data
nclade <- length(cl)
clade_sizes <- numeric()
extracted_trees <- lapply(1:nclade, function(i) extract.clade(tree, clades.s[i]))

clade_sizes <- c(clade_sizes, sapply(1:nclade, function(i) length(extracted_trees[[i]]$tip.label)))
extracted_states <- lapply(1:nclade, function(i) states[match(extracted_trees[[i]]$tip.label, names(states))])

# add outgroup rooted at HPI introduction time to each clades
extracted_trees2 <- list()
for (i in 1:length(extracted_trees)) {
	edgeL <- hpi_intro - max(getx(extracted_trees[[i]]))
	t1 <- add_root_edge2(extracted_trees[[i]], edgeL)
	tip <- list(edge=matrix(c(2,1),1,2), tip.label="outgroup", edge.length=hpi_intro, Nnode=1)
	class(tip)<-"phylo"
	extracted_trees2[[i]] <- bind.tree(t1, tip, where="root", position=edgeL)	
}

extracted_states2 <- list()
for (i in 1:length(extracted_states)) {
	extracted_states2[[i]] <- c(extracted_states[[i]], "outgroup"="no")	
}

# change status outgoup to 'yes' for clades with 'no' state only
select_no <- which(lapply(1:length(extracted_trees2), function (i) length(unique(extracted_states2[[i]]))) == 1) # no state only
for (i in select_no) {
	extracted_states2[[i]][which(names(extracted_states2[[i]]) == "outgroup")] <- "yes"	
}

# remove clades with only one state
select1 <- which(lapply(1:length(extracted_trees2), function (i) length(unique(extracted_states2[[i]]))) == 2)

# get final trees and states
extracted_states2 <- extracted_states2[select1]
extracted_trees2 <- extracted_trees2[select1]

### Run ARD model
# parameters
niter <- 1

# outputs
clade_res <- matrix(ncol=5, nrow=length(extracted_trees2))

for (m in 1:length(extracted_trees2)) {
	## parameter estimation
	mtrees <- make.simmap(extracted_trees2[[m]], extracted_states2[[m]], model="ARD", nsim=niter, pi=c(1,0), Q="empirical") # pi to contrain the root at state "no" / Q="empirical"

	# get Q matrix and change rates
	yn <- mtrees$Q[2,1]
	ny <- mtrees$Q[1,2]

	# get number of nodes and tips
	nbN <- mtrees$Nnode
	nbT <- length(mtrees$tip.label)

	# get clade age
	node_age <- getx(mtrees)

	# clade summary results
	clade_res[m,] <- c(m, ny, yn, max(node_age), nbT)
	
	print(m)
	}
	
clade_res <- as.data.frame(clade_res)
colnames(clade_res) <- c("cladeNb", "q01", "q10", "mrca", "Ntips")


# get parameters to constrain the Q matrix
mean(as.numeric(clade_res$q01)) 
mean(as.numeric(clade_res$q10)) 



