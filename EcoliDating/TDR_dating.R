library(ape)
library(RRphylo)
library(phytools)
library(TreeSim)
library(adephylo)
library(treeio)
library(dplyr)

# rate equation parameters to get
mu <- 1.778415e-05 # alpha
lambda <- 5.209e-03 # beta
k <- 1.821709e-09

# lower CI
mu <- 1.463065e-05 # alpha
lambda <- 6.784954e-03 # beta
k <- 1.821709e-09

# upper CI
mu <- 2.093766e-05 # alpha
lambda <- 3.633907e-03 # beta
k <- 1.821709e-09

## open date file and trees
a <- read.table("../data/typing.txt", header=TRUE)
tree <- read.tree("../data/core_ecoli_rooted.nwk")
tree.lsd <- read.nexus("lsd2_dating_core.date.nexus")

# check names
strains <- tree.lsd$tip.label
strains[startsWith(strains, "C0")] <- gsub("C0", "0", strains[startsWith(strains, "C0")])
strains[startsWith(strains, "X")] <- sub(".", "", strains[startsWith(strains, "X")])
strains <- gsub("[.]", "-", strains)
strains <- gsub("-", "-", strains)
strains[match(setdiff(strains, a$Strain), strains)] <- gsub("_", "-", strains[match(setdiff(strains, a$Strain), strains)])
strains <- gsub("-PETANC-allcontigs", "", strains)
tree.lsd$tip.label <- strains

setdiff(strains, a$Strain)

# rescale the lsd tree
tree.lsd$edge.length <- (tree.lsd$edge.length * max(getx(tree))) / max(getx(tree.lsd)) # same scale than fasttree
tree.lsd$edge.length[tree.lsd$edge.length==0] <- 1e-08

# apply TDR method
Ftosolve <- function(td, d) {
	return(
	(1/lambda) * (-mu *exp(-lambda *td) + k *lambda * td + mu) -d
	)
}

# to node age
tA <- tree.lsd
node.age <- integer(0)
for (i in (length(tA$tip.label) + (1:Nnode(tA)))) {
	
	# calculate d: node-to-tip distances
	tx <- extract.clade(tA, i)
	
	if (Nnode(tx) > 1) {
		d.all <- distNodes (tx, node = length(tx$tip.label) + 1)
		dx <- d.all[match(tx$tip.label, rownames(d.all)), 2] }
		
	if (Nnode(tx) == 1) {
		dx <- distRoot(tx) }

	# estimate node age
	tdsol <- integer(0)
	for (j in 1:length(dx)) {
	tdsol <- c(tdsol, uniroot(Ftosolve, upper= 1E9, c(1E-13, 1E7), d=dx[j])$root)
	}

	node.age <- rbind(node.age, c(i, mean(tdsol), sd(tdsol), mean(dx)/mean(tdsol)))
print(i)
}
node.age <- data.frame(node.age)
colnames(node.age) <- c("node_nb", "node_age", "age_sd", "rate")

# adjust edges
tA.d <- tA
tab.edge <- cbind(tA$edge, tA$edge.length)
new.edge <- numeric()

for (i in 1:dim(tab.edge)[1]) {
	
	n1 <- tab.edge[i,1]
	n2 <- tab.edge[i,2]
	
	age1 <- node.age[node.age[,1] == n1, 2]
	age2 <- node.age[node.age[,1] == n2, 2]
		if (length(age2) == 0) {age2 <- 0}
		
	edge.x <- age1-age2
	
	tA.d$edge.length [i] <- edge.x
	new.edge[i] <- edge.x
	}


write.tree(tA.d, "../data/TDR_tree_final.nwk")


