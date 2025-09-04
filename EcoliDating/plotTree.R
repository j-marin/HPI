library(ape)
library(treeio)
library(ggtree)
library(ggplot2)
library(TreeSim)

## plot TDR tree and color clades by phylogroup
tree1 <- read.tree("../data/TDR_tree.nwk")
tree <- drop.tip(tree1, c("001-014", "003-014", "CP054715-1", "CP028169-1", "CP023468-1", "ROAR354")) # E. coli s.s.

a <- read.table("../data/typing.txt", head=TRUE)
head(a)
max(getx(tree1)) 

## E coli sensus stricto 
tree2 <- tree
tip.n <- paste(a$Phylogroup[match(tree$tip.label, a$Strain)], tree2$tip.label, sep="_")
groupInfo <- split(tip.n, gsub("\\_..*", "", tip.n))
tree$tip.label <- tip.n
tree <- groupOTU (tree, groupInfo)

cols <- c("#000000", "#454fd4", "#18923f", "#ec2d18", "#579aeb", "#eeda12", "#7e1892", "#ec911f", "#e719df", "#9399a1", "black", "black", "black", "black")

ggtree(tree, mrsd='2016-01-01', aes(color= group)) + scale_color_manual(values = cols)  + theme_tree2() +
                theme(panel.grid.major   = element_line(color="black", size=.2),
                      panel.grid.minor   = element_line(color="grey", size=.2),
                      panel.grid.major.y = element_blank(),
                      panel.grid.minor.y = element_blank()) 

## E coli +  outgroup 
tree2 <- tree1
tip.n <- paste(a$Phylogroup[match(tree1$tip.label, a$Strain)], tree2$tip.label, sep="_")
groupInfo <- split(tip.n, gsub("\\_..*", "", tip.n))

tree1$tip.label <- tip.n
tree1 <- groupOTU (tree1, groupInfo)
attr(tree1, "group") <- factor(attr(tree1, "group"), levels = c("0","A","B1","B2","C","D","E","F","G","H","undescribed","-"))

cols <- c("#000000", "#454fd4", "#18923f", "#ec2d18", "#579aeb", "#eeda12", "#7e1892", "#ec911f", "#e719df", "#9399a1", "black", "black", "black", "black")

ggtree(tree1, mrsd='2016-01-01', aes(color= group)) + scale_color_manual(values = cols)  + theme_tree2() +
                theme(panel.grid.major   = element_line(color="black", size=.2),
                      panel.grid.minor   = element_line(color="grey", size=.2),
                      panel.grid.major.y = element_blank(),
                      panel.grid.minor.y = element_blank()) 

