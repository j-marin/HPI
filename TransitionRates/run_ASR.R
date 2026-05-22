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
print(tree)

# check names
setdiff(tree$tip.label, a$Strains)
tree$edge.length[which(tree$edge.length < 0)] <- 0.0000001
tree$node.label <- NULL

# prep the states file
states <- a$HPI[match(tree$tip.label, a$Strains)]
names(states) <- tree$tip.label

#### Split tree in subtrees
times <- getx(tree)
hpi_intro <- 4196.66 # HPI introduction
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
clades.s <- clades.s [-unique(r)] # nodes
	
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

# divide branch length by 10 reach convergence and /10 obtained rates
#for (i in 1:nclade) {
#	extracted_trees2[[i]]$edge.length = extracted_trees2[[i]]$edge.length/10
#	}
	
### Run ARD model (with constrained priors)
# parameters
alpha1 <- 1.40		#q10
beta1 <- 10

alpha2 <- 0.615		#q01
beta2 <- 10

# parameters
nsamp <- 1000

# outputs
clade_res <- matrix(ncol=5, nrow=length(extracted_trees2))
q10_list <- list()
q01_list <- list()
Tab_count <- data.frame(times=NA, count=NA, nb_br=0, nsamp_id=NA, cl_id=NA)

for (m in 1:length(extracted_trees2)) {
	# MCMC approach to sample character histories from their posterior probability distribution
	mtrees <- make.simmap(extracted_trees2[[m]], extracted_states2[[m]], model="ARD", nsim=nsamp, pi=c(1,0), Q="mcmc", prior=list(alpha=c(alpha1, alpha2), beta=c(beta1, beta2), use.empirical=FALSE), vQ=1e-02, samplefreq=100) # pi to contrain the root at state "no" / prior and vQ to constrain Q matrix (prior and variance)
		
	q10_mtrees <- sapply(1:nsamp, function(i) mtrees[[i]]$Q[2,1])
	q01_mtrees <- sapply(1:nsamp, function(i) mtrees[[i]]$Q[1,2])
	
	q10_list <- c(q10_list, list(q10_mtrees))
	q01_list <- c(q01_list, list(q01_mtrees))

	ESS <- effectiveSize(as.mcmc(cbind(q10_mtrees, q01_mtrees)))  
	print(ESS>200)
	print(ESS)
		
	## get frequency through time
	mtrees2 <- mtrees
	for (k in 1:nsamp) {	
		mtrees2[[k]] <- drop.tip(mtrees[[k]], "outgroup")
		tab_tree <- cbind(mtrees2[[k]]$edge, mtrees2[[k]]$edge.length)
		tab_tree <- cbind(getx(mtrees2[[k]])[match(tab_tree[,1], names(getx(mtrees2[[k]])))], tab_tree)
		colnames(tab_tree) <- c("node time", "ancestor", "descendant", "branch length")
	
		# counter - nb of branchs / time
		cp_b <- sort(getx(mtrees2[[k]]), decreasing=TRUE)
		yes_no <- mtrees2[[k]]$maps # state duration 
				
		# counter - nb of yes / time
		branch_yn <- list()
		for (z in 1:length(yes_no)) {
			timez <- tab_tree[z,1] - cumsum(yes_no[[z]])
			branch_yn[[z]] <- c(tab_tree[z,1], timez[-length(timez)])
			names(branch_yn[[z]]) <- names(yes_no[[z]])
		}
		
		tab_count <- data.frame(times=unique(sort(unlist(branch_yn), decreasing=TRUE)), count=0, nb_br=1, nsamp_id=k, cl_id=m)
		
		# branch counter
		br <- match((mtrees2[[k]]$Nnode+2):(mtrees2[[k]]$Nnode+1+mtrees2[[k]]$Nnode), tab_tree[,3]) # select branches
		br <- br[!is.na(br)]
		for (j in br) {
			cp_y <- names(branch_yn[[j]])
			cp_y <- ifelse(cp_y == "no", 0, 1)
			time_yn <- branch_yn[[j]]				
			tab_count$count[match(time_yn, tab_count$times)] <- tab_count$count[match(time_yn, tab_count$times)] + cp_y
		}
		
		# end of branches counter
		for (j in br) {
			cp_y <- names(branch_yn[[j]])
			cp_y <- ifelse(cp_y == "no", 0, 1)
			cp_y <- cp_y[length(cp_y)]
			time_yn <- branch_yn[[j]][length(branch_yn[[j]])]	
			br_end <- tab_tree[j,1] - tab_tree[j,4]
						
			tab_count$count[round(tab_count$times, 4) > round(br_end, 4) & round(tab_count$times, 4) < round(time_yn, 4)] <- tab_count$count[round(tab_count$times, 4) > round(br_end, 4) & round(tab_count$times, 4) < round(time_yn, 4)] + cp_y
		}
			
		# tip counter
		for (p in 1:(mtrees2[[k]]$Nnode+1)) {
			b <- unlist(branch_yn[which(tab_tree[,3]== p)])
			if(!is.null(b)){
				b <- b[length(b)]
				if (names(b) == "yes") {
					tab_count$count[which(tab_count$times <= b)] <- tab_count$count[which(tab_count$times <= b)] + 1
				}
			}
		}
		
		# nb of branchs
		for (n in 1:length(cp_b)) {
			tab_count$nb_br[which(tab_count$times <= cp_b[n])] <- tab_count$nb_br[which(tab_count$times <= cp_b[n])] + 1
		}
		
		# add stem branch
		yes_no2 <- mtrees[[k]]$maps[which.max(mtrees[[k]]$edge.length[-match(hpi_intro, mtrees[[k]]$edge.length)])]
		timez2 <- hpi_intro - cumsum(yes_no2[[1]])
		branch_yn2 <- c(hpi_intro, timez2[-length(timez2)])
		names(branch_yn2) <- names(yes_no2[[1]])
		
		tab_count2 <- data.frame(times=unique(sort(unlist(branch_yn2), decreasing=TRUE)), count=0, nb_br=1, nsamp_id=k, cl_id=m)
		tab_count2$count[which(names(branch_yn2)== "yes")] <- 1	
		
		# all in one table
		Tab_count <- rbind(Tab_count, tab_count, tab_count2)	
	}	
	
	# mean change rates
	yn <- mean(q10_mtrees)
	ny <- mean(q01_mtrees)

	# node age and tip nb
	nbN <- mtrees2[[1]]$Nnode
	nbT <- length(mtrees2[[1]]$tip.label)
	node_age <- getx(mtrees2[[1]])

	# clade summary results
	clade_res[m,] <- c(m, ny, yn, max(node_age), nbT)
	
	print(m)
	}
	
# prep outputs
Tab_count <- Tab_count[-1,]
Tab_count <- Tab_count[order(Tab_count$times, decreasing=TRUE),]	
write.table(Tab_count, "Tab_count.txt", quote= FALSE, row.names=FALSE, sep="	")

clade_res <- as.data.frame(clade_res)
colnames(clade_res) <- c("cladeNb", "q01", "q10", "mrca", "Ntips")
write.table(clade_res, "ace_summary_table_constrainedQ.txt", quote= FALSE, row.names=FALSE, sep="	")

write.table(q10_list, "q10_list.txt", quote= FALSE, row.names=FALSE, col.names=FALSE, sep="	")
write.table(q01_list, "q01_list.txt", quote= FALSE, row.names=FALSE, col.names=FALSE, sep="	")



##################################################################################################
# Open results
clade_res <- read.csv("ace_summary_table_constrainedQ.txt", header=TRUE, sep="	")

q10_list  <- read.table("q10_list.txt", header=FALSE)
q01_list  <- read.table("q01_list.txt", header=FALSE)

# Compute transition rates
q01 <- as.vector(as.matrix(q01_list))
q10 <- as.vector(as.matrix(q10_list))

n <- length(q01)
s01 <- sd(q01)
s10 <- sd(q10)

lowerI01 <- mean(q01) - (qt(0.975, df=n-1)*s01/sqrt(n))
upperI01 <- mean(q01) + (qt(0.975, df=n-1)*s01/sqrt(n))

lowerI10 <- mean(q10) - (qt(0.975, df=n-1)*s10/sqrt(n))
upperI10 <- mean(q10) + (qt(0.975, df=n-1)*s10/sqrt(n))

formatC(c(mean(q01), lowerI01, upperI01), format="e", digits=3)
formatC(c(mean(q10), lowerI10, upperI10), format="e", digits=3)


## Plot frequency through time
library(ggplot2)

# --> use make_freq_tab.r to get mean frequency through time
freq_tab <- read.table("hpi_frequency_time.txt", skip=1)

# get mrca of each clade
mrca_list <- integer()
for (i in 1:length(extracted_trees2)) {
	mrca_list[i] <- max(getx(drop.tip(extracted_trees2[[i]], "outgroup")))
	}	

# get collection frequencies
d_coll <- c(1980, 2000, 2001, 2002, 2010, 2016, 2016)
p_coll <- c(18/53, 26/50, 16/27, 13/27, 151/246, 11/20, 0.5866667)

d_coll <- c(1980, 2001, 2010, 2016, 2016)
p_coll <- c(18/53, c(26+16+13-6)/c(50+27+27-6), (151-9)/(246-9), 11/20, 0.5866667)

data_coll <- data.frame(date=2016-d_coll, prop=p_coll)
ratio_eq <- mean(q01) / (mean(q10) +  mean(q01)) # ratio attendu à l'équilibre

pdf(file="HPI_prop_time.pdf", width=6, height=4)

p1 <- ggplot(freq_tab, aes(V1, V4)) + 
	  geom_vline(aes(xintercept=mrca_list), data=data.frame(mrca_list), linetype="dashed", color="gray60", lwd=0.2) +
	  geom_hline(aes(yintercept=ratio_eq), linetype="dashed", color="black", lwd=0.5) +
      geom_point(aes(V1, V5), col="gray80", size=0.2) +
      geom_point(aes(V1, V6), col="gray80", size=0.2) +
      geom_point(size=0.2) +
      geom_point(data=data_coll, aes(date, prop), col=c(rep("blue", 4), "red")) +
      scale_x_sqrt(breaks=c(c(0,100,500), seq(1000,4000,1000))) +
	  xlab("Time (years ago since 2016)") + 
	  ylab("HPI proportion") + 
	  geom_segment(aes(x=0, y=0.1, xend=70, yend=0.1), size=0.2, arrow=arrow(ends="both", length=unit(0.08, "inches"))) +
	  annotate("text", x=18, y=0.14, label="Selection") +
	  geom_segment(aes(x=3300, y=0.5, xend=4200, yend=0.5), size=0.2, arrow=arrow(ends="both", length=unit(0.08, "inches"))) +
	  annotate("text", x=3600, y=0.54, label="Rapid stabilisation") +
	  theme_classic()
p1

dev.off()

############### ############ ###########
### Extract phylogroup results
b <- read.table("../data/typing.txt", header=TRUE)

# get clade phylogroup
phylogp <- sapply(1:length(extracted_states), function(i)(unique(b$Phylogroup[b$Strain %in% names(extracted_states[[i]])])))

f2 <- c("all", "A", "B1", "B2", "D", "F")
res <- matrix(ncol=5)
for (i in 1:length(f2)) {
	# select phylogroup results
	q10_list_p <- q10_list[,grep(f2[i], phylogp)]
	q01_list_p <- q01_list[,grep(f2[i], phylogp)]

	if (f2[i] == "all") {
		q10_list_p <- q10_list
		q01_list_p <- q01_list
		}

	print(length(extracted_states[grep(f2[i], phylogp)]))
	print(length(unlist(extracted_states[grep(f2[i], phylogp)])))
	
	# visualize freq
	#plot(freq_tab2$V1, freq_tab2$V4)

	# changing rates
	q01 <- mean(as.vector(as.matrix(q01_list_p)))
	q10 <- mean(as.vector(as.matrix(q10_list_p)))

	n <- length(as.vector(as.matrix(q01_list_p)))
	s01 <- sd(as.vector(as.matrix(q01_list_p)))
	s10 <- sd(as.vector(as.matrix(q10_list_p)))

	lowerI01 <- q01 - (qt(0.975, df=n-1)*s01/sqrt(n))
	upperI01 <- q01 + (qt(0.975, df=n-1)*s01/sqrt(n))

	lowerI10 <- q10 - (qt(0.975, df=n-1)*s10/sqrt(n))
	upperI10 <- q10 + (qt(0.975, df=n-1)*s10/sqrt(n))

	res01 <- formatC(c(mean(q01), lowerI01, upperI01), format="e", digits=3)
	res10 <- formatC(c(mean(q10), lowerI10, upperI10), format="e", digits=3)

	
	# get formated results
	print(paste(res01[1], " [", res01[2], "-", res01[3], "]", sep="")) # q01
	print(paste(res10[1], " [", res10[2], "-", res10[3], "]", sep="")) # q10

	# get res as data.frame 	
	value <- c(mean(q01), mean(q10)) # q01 q10
	lower <- c(lowerI01, lowerI10)
	upper <- c(upperI01, upperI10)
	parameter <- c("q01", "q10")
	
	res_f <- cbind(f2[i], value, lower, upper, parameter)
	
	res <- rbind(res, res_f)
	}

res <- res[-1,]
colnames(res) <- c("group", "value", "lower", "upper", "parameter")
res_df <- as.data.frame(res)
res_df[,2:4] <- apply(res_df[,2:4], 2, as.numeric)
apply(res_df[,2:4], 2, formatC, format="e", digits=3)

# plot results
library(ggplot2)

res_df$group <- factor(res_df$group, levels=c("all", "A", "B1", "B2", "D", "F"))

p2 <- ggplot(res_df, aes(x=group, y=value)) +
     geom_errorbar(aes(ymin=lower, ymax=upper, color=parameter), position = position_dodge(0.3)) + 
     geom_point(aes(color=parameter), position=position_dodge(0.3)) + 
     scale_y_continuous(trans="log10") + 
     scale_colour_manual(values=c("gray70", "gray30")) +
     theme_classic() + 
     ylab("Estimate") +
     xlab("Phylogroup")


pdf(file="changing_rates_phylogroup.pdf", width=6, height=4)
p2
dev.off()


##### Combine plots
library(gridExtra)

pdf(file="Transition_plots.pdf", width=12, height=4)
grid.arrange(p1, p2, nrow=1)
dev.off()





