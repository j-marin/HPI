library(ape)
library(adephylo)
library(phytools)
library(seqinr)
library(stringr)
library(Biostrings)
library(rlist)
library(BactDating)
library(nlstools)

rm(list=ls())

# open data
a <- read.table("../data/typing.txt", header=TRUE)
tree <- read.tree("../data/core_ecoli_rooted.nwk")

# check names
strains <- tree$tip.label
strains[startsWith(strains, "C0")] <- gsub("C0", "0", strains[startsWith(strains, "C0")])
strains[startsWith(strains, "X")] <- sub(".", "", strains[startsWith(strains, "X")])
strains <- gsub("[.]", "-", strains)
strains <- gsub("-", "-", strains)
strains[match(setdiff(strains, a$Strain), strains)] <- gsub("_", "-", strains[match(setdiff(strains, a$Strain), strains)])
strains <- gsub("-PETANC-allcontigs", "", strains)
tree$tip.label <- strains

setdiff(strains, a$Strain)

# plot the tree
plot(tree)
nodelabels()

### Molecular clock signal
# check for molecular clock signal
tA <- extract.clade(tree, 457)
plot(tA) # 444 tips

root.tip <- diag(vcv.phylo(tA))
year <- as.numeric(a$Date[match(names(diag(vcv.phylo(tA))), a$Strain)])
cor.test(year, root.tip)

# select STs with positive relationship (root-to-tip / sampling date)
infos <- a[match(tA$tip.label, a$Strain),]
infos <- infos[!is.na(infos$Strain),]
sts <- names(which(table(infos$ST) > 6) [-1])

pvals <- rep(NA, length(sts))
rate <- rep(NA, length(sts))
for (i in 1:length(sts)) {
	x <- infos[infos$ST == sts[i],]
	tx <- keep.tip(tA, x$Strain)
	r <- roottotip(tx, x$Date, showFig=T)
	pvals[i] <- r$pvalue
	rate[i] <- r$rate
	}

names(rate) <- sts
names(pvals) <- sts
length(which((rate>0) == TRUE))

# permutation test : sample within ST
st <- names(which((rate>0) == TRUE)) # 5 ST

nb.rate <- integer(0)
for(j in 1:1000) {
	res.st <- integer(0)
	
	for(i in 1:length(st)){
		x <- infos[infos$ST == st[i],]
		tx <- keep.tip(tA, x$Strain)
				
		root.tip <- c(diag(vcv.phylo(tx)))
		year <- c(as.numeric(a$Date[match(names(diag(vcv.phylo(tx))), a$Strain)]))
		if (j > 1) {year <- sample(year)}

		m1 <- lm(root.tip ~ year)
		mrca.age <- - coef(m1) [1] / coef(m1) [2]
		res.st <- rbind(res.st, c(length(x$Strain), coef(m1), mrca.age, mean(root.tip)))
		}

	rates <- as.numeric(res.st[res.st[,3] > 0, 3])
	calib <- as.numeric(res.st[res.st[,3] > 0, 4])

	if (length(rates) > 1) {
		nb.rate <- c(nb.rate, length(rates))
		}

	if (length(rates) == 0) {
		nb.rate <- c(nb.rate, length(rates))
		}

	print(j)
	}

hist(nb.rate, main="Number of STs with molecular clock signal", xlim=c(range(nb.rate)), breaks=5, xlab="")
abline(v=nb.rate[1], col="red", lwd=2)
ci <- quantile(nb.rate[-1], c(0.025, 0.975), na.rm=TRUE)
abline(v=ci[1], lty=2, col="blue")
abline(v=ci[2], lty=2, col="blue")

length(which(nb.rate[-1]==5))/1000*100


### Get rate of change and calibration points (mrca)
# root calibration
root.tip <- diag(vcv.phylo(tree))
year <- as.numeric(a$Date[match(names(diag(vcv.phylo(tree))), a$Strain)])

root.tip <- c(root.tip, 0)
year <- c(year, -102E6)

m0 <- lm(root.tip ~ year)
y.root <- -102E6
x.root <- m0$coef[2] 

# STs
res.st <- integer(0)
for(i in 1:length(st)){
	x <- infos[infos$ST == st[i],]
	tx <- keep.tip(tA, x$Strain)
				
	root.tip <- c(diag(vcv.phylo(tx)))
	year <- c(as.numeric(a$Date[match(names(diag(vcv.phylo(tx))), a$Strain)]))

	m1 <- lm(root.tip ~ year)
	mrca.age <- - coef(m1) [1] / coef(m1) [2]
	res.st <- rbind(res.st, c(length(x$Strain), coef(m1), mrca.age, mean(root.tip)))
	}

colnames(res.st) <- c("nb_tips", "intercept", "rate", "mrca_age", "mean_root_to_tips")
res.st <- cbind(res.st, ST=st)
res.st <- as.data.frame(res.st)

rates <- c(res.st$rate , x.root)
calib <- c(res.st$mrca_age , y.root)

data.df <- data.frame(rates, calib)


### Estimation of parameters
x2 <- as.numeric(data.df$rates)
y2 <- max(infos$Date) - as.numeric(data.df$calib)
k <- x.root

# Estimate parameters using a linear model
model.0 <- lm(log(x2) ~ y2)  
alpha.0 <- exp(coef(model.0)[1])
beta.0 <- coef(model.0)[2]

# starting parameters
start <- list(alpha = alpha.0, beta = beta.0)
upper <- list(alpha = +Inf, beta = 0)
lower <- list(alpha = -Inf, beta = -0.01)
model <- nls(x2 ~  alpha  * exp(beta * y2) + k, start = start, control = nls.control(maxiter=1000), 
lower = lower, upper = upper, algorithm = "port")
summary(model)

# plot results
plot(c(y2), c(x2), xlab = "age of calibration (years ago)", ylab="rate of change (changes/site/year)", pch=19, type="p", col="black", log="yx", ylim=c(min(x2), 0.00005))
d <- sort(c(y2, seq(0, 1000, 1), seq(0, 100000, 200), seq(min(y2), max(y2), 1000000)))
p <- coef(model)[1] * exp(coef(model)[2] * d) + k
lines(d,p, type="l", col = 'coral1', lwd = 3)
text(c(y2), c(x2), c(st,"root"), pos=3)

# estimate confidence interval
ci.a <- confint2 (model, level = 0.95, method = "asymptotic")

p.lower <- ci.a[1,1] * exp(ci.a[2,1] * d) + k
p.upper <- ci.a[1,2] * exp(ci.a[2,2] * d) + k

lines(d,p.lower, type="l", col = 'black', lwd = 3, lty = 2)
lines(d,p.upper, type="l", col = 'black', lwd = 3, lty = 2)

