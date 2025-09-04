library(phytools)
library(TreeSim)
library(zoo)
library(coda)

# open data
Tab_count <- read.csv("Tab_count.txt", header=TRUE, sep="	")

# get sampling times
Tab_count["times_round"] <- round(Tab_count$times, 1)
s_times <- unique(Tab_count$times_round)
print(length(s_times))

# add pos id
n <- length(unique(Tab_count $nsamp_id))
Tab_count["pos_id"] <- ((Tab_count$cl_id-1)*n)+ Tab_count$nsamp_id

# parameters
nsamp <- max(Tab_count$nsamp_id)
nclade <- max(Tab_count$cl_id)
print(nclade)
print(nsamp)

# get frequency through time
count_matrix <- matrix(0, ncol=nsamp*nclade, nrow=length(s_times))
for (z in 1:length(s_times)) {
	tab_z <- Tab_count [which(Tab_count$times_round == s_times[z]),]
	count_matrix[z:length(s_times),tab_z$pos_id] <- matrix(tab_z$count, ncol=length(tab_z$pos_id), nrow=length(z:length(s_times)), byrow=TRUE)
	print(z)
}

branch_matrix <- matrix(0, ncol=nsamp*nclade, nrow=length(s_times))
for (z in 1:length(s_times)) {
	tab_z <- Tab_count [which(Tab_count$times_round == s_times[z]),]
	branch_matrix[z:length(s_times),tab_z$pos_id] <- matrix(tab_z$nb_br, ncol=length(tab_z$pos_id), nrow=length(z:length(s_times)), byrow=TRUE)
	print(z)
}

freq_tab2 <- cbind(s_times, rowSums(count_matrix, na.rm = TRUE), rowSums(branch_matrix, na.rm = TRUE), rowSums(count_matrix, na.rm = TRUE)/rowSums(branch_matrix, na.rm = TRUE))


# 95%CI
n <- dim(count_matrix)[2]
xbar <- rowMeans(count_matrix, na.rm = TRUE)/rowMeans(branch_matrix, na.rm = TRUE)
xsd <- rowMeans(apply(count_matrix, 1, sd, na.rm=TRUE)/branch_matrix)

margin <- qt(0.975, df=n-1)*xsd/sqrt(n)

lowerCI <- xbar - margin
upperCI <- xbar + margin

freq_tab2 <- cbind(freq_tab2, lowerCI, upperCI)

write.table(freq_tab2, "hpi_frequency_time.txt", quote= FALSE, row.names=FALSE, sep="	")












