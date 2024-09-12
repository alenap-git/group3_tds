### PCA of population based on protein expression

# Load packages
library(tidyverse)
library(FactoMineR)
library(colorspace)

# Load data
proteins = readRDS("../Data/Proteins_nd_imputed.rds")
print(ncol(proteins))

mypca=PCA(proteins, graph = FALSE)

saveRDS(mypca, "../Data/Proteins_PCA.rds")

ev=mypca$eig[,2]/100

# Scree plot
pdf(paste0("../Figures/PCA_scree_plot.pdf"), width=5, height=5)
par(mar=c(5,5,1,1))
plot(ev, pch=19, col="#0019a8", las=1, type="b", ylim=c(0,1),
     ylab="Proportion of explained variance", xlab="PCs")
points(cumsum(ev), pch=19, col="#dc251f", type="b")
legend("right", pch=19, col=c("#0019a8", "#dc251f"),
       legend=c("Proportion of e.v.", "Cumulative proportion of e.v."))
dev.off()
