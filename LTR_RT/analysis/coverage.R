# ---------------------------------------------------------------------------------#
# Analisys of distribution pattern of LTR-RTs -------------------------------------#
# Marius Hodel --------------------------------------------------------------------#
# 2019-06 -------------------------------------------------------------------------#
#----------------------------------------------------------------------------------#


# setwd("P:/mahodel")
library(ggplot2)
library(dplyr)
library(gridExtra)


for (i in 1:nlevels(cov$V1)){
  cov1 <- cov[cov$V1 == levels(cov$V1)[i],]
  plot(cov1$V3, cov1$V7, type = "l", main = levels(cov$V1)[i], ylim = c(0,1))
}

cov <- read.table("LTR.cov")
cov.cds <- read.table("LTR.cds.cov")

par(mfrow=c(2,2))
for (i in 1:nlevels(cov$V1)){
  cov1 <- head(cov[cov$V1 == levels(cov$V1)[i],], -49)
  plot(cov1$V3, cov1$V7, type = "l", main = levels(cov$V1)[i], ylim = c(0, 1))
  cov.cds1 <- head(cov.cds[cov.cds$V1 == levels(cov$V1)[i],], -49)
  plot(cov.cds1$V3, cov.cds1$V7, type = "l", main = levels(cov.cds$V1)[i], ylim = c(0, 0.05))
}



for (i in 1:nlevels(cov.cds$V1)){
  cov.cds1 <- head(cov.cds[cov.cds$V1 == levels(cov$V1)[i],], -49)
  plot(cov.cds1$V3, cov.cds1$V7, type = "l", main = levels(cov.cds$V1)[i], ylim = c(0, 0.1))
}

# Plot coverage on linkage groups

cov <- read.table("LTR.LG.cov")
names(cov) <- c("LG", "start", "end", "nr", "lencov", "lenLG", "cov")
cov <- cov[cov$lenLG == 30000000,]
cov$start <- cov$start/10^6
cov$end <- cov$end/10^6
for (i in c(1:7)){
  cov1 <- cov[cov$LG == levels(cov$LG)[i],]
  print(
  ggplot(cov1, aes(x=start, y=lencov, fill = cov)) +
    geom_bar(stat = "identity") +
    scale_fill_gradient(low = "darkred", high = "yellow")
  )
}


# Per family

cov <- read.table("LTR.LG.fam.cov")
names(cov) <- c("LG", "start", "end", "nr", "lencov", "lenLG", "cov", "fam")
cov <- cov[cov$lenLG == 30000000,]
cov$start <- cov$start/10^6
cov$end <- cov$end/10^6

# Standardize coverage
cov2 <- cov[0,]
cov2$stdcov <- rep(0,0)
for (i in c(1:30)){
  covfam <- cov[cov$fam == levels(cov$fam)[i],]
  covfam$stdcov <- (covfam$cov - mean(covfam$cov))/sd(covfam$cov)
  cov2 <- rbind(covfam, cov2)
}

# Coverage per LG and family (one plot per LG)
cov <- na.omit(cov2)
pdf("Rplots/Families/coverage.analysis2.pdf", height = 8, width = 12)
for (i in c(1:7)){
  cov1 <- cov[cov$LG == levels(cov$LG)[i],]
  print(
    ggplot(cov1, aes(x=start, y=fam)) +
      scale_fill_gradient(low = "yellow", high = "darkred") +
      geom_tile(aes(fill = stdcov)) +
      ggtitle(levels(cov$LG)[i]) +
      theme_classic() +
      xlab("Position [Mb]")
  )
}
dev.off()

# Coverage per LG and family (one plot per family)
pdf("Rplots/Families/coverage.analysis.pdf", height = 8, width = 12)
for (i in 1:nlevels(cov$fam)){
  cov1 <- cov[cov$fam == levels(cov$fam)[i],]
  if (nrow(cov1) > 0) {
    cov1$covperc <- cov1$cov*100
    print(
      ggplot(cov1, aes(y = LG, x = start+15)) +
        scale_fill_gradient(low = "mistyrose", high = "darkred") +
        geom_tile(aes(fill = covperc), height = 0.9) +
        ggtitle(levels(cov$fam)[i]) +
        theme_classic() +
        xlab("Position [Mb]") +
        ylab("Linkage group") +
        labs(fill = "Coverage [%]")
    )
  }
}
dev.off()

# Plot Castor and Camilla in one plot

covCastor <- cov[cov$fam == "Castor",]
covCamilla <- cov[cov$fam == "Camilla",]

covCastor$covperc <- covCastor$cov*100
covCamilla$covperc <- covCamilla$cov*100

covCastor <- covCastor %>% mutate(LG = factor(LG), LG = factor(LG, levels = rev(levels(LG))))
covCamilla <- covCamilla %>% mutate(LG = factor(LG), LG = factor(LG, levels = rev(levels(LG))))

p1 <- ggplot(covCamilla, aes(y = LG, x = start+15)) +
  scale_fill_gradient(low = "mistyrose", high = "darkred") +
  geom_tile(aes(fill = covperc), height = 0.9) +
  ggtitle("Camilla") +
  theme_classic() +
  xlab("Position [Mb]") +
  ylab("Linkage group") +
  labs(fill = "Coverage [%]") +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
  theme(legend.position = c(0.85,0.85))

p2 <- ggplot(covCastor, aes(y = LG, x = start+15)) +
  scale_fill_gradient(low = "mistyrose", high = "darkred") +
  geom_tile(aes(fill = covperc), height = 0.9) +
  ggtitle("Castor") +
  theme_classic() +
  xlab("Position [Mb]") +
  ylab("Linkage group") +
  labs(fill = "Coverage [%]") +
  theme(legend.position = c(0.85,0.85))
lay <- rbind(1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2)
pdf("Rplots/Families/coverage_finalplot.pdf", height = 8, width = 12)
grid.arrange(p1, p2, layout_matrix = lay)
dev.off()

# All elements together (still only big families)
fullcov <- cov %>% group_by(LG, start) %>% summarise(cov = sum(cov))
pdf("Rplots/Families/coverage.analysis3.pdf", height = 8, width = 12)
ggplot(fullcov, aes(x=start+15, y=LG)) +
  scale_fill_gradient(low = "mistyrose", high = "darkred") +
  geom_tile(aes(fill = cov*100), height = 0.9) +
  ggtitle("Proportion of LTR-RTs") +
  theme_classic() +
  xlab("Position [Mb]") +
  ylab("Linkage group") +
  labs(fill = "Coverage [%]")
dev.off()
