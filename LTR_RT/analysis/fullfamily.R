# ---------------------------------------------------------------------------------#
# Analysis about differences between (full)families --------------------------------#
# Marius Hodel --------------------------------------------------------------------#
# 2019-06-07 ----------------------------------------------------------------------#
#----------------------------------------------------------------------------------#

setwd("P:/mahodel")

# Required packages
library(ggplot2)
library(ape)
library(dplyr)
library(gridExtra)
library(plyr)
library(RColorBrewer)
library(treemapify)

# Genome size
genome.len <- 4531731679

ltr <- read.table("LTR-RT.FL.FAM.attr.lst")
colnames(ltr) <- c("id", "motif", "tsd", "ltridentity", "superfam", "scaffold", "fam", "cmotif", "ctsd",
                   "cnesting", "cgene", "cintronic", "cutr", "protdomains", "LTRsize", "size")

# Add family instead of (sub)family
ltr$fam <- sub("\\_.*", "",ltr$fam)

# Add insertiontime from extra file
age <- read.table("agev2.lst")
colnames(age) <- c("id", "insertiontime")
ltr <- left_join(ltr, age)
ltr$insertiontime <- ltr$insertiontime/10^6

# Big families (containing more than 50 copies)
big <- as.data.frame(table(ltr$fam))
big <- as.vector(big[big$Freq > 50, 1])

big.ltr <- ltr[ltr$fam %in% big,]

## Read in (coverage) gff file
# Produced by Parse_RepeatMaskerOutput.V2 (and shell script)
ltr2 <- read.gff("LTR.copies.all.v3.gff")
ltr2$family <- gsub("ID=Trunc.", "", ltr2$attributes)
ltr2$family <- gsub("ID=Solo.", "", ltr2$family)
ltr2$family <- gsub("ID=", "", ltr2$family)
ltr2$family <- gsub("\\..*", "", ltr2$family)
ltr2$family <- sub("\\_.*", "",ltr2$family)
ltr2$family <- as.factor(ltr2$family)

ltr.cop <- ltr2[ltr2$type != "part",]
FL <- as.data.frame(table(ltr.cop$family[ltr.cop$type == "FL-LTR-RT"]))


ltr.counts <- ddply(ltr.cop, .(ltr.cop$type, ltr.cop$family), nrow)
names(ltr.counts) <- c("type", "family", "freq")

# Read in superfam file (extra file to correct falsely annotate full-length elements)
superfam <- read.table("superfamv2.lst")
names(superfam) <- c("family", "superfamily")
ltr.counts <- left_join(ltr.counts, superfam)

# Ordering
ltr.FL <- ltr.counts[ltr.counts$type == "FL-LTR-RT",]
ltr.FL <- ltr.FL[order(-ltr.FL$freq),]
pos <- ltr.FL$family[ltr.FL$freq > 100]

# Plot number of identified FL-element
ggplot(data = ltr.FL[ltr.FL$freq > 100,], aes(x = family, y = freq, fill = superfamily)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.75)) +
  scale_x_discrete(limits = pos) +
  ggtitle("FL-copies") +
  ylab("Number of copies") +
  xlab("Family") +
  labs(fill = "Superfamily")

### Coverage on full genome ####
# ------------------------------

ltr2 <- left_join(ltr2, superfam)
ltrcov <- ltr2[ltr2$type == "part" | ltr2$type == "Solo-LTR",]
ltrcov$len <- ltrcov$end - ltrcov$start

supercov <- as.data.frame(ltrcov %>% group_by(superfamily) %>% dplyr::summarise(sum = sum(len)))
supercov$superfamily <- as.character(supercov$superfamily)
supercov[4,] <- c("non-LTR", genome.len - sum(supercov$sum[1:3]))
supercov$sum <- as.numeric(supercov$sum)
supercov$cov <- supercov$sum/genome.len
ltrsum <- sum(supercov$sum[1:3])

mycols <- c("#868686FF","#0073C2FF", "darkorange3", "darkgreen")

RLplot<- ggplot(supercov, aes(x="", y=cov, fill=superfamily))+
  geom_bar(width = 1, stat = "identity", color = "black", lwd = 0.1) +
  coord_polar("y", start=0) +
  theme_void() +
  scale_fill_manual(values = mycols) +
  theme(legend.position = "left", legend.margin = margin(0,-40,0,10)) +
  labs(fill = "Superfamily") +
  guides(fill=guide_legend(nrow=4,byrow=TRUE)) +
  theme(plot.margin = margin(0,-40,0,0))

# RLX
RLXcov <- as.data.frame(ltrcov[ltrcov$superfamily == "RLX",] %>% group_by(family) %>% dplyr::summarise(sum = sum(len)))
RLXcov$cov <- RLXcov$sum/supercov$sum[supercov$superfamily == "RLX"]
RLXcov <- RLXcov[order(-RLXcov$cov),]
pos <- as.character(RLXcov$family)
RLXcov$family <- factor(RLXcov$family, levels = pos)

greens <- brewer.pal(n = 8, name = "Greens")
greens <- rev(greens)

RLXplot <- ggplot(RLXcov, aes(x="", y=cov, fill=family),  order=-as.numeric(cov)) +
  geom_bar(width = 1, stat = "identity", color = "black", lwd = 0.1) +
  ylim(c(0,1)) +
  coord_polar("y", start=0) +
  theme_void() +
  scale_fill_manual(values = greens) +
  theme(legend.position = "right", legend.margin = margin(0, -10, 0,-10)) +
  labs(fill = "RLX families") +
  guides(fill=guide_legend(ncol=3,byrow=TRUE)) +
  theme(plot.margin = margin(0,0,0,-40))


# RLG
RLGcov <- as.data.frame(ltrcov[ltrcov$superfamily == "RLG",] %>% group_by(family) %>% dplyr::summarise(sum = sum(len)))
RLGcov$cov <- RLGcov$sum/supercov$sum[supercov$superfamily == "RLG"]
part1 <- RLGcov[RLGcov$cov >= 0.005,]
part1$family <- as.character(part1$family)
nrsmallfam <- nrow(RLGcov[RLGcov$cov < 0.005,])
smallfam <- paste("Rest (", nrsmallfam, ")", sep = "")
part2 <- c(smallfam, sum(RLGcov[RLGcov$cov < 0.005,]$sum), sum(RLGcov[RLGcov$cov < 0.005,]$cov))
RLGcov <- rbind(part1, part2)
RLGcov[,2] <- as.numeric(RLGcov[,2])
RLGcov[,3] <- as.numeric(RLGcov[,3])
RLGcov <- RLGcov[order(-RLGcov$cov),]
pos <- as.character(RLGcov$family)
RLGcov$family <- factor(RLGcov$family, levels = c(pos[pos != smallfam], smallfam))

Oranges <- colorRampPalette(brewer.pal(8, "Oranges"))(nlevels(RLGcov$family))
Oranges <- rev(Oranges)

RLGplot <- ggplot(RLGcov, aes(x="", y=cov, fill=family),  order=-as.numeric(cov)) +
  geom_bar(width = 1, stat = "identity", color = "black", lwd = 0.1) +
  coord_polar("y", start=0) +
  theme_void() +
  scale_fill_manual(values = Oranges) +
  theme(legend.position = "right", legend.margin = margin(0, -10, 0,-10)) +
  labs(fill = "RLG families") +
  guides(fill=guide_legend(ncol=3,byrow=TRUE)) +
  theme(plot.margin = margin(0,0,0,-40))

# RLC
RLCcov <- as.data.frame(ltrcov[ltrcov$superfamily == "RLC",] %>% group_by(family) %>% dplyr::summarise(sum = sum(len)))
RLCcov$cov <- RLCcov$sum/supercov$sum[supercov$superfamily == "RLC"]
part1 <- RLCcov[RLCcov$cov >= 0.005,]
part1$family <- as.character(part1$family)
nrsmallfam <- nrow(RLCcov[RLCcov$cov < 0.005,])
smallfam <- paste("Rest (", nrsmallfam, ")", sep = "")
part2 <- c(smallfam, sum(RLCcov[RLCcov$cov < 0.005,]$sum), sum(RLCcov[RLCcov$cov < 0.005,]$cov))
RLCcov <- rbind(part1, part2)
RLCcov[,2] <- as.numeric(RLCcov[,2])
RLCcov[,3] <- as.numeric(RLCcov[,3])
RLCcov <- RLCcov[order(-RLCcov$cov),]
pos <- as.character(RLCcov$family)
RLCcov$family <- factor(RLCcov$family, levels = c(pos[pos != smallfam], smallfam))

Blues <- colorRampPalette(brewer.pal(8, "Blues"))(nlevels(RLCcov$family))
Blues <- rev(Blues)

RLCplot <- ggplot(RLCcov, aes(x="", y=cov, fill=family)) +
  geom_bar(width = 1, stat = "identity", color = "black", lwd = 0.1) +
  ylim(c(0,1)) +
  coord_polar("y", start=0) +
  theme_void() +
  scale_fill_manual(values = Blues) +
  theme(legend.position = "right", legend.margin = margin(0, -10, 0,-10)) +
  labs(fill = "RLC families") +
  guides(fill=guide_legend(ncol=3,byrow=TRUE)) +
  theme(plot.margin = margin(0,0,0,-40))

# put all plots together

lay <- rbind(c(NA,2,2,3,3),
             c(1,1,3,3,3),
             c(1,1,1,3,3),
             c(1,1,1,4,4),
             c(1,1,4,4,4))
pdf("Rplots/Families/fam.coverage.pdf", height = 8, width = 15)
grid.arrange(RLplot,RLXplot, RLGplot, RLCplot, layout_matrix = lay)
dev.off()

#### Tremap plot #####
#---------------------
mycols <- c("#0073C2FF", "darkorange3", "darkgreen")

cov <- as.data.frame(ltrcov %>% group_by(superfamily, family) %>% dplyr::summarise(sum = sum(len)))

treemap <- ggplot(cov, aes(area = sum, fill = superfamily, subgroup = superfamily, label = family)) +
  geom_treemap(color = "white") +
  scale_fill_manual(values = mycols) +
  geom_treemap_text(place = "centre", min.size = 8, colour = "white", padding.y = grid::unit(4, "mm"), padding.x = grid::unit(2, "mm")) +
  theme(legend.position = c(0.2,0.2))

pdf("Rplots/Families/fam.coverage.treemap.pdf", height = 8, width = 15)
treemap
dev.off()

#### Compare number of copies (LTR_retriever and Repeat_masker output) ####
# -------------------------------------------------------------------------

orig <- as.data.frame(table(ltr$fam))
names(orig) <- c("family", "cp.number.orig")

comp <- left_join(ltr.FL, orig)
plot(comp$freq, comp$cp.number.orig)
plot(comp$freq, comp$cp.number.orig, xlim = c(0,2000), ylim = c(0,500))

comp <- left_join(ltr.counts, orig) 
complex <- comp[comp$type == "FL-LTR-RT-complex",]
plot(complex$freq, complex$cp.number.orig)
plot(complex$freq, complex$cp.number.orig, xlim = c(0,2000), ylim = c(0,2000))


## Number of bp that get masked (proportion)
genome.len <- 4531731679 

ltrpart <- ltr2[ltr2$type == "part" | ltr2$type == "Solo-LTR",]

proportion <- as.data.frame(matrix(nrow = 0, ncol = 2))
names(proportion) <- c("family", "proportion")

for (i in 1:nlevels(ltrpart$family)){
  tmp <- ltrpart[ltrpart$family == levels(ltrpart$family)[i],]
  tmp$len <- tmp$end - tmp$start
  proportion[i,] <- NA
  proportion$family[i] <- levels(ltrpart$family)[i]
  proportion$proportion[i]<- sum(tmp$len)/genome.len*100
}

proportion <- left_join(proportion, superfam)
proportion <- proportion[order(-proportion$proportion),]
pos <- proportion$family[proportion$proportion > 0.1]

# Plot proportion
ggplot(data = proportion[proportion$proportion > 0.1,], aes(x = family, y = proportion, fill = superfamily)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.75)) +
  scale_x_discrete(limits = pos) +
  ggtitle("Proportion of the assembly") +
  ylab("Proportion [%]")


ltr <- left_join(ltr, superfam, by = c("fam" = "family"))
ltr <- ltr[,-5]

big.ltr <- ltr[ltr$fam %in% big,]

ggplot(big.ltr, aes(x=fam, y=insertiontime, fill = superfamily)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5))

### Lifespan / ltridentity (plot violin plot ####
# -----------------------------------------------

mycols <- c("#0073C2FF", "darkorange3", "darkgreen")

# Calculate median
median.big.ltr <- big.ltr %>% group_by(fam) %>% dplyr::summarise(median = median(ltridentity))
# Ordering
median.big.ltr <- median.big.ltr[order(-median.big.ltr$median),]
pos <- median.big.ltr$fam

p1 <- ggplot(big.ltr, aes(x=fam, y=ltridentity)) +
  geom_violin(aes(fill = superfamily), linetype = "blank") +
  theme_minimal() +
  scale_x_discrete(limits = pos) +
  coord_flip() +
  xlab("family") +
  ylab("LTR identity [%]") +
  theme(legend.position = c(0.15,0.85)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(axis.line=element_blank(), axis.ticks.x = element_line())

fam.size <- as.data.frame(table(big.ltr$fam))
names(fam.size) <- c("fam", "size")
p2 <- ggplot(data = fam.size, aes(y = size, x = fam)) +
  geom_bar(stat = "identity") +
  scale_x_discrete(limits = pos) +
  coord_flip() +
  theme_minimal() +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(axis.line=element_blank(),axis.text.y=element_blank(), axis.title.y=element_blank(), 
        axis.ticks.x = element_line()) +
  ylab("Number of copies")

lay <- rbind(c(1,1,1,1,2))

grid.arrange(p1, p2, layout_matrix = lay)
big.ltr$fam <- as.factor(big.ltr$fam)
max.big.ltr <- big.ltr %>% group_by(fam) %>% dplyr::summarise(max = max(insertiontime), min = min(insertiontime), median = median(insertiontime))

# Ordering
max.big.ltr <- max.big.ltr[order(max.big.ltr$median),]
pos <- max.big.ltr$fam

# Plot lifespan per family (violin plot)

median.quartile <- function(x){
  out <- quantile(x, probs = c(0.05,0.5,0.95))
  names(out) <- c("ymin","y","ymax")
  return(out) 
}

p1 <- ggplot(big.ltr, aes(x=fam, y=insertiontime)) +
  geom_violin(aes(fill = superfamily), linetype = "blank") +
  stat_summary(fun.data = median.quartile, geom = "pointrange", color = "black", size = 0.25, linetype = "dashed", fatten = 1) +
  theme_minimal() +
  scale_x_discrete(limits = pos) +
  coord_flip() +
  xlab("Family") +
  ylab("Insertion time [myr]") +
  theme(legend.position = c(0.15,0.85)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(axis.line=element_blank(), axis.ticks.x = element_line()) +
  scale_y_reverse(c(0,1)) +
  ylim(c(2.7,0)) +
  scale_fill_manual(values = mycols)

p2 <- ggplot(data = fam.size, aes(y = size, x = fam)) +
  geom_bar(stat = "identity", fill = "black") +
  scale_x_discrete(limits = pos) +
  coord_flip() +
  theme_minimal() +
  theme(axis.line=element_blank(),axis.text.y=element_blank(), axis.title.y=element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(axis.line=element_blank(),axis.text.y=element_blank(), axis.title.y=element_blank(), 
        axis.ticks.x = element_line()) +
  geom_text(aes(label=size), position=position_dodge(width=0.9), hjust=-0.5, size = 2, vjust = 0.25) +
  annotate("text", x = 27, y = 10500, label="15670", color = "white",  hjust=-0.5, size = 2, vjust = 0.25) +
  scale_y_continuous(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
  ylab("Number of flLTR-RTs")

lay <- rbind(c(1,1,1,1,2))

pdf("Rplots/Families/Lifespanfullfam2.pdf", width = 9, height = 5.85)
grid.arrange(p1, p2, layout_matrix = lay)
dev.off()

# Plot with coverage instead of number of copies
fam.size <- left_join(fam.size, proportion, by = c("fam" = "family"))

p2 <- ggplot(data = fam.size, aes(y = proportion, x = fam)) +
  geom_bar(stat = "identity", fill = "black") +
  scale_x_discrete(limits = pos) +
  coord_flip() +
  theme_minimal() +
  theme(axis.line=element_blank(),axis.text.y=element_blank(), axis.title.y=element_blank()) +
  ylab("Assembly coverage") +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(axis.line=element_blank(),axis.text.y=element_blank(), axis.title.y=element_blank(), 
        axis.ticks.x = element_line())

grid.arrange(p1, p2, layout_matrix = lay)

# Coverage vs. Number of clean FL copies
fam.size$cpnrvscoverage <- fam.size$size/fam.size$proportion

p2 <- ggplot(data = fam.size, aes(y = cpnrvscoverage, x = fam)) +
  geom_bar(stat = "identity", fill = "black") +
  scale_x_discrete(limits = pos) +
  coord_flip() +
  theme_minimal() +
  theme(axis.line=element_blank(),axis.text.y=element_blank(), axis.title.y=element_blank()) +
  ylab("Assembly coverage") +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(axis.line=element_blank(),axis.text.y=element_blank(), axis.title.y=element_blank(), 
        axis.ticks.x = element_line())
grid.arrange(p1, p2, layout_matrix = lay)

# Coverage clean FL vs coverage Repeatmasker
bp <- big.ltr %>% group_by(fam) %>% dplyr::summarise(bp = sum(size))
fam.size <- left_join(fam.size, bp)
fam.size$bp2 <- fam.size$proportion * genome.len
fam.size$comp <- fam.size$bp2/fam.size$bp

p2 <- ggplot(data = fam.size, aes(y = comp, x = fam)) +
  geom_bar(stat = "identity", fill = "black") +
  scale_x_discrete(limits = pos) +
  coord_flip() +
  theme_minimal() +
  theme(axis.line=element_blank(),axis.text.y=element_blank(), axis.title.y=element_blank()) +
  ylab("Assembly coverage") +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(axis.line=element_blank(),axis.text.y=element_blank(), axis.title.y=element_blank(), 
        axis.ticks.x = element_line())

#### Relationship between cv of length and age ####
# -------------------------------------------------

# Calculate CV
# Read in TE size
for (i in 1:nlevels(big.ltr$fam)) {
  fam.size$sdlen[fam.size$fam == levels(big.ltr$fam)[i]] <- sd(big.ltr[big.ltr$fam == levels(big.ltr$fam)[i],]$size)/mean(big.ltr[big.ltr$fam == levels(big.ltr$fam)[i],]$size)
  fam.size$meanage[fam.size$fam == levels(big.ltr$fam)[i]]<- mean(big.ltr[big.ltr$fam == levels(big.ltr$fam)[i],]$insertiontime)
}

plot(fam.size$sdlen, fam.size$meanage)
agevssdlen <- lm(meanage ~ sdlen, data = fam.size)
summary(agevssdlen)

# Plot with ggplot
int <- coefficients(agevssdlen)[1]
slo <- coefficients(agevssdlen)[2]
pdf("Rplots/Families/TEsize_VS_Age.pdf", height = 5, width = 5)
ggplot(data = fam.size, aes(y = meanage, x = sdlen)) +
  geom_point() +
  xlab("Coefficient of variation of TE size per family") +
  ylab("Mean insertion time per family [mya]") +
  geom_abline(intercept = int, slope = slo, linetype = "dashed") +
  theme_classic() +
  geom_text(x = 0.1, y = 0.8, label = paste("R^2 == ", round(summary(agevssdlen)$r.squared,2) , sep = ""), parse = T)
dev.off()

# Check with median
for (i in 1:nlevels(big.ltr$fam)) {
  fam.size$sdlen[fam.size$fam == levels(big.ltr$fam)[i]] <- sd(big.ltr[big.ltr$fam == levels(big.ltr$fam)[i],]$size)/mean(big.ltr[big.ltr$fam == levels(big.ltr2$fam)[i],]$size)
  fam.size$medianage[fam.size$fam == levels(big.ltr$fam)[i]]<- median(big.ltr[big.ltr$fam == levels(big.ltr$fam)[i],]$insertiontime)
}

plot(fam.size$sdlen, fam.size$medianage)
agevssdlen <- lm(sdlen ~ medianage, data = fam.size)
summary(agevssdlen)

#### Histogram of TE size len ####
#---------------------------------

for (i in levels(big.ltr$fam)){
  tmp <- big.ltr[big.ltr$fam == i,]
  if (nrow(tmp) > 500){
    message(paste(i,nrow(tmp), sep =" "))
    tmp$subfamily <- gsub("[0-9]+", "", tmp$id)
    tmp$subfamily <- gsub("_$", "", tmp$subfamily)
    nam <- paste("gg", i, sep = "_")
    p <- ggplot(data = tmp, aes(size, fill = cnesting)) +
      geom_histogram(binwidth = 100) +
      annotate("text", label = i, x = 4000, y = Inf, vjust = 2, size = 5) +
      xlab("TE size [nt]") +
      ylab("Copy number") +
      scale_fill_discrete(name = "", labels = c("Not nested", "Insert", "Host")) +
      theme_bw() +
      scale_x_continuous(labels = function(x) format(x, big.mark = ",", scientific = FALSE), limits = c(2000, 23000)) +
      theme(legend.position = c(0.8, 0.7))
    assign(nam, p)
    print(p)
  }
}

gg_Dom <- gg_Dom +
  theme(axis.title.y = element_blank(), legend.position = "blank")
gg_Parrot <- gg_Parrot +
  theme(axis.title.y = element_blank(), legend.position = "blank")
lay <- rbind(c(1,2,3),c(1,2,3))
pdf("Rplots/Families/elementSize.pdf", height = 5, width = 11)
grid.arrange(gg_Camilla, gg_Dom, gg_Parrot, layout_matrix = lay)
dev.off()

#### Test relationship between categorical variables (Collocalization) ####
# -------------------------------------------------------

# Binominal tests for intron, utr and genes
features <- as.data.frame(matrix(nrow = 29, ncol = 4))
tabfeat2 <- data.frame(matrix(nrow = 29, ncol = 0))
expected <- vector()
for (j in 1:3) {
  tabfeat  <- table(big.ltr$fam, big.ltr[,j+9])
  tabfeat2 <- cbind(tabfeat2, as.data.frame(tabfeat[,2]/tabfeat[,1])[,1])
  expected[j] <- (sum(tabfeat[,2])/sum((tabfeat[,1]+tabfeat[,2]))*100)
  for (i in 1:nrow(tabfeat)){
    bi <- binom.test(x = tabfeat[i,2], tabfeat[i,1]+tabfeat[i,2], p = sum(tabfeat[,2])/sum((tabfeat[,1]+tabfeat[,2])),
                     alternative = "two.sided")
    features[i,1] <- rownames(tabfeat)[i]
    features[i,j+1]<- bi$p.value
  }
}
tabfeat2 <- cbind(rownames(tabfeat), tabfeat2)
colnames(tabfeat2) <- c("family", "propgene", "propintron", "proputr")
tabfeat2 <- left_join(tabfeat2, superfam)
colnames(features) <- c("family", "p.gene", "p.intron", "p.utr")

# Add asteriks to dataframe
features <- cbind(features, as.data.frame(matrix(nrow = 29, ncol = 3)))
colnames(features)[5:7] <- c("a.gene", "a.intron", "a.utr")
for (i in 2:4){
  features[features[,i] < 0.05,i+3] <- "*"
  features[features[,i] < 0.01,i+3] <- "**"
  features[features[,i] < 0.001,i+3] <- "***"
}


# Plot results
# Dont plot p1
p1 <- ggplot(data = tabfeat2, aes(x = family, y = propgene*100, fill = superfamily)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5)) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  ylab("TEs containing a gene [%]") +
  xlab("Family") +
  annotate("text", x = features$family, y = tabfeat2$propgene*100+max(tabfeat2$propgene)*10, label = features$a.gene) +
  theme(axis.text.x=element_blank(), axis.title.x=element_blank(), 
        axis.ticks.x = element_line()) +
  theme(legend.position = c(0.1, 0.7)) +
  scale_fill_manual(values = mycols)


p2 <- ggplot(data = tabfeat2, aes(x = family, y = propintron*100, fill = superfamily)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.25)) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  ylab("Insertions located inside introns [%]") +
  xlab("Family") +
  annotate("text", x = features$family, y = tabfeat2$propintron*100+max(tabfeat2$propintron)*10, label = features$a.intron) +
  theme(axis.text.x=element_blank(), axis.title.x=element_blank(), 
        axis.ticks.x = element_line()) +
  theme(legend.position=c(0.1, 0.7)) +
  scale_fill_manual(values = mycols) +
  geom_hline(yintercept = expected[2], linetype = "dashed", size = 0.25)
  

p3 <- ggplot(data = tabfeat2, aes(x = family, y = proputr*100, fill = superfamily)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 0.25)) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  ylab("Insertions overlapping with UTRs of genes [%]") +
  xlab("Family") +
  annotate("text", x = features$family, y = tabfeat2$proputr*100+max(tabfeat2$proputr)*10, label = features$a.utr) +
  theme(legend.position="none") +
  theme(axis.title.x = element_text(margin = margin(t = 10)), axis.text.x = element_text(margin = margin(t = -10))) +
  scale_fill_manual(values = mycols) +
  geom_hline(yintercept = expected[3], linetype = "dashed", size = 0.25)


lay <- rbind(c(2,2),c(2,2),c(2,2),c(2,2),c(2,2),c(3,3),c(3,3),c(3,3),c(3,3),c(3,3),c(3,3))
pdf("Rplots/Families/fullfamily.features.pdf", height = 8.27, width = 11.69)
grid.arrange(p2, p3, layout_matrix = lay)
dev.off()

# Test same differences between superfamilies
features <- as.data.frame(matrix(nrow = 3, ncol = 4))
tabfeat2 <- data.frame(matrix(nrow = 3, ncol = 0))
expected <- vector()
for (j in 1:3) {
  tabfeat  <- table(ltr$superfamily, ltr[,j+9])
  tabfeat2 <- cbind(tabfeat2, as.data.frame(tabfeat[,2]/tabfeat[,1])[,1])
  for (i in 1:nrow(tabfeat)){
    bi <- binom.test(x = tabfeat[i,2], tabfeat[i,1]+tabfeat[i,2], p = sum(tabfeat[,2])/sum((tabfeat[,1]+tabfeat[,2])),
                     alternative = "two.sided")
    print(bi)
    features[i,1] <- rownames(tabfeat)[i]
    features[i,j+1]<- bi$p.value
    expected[j] <- bi$null.value
  }
}
tabfeat2 <- cbind(rownames(tabfeat), tabfeat2)
colnames(tabfeat2) <- c("superfamily", "propgene", "propintron", "proputr")
colnames(features) <- c("superfamily", "p.gene", "p.intron", "p.utr")

# Add asteriks to dataframe
features <- cbind(features, as.data.frame(matrix(nrow = 3, ncol = 3)))
colnames(features)[5:7] <- c("a.gene", "a.intron", "a.utr")
for (i in 2:4){
  features[features[,i] < 0.05,i+3] <- "*"
  features[features[,i] < 0.01,i+3] <- "**"
  features[features[,i] < 0.001,i+3] <- "***"
}

p4 <- ggplot(data = tabfeat2, aes(x = superfamily, y = propintron*100, fill = superfamily)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  xlab("Superfamily") +
  annotate("text", x = features$superfamily, y = tabfeat2$propintron*100+max(tabfeat2$propintron)*10, label = features$a.intron) +
  theme(legend.position="none", axis.text.x=element_blank(), axis.title.x=element_blank(), 
        axis.ticks.x = element_line(), axis.title.y=element_blank()) +
  scale_fill_manual(values = mycols) +
  geom_hline(yintercept = expected[2]*100, linetype = "dashed", size = 0.25)

p5 <- ggplot(data = tabfeat2, aes(x = superfamily, y = proputr*100, fill = superfamily)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  xlab("Superfamily") +
  annotate("text", x = features$superfamily, y = tabfeat2$proputr*100+max(tabfeat2$proputr)*10, label = features$a.utr) +
  theme(legend.position="none", axis.title.y=element_blank()) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 0.75)) +
  scale_fill_manual(values = mycols) +
  geom_hline(yintercept = expected[3]*100, linetype = "dashed", size = 0.25)

lay <- rbind(c(2,2,2,2,4),c(2,2,2,2,4),c(2,2,2,2,4),c(2,2,2,2,4),c(2,2,2,2,4),c(3,3,3,3,5),c(3,3,3,3,5),c(3,3,3,3,5),c(3,3,3,3,5),c(3,3,3,3,5))
pdf("Rplots/Families/superfamily_family_features.pdf", height = 8.27, width = 11.69)
grid.arrange(p2, p3 ,p4, p5, layout_matrix = lay)
dev.off()

#### Table of features (comparisons between superfamilies) ####
# -------------------------------------------------------------

# (Final table created by hand, but data was analysed here)
# Change superfamily according to family
ltr <- left_join(ltr, superfam, by = c("fam" = "family"))
ltr <- ltr[,-5]
ltr$LTRsize <- as.numeric(ltr$LTRsize)
ltr$size <- as.numeric(ltr$size)


# Nuber of flElements
nr <- table(ltr$superfamily)

# Mean and sd of LTRsize
#mean
ltr %>%
  group_by(superfamily) %>%
  select_if(is.numeric) %>% 
  summarise_at(vars(-ltridentity),funs(mean(., na.rm=TRUE))) %>% select(LTRsize)
#sd
ltr %>%
  group_by(superfamily) %>%
  select_if(is.numeric) %>% 
  summarise_at(vars(-ltridentity),funs(sd(., na.rm=TRUE))) %>% select(LTRsize)

# Mean and sd of LTRsize
#mean
ltr[ltr$cnesting != "N",] %>%
  group_by(superfamily) %>%
  select(size) %>%
  summarise_at(vars(size),funs(mean(., na.rm=TRUE)))
#sd
ltr[ltr$cnesting != "N",] %>%
  group_by(superfamily) %>%
  select(size) %>%
  summarise_at(vars(size),funs(sd(., na.rm=TRUE)))

# Autonomous elements
prop.table(table(ltr[ltr$superfamily == "RLG",]$protdomains, useNA="always"))
prop.table(table(ltr[ltr$superfamily == "RLC",]$protdomains, useNA="always"))

# Motif
ltr %>%
  group_by(superfamily) %>%
  select(cmotif) %>% 
  table
# Calculate percentage by hand ...

# TSD
ltr %>%
  group_by(superfamily) %>%
  select(ctsd) %>% 
  table
# Calculate percentage by hand ...

#### Check same stuff for biggest families in each superfamily ####
# -----------------------------------------------------------------

head(sort(table(ltr[ltr$superfamily == "RLC",]$fam), decreasing = T),3)
head(sort(table(ltr[ltr$superfamily == "RLG",]$fam), decreasing = T),3)
head(sort(table(ltr[ltr$superfamily == "RLX",]$fam), decreasing = T),3)

# mean LTR size
ltr %>%
  group_by(fam) %>%
  select_if(is.numeric) %>% 
  summarise_at(vars(-ltridentity),funs(mean(., na.rm=TRUE))) %>% 
  filter(fam == "Camilla" | fam == "Dom" | fam == "Lis"
         | fam == "Tasch" | fam == "Leojyg" | fam == "Roseg" 
         | fam == "Monch" | fam == "Allalin" | fam == "Dirru")
# SD LTR size
ltr %>%
  group_by(fam) %>%
  select_if(is.numeric) %>% 
  summarise_at(vars(-ltridentity),funs(sd(., na.rm=TRUE))) %>% 
  filter(fam == "Camilla" | fam == "Dom" | fam == "Lis" 
         | fam == "Tasch" | fam == "Leojyg" | fam == "Roseg" 
         | fam == "Monch" | fam == "Allalin" | fam == "Dirru")

# mean TE size
ltr[ltr$cnesting != "N",] %>%
  group_by(fam) %>%
  select(size) %>%
  summarise_at(vars(size),funs(mean(., na.rm=TRUE))) %>% 
  filter(fam == "Camilla" | fam == "Dom" | fam == "Lis" 
         | fam == "Tasch" | fam == "Leojyg" | fam == "Roseg" 
         | fam == "Monch" | fam == "Allalin" | fam == "Dirru")
# sd TE size
ltr[ltr$cnesting != "N",] %>%
  group_by(fam) %>%
  select(size) %>%
  summarise_at(vars(size),funs(sd(., na.rm=TRUE))) %>% 
  filter(fam == "Camilla" | fam == "Dom" | fam == "Lis" 
         | fam == "Tasch" | fam == "Leojyg" | fam == "Roseg" 
         | fam == "Monch" | fam == "Allalin" | fam == "Dirru")

big.fam.list <- c("Tasch", "Leojyg", "Roseg", "Camilla", "Dom", "Lis", "Monch", "Allalin", "Dirru")

#Protein domains
for (i in big.fam.list) {
  message(i)
  print(
  prop.table(table(ltr[ltr$fam == i,]$protdomains, useNA="always")))
}

# Motif
for (i in big.fam.list) {
  message(i)
  print(
  ltr %>%
    group_by(fam) %>%
    select(cmotif) %>% 
    filter(fam == i) %>% 
    table %>%
    prop.table)
}

for (i in big.fam.list) {
  message(i)
  print(
    ltr %>%
      group_by(fam) %>%
      select(ctsd) %>% 
      filter(fam == i) %>% 
      table %>%
      prop.table)
}
