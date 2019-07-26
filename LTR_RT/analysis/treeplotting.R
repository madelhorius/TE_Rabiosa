# ---------------------------------------------------------------------------------#
# Phylogenetic Trees --------------------------------------------------------------#
# Visualization of clustalw2 output -----------------------------------------------#
# Marius Hodel --------------------------------------------------------------------#
# 2019-06-07-----------------------------------------------------------------------#
#----------------------------------------------------------------------------------#

setwd("P:/mahodel")

# Packages

library(ggtree)
library(ggplot2)
library(phytools)
library(gridExtra)

# Color palette
col <- c("#00b34e","#0235da","#72e91b","#f900cb","#71ac00","#c276ff","#fbff66","#002897","#ffc843","#0188f5","#ff3616","#01d6ee","#ff4146",
         "#8effdb","#f80083","#006d1c","#8f0079","#018972","#a30050","#d8edff","#90000a","#eab2ff","#001a00","#ff8053","#002650","#a26200",
         "#190034","#ffe6d7","#003c37", "#583600", "black")

# Gypsy and Copia combined
# Biggest families
# ------------------------

tree1 <- read.tree("tree1.ph")

groupInfo <- split(tree1$tip.label, as.data.frame(matrix(na.omit(unlist(strsplit(as.character(tree1$tip.label), "_"))), ncol = 3, byrow = T))[,3])
tree1 <- groupOTU(tree1, groupInfo)

ggtr <- ggtree(tree1, layout = "circular") +
  theme(legend.position = "bottom") +
  geom_tippoint(aes(color=group)) +
  scale_color_manual(values = col) +
  ggtitle(paste("TREE1")) +
  geom_text2(aes(subset = !isTip, label=node), hjust=1.1, vjust = -0.3)
print(ggtr)

# Reroot tree (so that Gypsy and Copia are in different branches)

position <- 0.5*tree1$edge.length[which(tree1$edge[,2]== 1254)]
tree1.reroot <- reroot(tree = tree1, node.number = 1254, position = position)

groupInfo <- split(tree1.reroot$tip.label, as.data.frame(matrix(na.omit(unlist(strsplit(as.character(tree1.reroot$tip.label), "_"))), ncol = 3, byrow = T))[,3])
tree1.reroot <- groupOTU(tree1.reroot, groupInfo)

ggtr <- ggtree(tree1.reroot, layout = "fan", open.angle = 120) +
  theme(legend.position = c(0.5,0.8)) +
  geom_tippoint(aes(color=group)) +
  scale_color_manual(values = col, name = "Family") +
  geom_cladelabel(node=1269, label="Copia", align=T, color='red', offset.text = 0.1) +
  geom_cladelabel(node=709, label="Gypsy", align=T, color='blue', offset.text = 0.01) +
  guides(color=guide_legend(nrow=4,byrow=TRUE))
pdf("Rplots/Trees/final.fam.tree.both.pdf", width = 10, height = 10)
print(rotate_tree(ggtr, 150))
dev.off()

# Additional tree
# ---------------

tree2 <- read.tree("tree2.trep.ph")
lst <- strsplit(as.character(tree2$tip.label), "_")
df <- as.data.frame(matrix(nrow = 0, ncol = 5))
for (i in 1:length(lst)){
  if (length(lst[[i]]) == 5) {
    df[nrow(df)+1,] <- as.vector(lst[[i]])
  } else {
    df[nrow(df)+1,] <- as.vector(c(lst[[i]][1], lst[[i]]))
  }
}

groupInfo <- split(tree2$tip.label, df[,5])
tree2 <- groupOTU(tree2, groupInfo)

ggtr <- ggtree(tree2, layout = "circular") +
  theme(legend.position = "bottom") +
  geom_tippoint(aes(color=group)) +
  scale_color_manual(values = col) +
  ggtitle(paste("TREE2")) +
  geom_text2(aes(subset = !isTip, label=node), hjust=1.1, vjust = -0.3)
print(ggtr)

position <- 0.5*tree2$edge.length[which(tree2$edge[,2]== 797)]
tree2.reroot <- reroot(tree = tree2, node.number = 797, position = position)
lst <- strsplit(as.character(tree2.reroot$tip.label), "_")
df <- as.data.frame(matrix(nrow = 0, ncol = 5))
for (i in 1:length(lst)){
  if (length(lst[[i]]) == 5) {
    df[nrow(df)+1,] <- as.vector(lst[[i]])
  } else {
    df[nrow(df)+1,] <- as.vector(c(lst[[i]][1], lst[[i]]))
  }
}

groupInfo <- split(tree2.reroot$tip.label, df[,4])
tree2.reroot <- groupOTU(tree2.reroot, groupInfo)

ggtr <- ggtree(tree2.reroot, layout = "circular") +
  theme(legend.position = "bottom") +
  geom_tippoint(aes(color=group)) +
  scale_color_manual(values = col) +
  ggtitle(paste("TREE2"))
print(ggtr)

# Gypsy tree
# ----------

tree1 <- read.tree("RLG.subsample.RT.INT.domain.bigfam.ph")
tree1$tip.label <- gsub("_[a-zA-Z]_", "_", tree1$tip.label)
groupInfo <- split(tree1$tip.label, as.data.frame(matrix(na.omit(unlist(strsplit(as.character(tree1$tip.label), "_"))), ncol = 4, byrow = T))[,3])
tree1 <- groupOTU(tree1, groupInfo)

ggtr <- ggtree(tree1, layout = "circular") +
  theme(legend.position = "bottom") +
  geom_tippoint(aes(color=group)) +
  scale_color_manual(values = col) +
  ggtitle(paste("TREE1")) +
  geom_text2(aes(label=node), hjust=1.1, vjust = -0.3)
print(ggtr)

# Reroot tree (according to RLC outgroup)

position <- 0.01*tree1$edge.length[which(tree1$edge[,2]== 1)]
tree1.reroot <- reroot(tree = tree1, node.number = 1, position = position)

groupInfo <- split(tree1.reroot$tip.label, as.data.frame(matrix(na.omit(unlist(strsplit(as.character(tree1.reroot$tip.label), "_"))), ncol = 4, byrow = T))[,3])
tree1.reroot <- groupOTU(tree1.reroot, groupInfo)

ggtr <- ggtree(tree1.reroot, layout = "fan", open.angle = 120) +
  theme(legend.position = c(0.5,0.8)) +
  geom_tippoint(aes(color=group), size = 1) +
  scale_color_manual(values = col, name = "Family") +
  geom_text2(aes(label=node), hjust=1.1, vjust = -0.3) +
  geom_tiplab2()
print(rotate_tree(ggtr, 150))

RLG.tree1.reroot <- tree1.reroot

# Add Neumann lineages to tree (Respective nodes were determined using blastp)

RLG.tree <- ggtree(RLG.tree1.reroot) +
  geom_tippoint(aes(color=group), size = 1) +
  scale_color_manual(values = col, name = "Family") +
  geom_hilight(node=884, fill="firebrick4", alpha = 0.25) + 
  geom_hilight(node=848, fill="seagreen", alpha = 0.25) +
  geom_hilight(node = 898, fill = "blue", alpha = 0.25) +
  geom_hilight(node = 982, fill = "deeppink", alpha = 0.25) +
  geom_hilight(node = 1159, fill = "aquamarine4", alpha = 0.25) +
  geom_hilight(node = 843, fill = "darkgoldenrod3", alpha = 0.25) +
  geom_cladelabel(node=884, label='paste(italic("Reina"))', geom = "text", color = c("white", "black"), fontsize = 5, parse = T) + 
  geom_cladelabel(node=848, label='paste(italic("Tekay"))', geom = "text", color = c("white", "black"), fontsize = 5, parse = T) +
  geom_cladelabel(node=898, label='paste(italic("CRM"))', geom = "text", color = c("white", "black"), fontsize = 5, parse = T) +
  geom_cladelabel(node=982, label='paste(italic("Athila"))', geom = "text", color = c("white", "black"), fontsize = 5, parse = T) +
  geom_cladelabel(node=1159, label='paste(italic("Ogre"))', geom = "text", color = c("white", "black"), fontsize = 5, parse = T) +
  geom_cladelabel(node=843, label='paste(italic("Retand"))', geom = "text", color = c("white", "black"), fontsize = 5, parse = T) +
  geom_cladelabel(node=1, label="RLC", geom = "text", color = c("white", "black"), hjust = 1, vjust = 1.5, fontsize = 5, parse = T) +
  theme(legend.position = "none") +
  ggtitle("Gypsy") +
  theme(plot.title = element_text(margin = margin(t = 10, b = -40), size = 25, face = "italic"))



# Copia tree
# ----------

tree1 <- read.tree("RLC.subsample.RT.INT.domain.bigfam.ph")
tree1$tip.label <- gsub("_[a-zA-Z]_", "_", tree1$tip.label)
groupInfo <- split(tree1$tip.label, as.data.frame(matrix(na.omit(unlist(strsplit(as.character(tree1$tip.label), "_"))), ncol = 4, byrow = T))[,3])
tree1 <- groupOTU(tree1, groupInfo)

ggtr <- ggtree(tree1, layout = "circular") +
  theme(legend.position = "bottom") +
  geom_tippoint(aes(color=group)) +
  scale_color_manual(values = col) +
  ggtitle(paste("TREE1"))
print(ggtr)

# Reroot tree (according to RLG outgroup)

position <- 0.1*tree1$edge.length[which(tree1$edge[,2]== 1)]
tree1.reroot <- reroot(tree = tree1, node.number = 1, position = position)

groupInfo <- split(tree1.reroot$tip.label, as.data.frame(matrix(na.omit(unlist(strsplit(as.character(tree1.reroot$tip.label), "_"))), ncol = 4, byrow = T))[,3])
tree1.reroot <- groupOTU(tree1.reroot, groupInfo)

ggtr <- ggtree(tree1.reroot, layout = "fan", open.angle = 120) +
  theme(legend.position = c(0.5,0.8)) +
  geom_tippoint(aes(color=group), size = 1) +
  scale_color_manual(values = col, name = "Family") +
  geom_text2(aes(label=node), hjust=1.1, vjust = -0.3) +
  geom_tiplab2()
print(rotate_tree(ggtr, 150))
ggtr

ggtree(tree1.reroot)
RLC.tree1.reroot <- tree1.reroot

# Add Neumann lineages to tree (respective nodes were determined using blastp)

RLC.tree <- ggtree(RLC.tree1.reroot) +
  geom_tippoint(aes(color=group), size = 1) +
  scale_color_manual(values = col, name = "Family") +
  geom_hilight(node=900, fill="gold", alpha = 0.25) + 
  geom_hilight(node=828, fill="purple", alpha = 0.25) +
  geom_hilight(node = 909, fill = "grey", alpha = 0.25) +
  geom_hilight(node = 1100, fill = "red", alpha = 0.25) +
  geom_hilight(node = 1130, fill = "blue", alpha = 0.25) +
  geom_hilight(node = 1202, fill = "orange", alpha = 0.25) +
  geom_hilight(node = 822, fill = "pink", alpha = 0.25) +
  geom_cladelabel(node=900, label='paste(italic("Ivana"))', geom = "text", color = c("white", "black"), fontsize = 5, parse = T) + 
  geom_cladelabel(node=828, label='paste(italic("SIRE"))', geom = "text", color = c("white", "black"), fontsize = 5, parse = T) +
  geom_cladelabel(node=909, label='paste(italic("Ale"))', geom = "text", color = c("white", "black"), fontsize = 5, parse = T) +
  geom_cladelabel(node=1100, label='paste(italic("Tork"))', geom = "text", color = c("white", "black"), fontsize = 5, parse = T) +
  geom_cladelabel(node=1130, label='paste(italic("TAR"))', geom = "text", color = c("white", "black"), fontsize = 5, parse = T) +
  geom_cladelabel(node=1202, label='paste(italic("Ikeros"))', geom = "text", color = c("white", "black"), fontsize = 5, parse = T) +
  geom_cladelabel(node=822, label='paste(italic("Angela"))', geom = "text", color = c("white", "black"), fontsize = 5, parse = T) +
  geom_cladelabel(node=1, label="RLG", geom = "text", color = c("white", "black"), hjust = 1, vjust = 1.5, fontsize = 5, parse = T) +
  theme(legend.position = "none") +
  ggtitle('Copia') + 
  theme(plot.title = element_text(margin = margin(t = 10, b = -40), size = 25, face = "italic"))

pdf("Rplots/Trees/finalTree.pdf", height = 16.5, width = 11.7)
grid.arrange(RLC.tree, RLG.tree, layout_matrix = rbind(c(1,2),c(1,2)))
dev.off()

# Plots for ppt
RLC.pres <-  ggtree(RLC.tree1.reroot) +
  coord_flip() +
  geom_tippoint(aes(color=group), size = 1) +
  scale_color_manual(values = col, name = "Family") +
  geom_hilight(node=900, fill="gold", alpha = 0.25) + 
  geom_hilight(node=828, fill="purple", alpha = 0.25) +
  geom_hilight(node = 909, fill = "grey", alpha = 0.25) +
  geom_hilight(node = 1100, fill = "red", alpha = 0.25) +
  geom_hilight(node = 1130, fill = "blue", alpha = 0.25) +
  geom_hilight(node = 1202, fill = "orange", alpha = 0.25) +
  geom_hilight(node = 822, fill = "pink", alpha = 0.25) 
RLG.pres <- ggtree(RLG.tree1.reroot) +
  geom_tippoint(aes(color=group), size = 1) +
  scale_color_manual(values = col, name = "Family") +
  geom_hilight(node=884, fill="firebrick4", alpha = 0.25) + 
  geom_hilight(node=848, fill="seagreen", alpha = 0.25) +
  geom_hilight(node = 898, fill = "blue", alpha = 0.25) +
  geom_hilight(node = 982, fill = "deeppink", alpha = 0.25) +
  geom_hilight(node = 1159, fill = "aquamarine4", alpha = 0.25) +
  geom_hilight(node = 843, fill = "darkgoldenrod3", alpha = 0.25) +
  coord_flip()

grid.arrange(RLC.pres, RLG.pres, layout_matrix = rbind(c(1,1),c(2,2)))
