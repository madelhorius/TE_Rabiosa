# ---------------------------------------------------------------------------------#
# Consensus sequence from clustalw2 alignment file --------------------------------#
# Marius Hodel --------------------------------------------------------------------#
# 2019-04-24 ----------------------------------------------------------------------#
#----------------------------------------------------------------------------------#

args <- commandArgs(trailingOnly = TRUE)
if(length(args)<1) stop("didn't receive 2 arguments")
alnfile <- args[1]
outfile <- args[2]
if(!file.exists(alnfile)) stop("1st argument isn't an existing file")


# setwd("P:/mahodel")

# Load/install required packages
.libPaths("/home/mahodel/bin/R")
if(! require("dplyr")) {
  install.packages("dplyr", repos="https://stat.ethz.ch/CRAN/", lib = "/home/mahodel/bin/R")
  require("dplyr")
}
if(! require("tidyr")) {
  install.packages("tidyr", repos="https://stat.ethz.ch/CRAN/", lib = "/home/mahodel/bin/R")
  require("tidyr")
}
# Read in aln file

aln <- read.table(alnfile, skip = 2)
colnames(aln) <- c("TE", "aln")

# Rearrange data
aln <- aln %>% 
  group_by(TE)  %>%
  summarise(seq = paste(aln, collapse = ""))

aln <- separate(aln, seq, as.character(c(1:nchar(aln$seq[1]))), sep = c(1:(nchar(aln$seq[1])-1)))
aln <- as.data.frame(aln)
cons <- vector()
for (i in 2:ncol(aln)) {
  freq <- as.data.frame(table(aln[,i]))
  freq$Freq <- freq$Freq/sum(freq$Freq)
  freq <- freq[order(freq$Freq, decreasing = T),]
  base <- as.character(freq$Var1[1])
  cons[i-1] <- base
}

cons.seq <- as.data.frame(paste(cons, collapse = ""))
name <- paste(">", outfile, sep = "")
colnames(cons.seq) <- name

f <- file(outfile, open="wb")
write.table(cons.seq, file = f, quote = F, eol = "\n", row.names = F, col.names = T, sep="\t")
close(f)

