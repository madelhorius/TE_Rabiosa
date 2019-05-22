#------------------------------------#
#-----Find nested elements-----------#
#-----Marius Hodel-------------------#
#-----2019-02-18---------------------#
#------------------------------------#

# Set wd
setwd("P:/mahodel")
library(ggplot2)

# Read data
gff <- read.table("FL_PARSED_FL_LTR_RT.pass.code.gff", sep = "\t")

# Find nested elements

# create data.frame
nested <- gff[0,]

# loops to find nested elements
for (i in 1:nlevels(gff[,1])){
  tmp <- gff[gff[,1] == levels(gff[,1])[i],]
  if (nrow(tmp) > 1){
    for (j in 1:(nrow(tmp)-1)){
      diff <- tmp[j+1,4] - tmp[j,5]
      if (diff <= 0){
        if (j < (nrow(tmp)-1)){
          diff2 <- tmp[j+2,4] - tmp[j,5]
          if (diff2 <= 0){
            if (j < (nrow(tmp)-2)){
              diff3 <- tmp[j+3,4] - tmp[j,5]
              if (diff3 <=0) {
                message("3")
                message(j)
                nested <- rbind(nested, tmp[c(j,j+1,j+2,j+3),])
              } else {
                nested <- rbind(nested, tmp[c(j,j+1,j+2),])
              }
            } else {
              nested <- rbind(nested, tmp[c(j,j+1,j+2),])
            }
          } else {
            nested <- rbind(nested, tmp[c(j,j+1),])
          }
        } else {
          nested <- rbind(nested, tmp[c(j,j+1),])
        }
      }
    } 
  }
}

nested <- unique(nested)

# Get host TEs
d <- 0
nested$position <- NA
nested$nrIns <- 0
nested[,1] <- as.factor(as.character(nested[,1]))
nested2 <- nested[1,]
nested2 <- nested2[-1,]

for (j in 1:nlevels(nested[,1])){
  d <- 0
  tmp <- nested[nested[,1] == levels(nested[,1])[j],]
  for (i in 1:(nrow(tmp)-1)){
    if (i < d) next
    diff <- -1
    nrIns <- 1
    tmp$position[i] <- "host"
    while (diff < 0 & (i+nrIns <= nrow(tmp))){
      diff <- tmp[i+nrIns,4] - tmp [i,5]
      if (diff < 0){
        tmp$position[i+nrIns] <- "insert"
        tmp$nrIns[i] <- nrIns
        nrIns <- nrIns+1
        d <- i+nrIns
      }
    }
  }
  nested2 <- rbind(nested2, tmp)
}

nested <- nested2
rm(nested2)

nested <- within(nested, V9 <- data.frame(do.call('rbind', strsplit(as.character(V9), ';', fixed=TRUE))))
nested$V9$X4 <- gsub("ltr_identity=","", nested$V9$X4)
nested$V9$X4 <- as.numeric(nested$V9$X4)

# Calcualte insertion time (according to supplemental material of Ou et al. 2018)
nested$insertiontime <- 1 - nested$V9$X4 # d = 1 - similarity
nested$insertiontime <- -3/4*log(1-nested$insertiontime*4/3) # K = -3/4*log(1-d*4/3)
nested$insertiontime <- nested$insertiontime/(2*(1.3*10^-8)) # T = K/(2*(1.3*10^-8))
nested$insertiontime <- round(nested$insertiontime/10^6,2) # in myrs

# Calculate difference in insertion times between host and insertion
nested$diff <- 0
nested2 <- nested[0,]
for (j in 1:nlevels(nested[,1])){
  tmp <- nested[nested[,1] == levels(nested[,1])[j],]
  for (i in 1:(nrow(tmp))){
    if (tmp$nrIns[i] == 0) next
    for (k in 1:tmp$nrIns[i]){
      tmp$diff[i+k] <- tmp$insertiontime[i]-tmp$insertiontime[i+k]
    }
  }
  nested2 <- rbind(nested2, tmp)
}

nested <- nested2
rm(nested2)

ggplot(nested[nested$nrIns == 0,], aes(diff)) +
  geom_histogram(binwidth = 0.05, col = "black", boundary = 0, closed = "left") +
  xlab("Time difference between host and inserted element [myrs]") +
  theme_bw() +
  geom_vline(xintercept = 0, col = "red")

nrow(nested[nested$nrIns == 0 & nested$diff >= 0,])


# Define cases
# Case 1: Inserted element is completely within the host; at least 100 bp on each side:
# --------------------------------------------------------------------------------------
# ....> 100 bp ....---------------------------------------------........>100 bp.........

case1 <- nested[0,]
for (i in 1:nrow(nested)){
  if (nested$position[i] == "host" & nested$nrIns[i] == 1){
    diffLeft <- nested[i+1,4] - nested[i,4]
    diffRight <- nested[i,5] - nested[i+1,5]
    if (diffLeft > 99 & diffRight > 99){
      case1 <- rbind(case1, nested[c(i,(i+1)),])
    }
  } else {
    next
  }
}

ggplot(case1[case1$nrIns == 0,], aes(diff)) +
  geom_histogram(binwidth = 0.05, col = "black", boundary = 0, closed = "left") +
  xlab("Time difference between host and inserted element [myrs]") +
  theme_bw() +
  geom_vline(xintercept = 0, col = "red")

nrow(case1[case1$nrIns == 0 & case1$diff >= 0,])/nrow(case1[case1$nrIns == 0,])

# Export cases to split fasta files on server

# Export all nested elements
nestedTE <- nested$V9$X1 
nestedTE <- gsub("ID=LTR_retrotransposon","LTR_RT", nestedTE)
nestedTE <- paste(nestedTE, "_", sep = "")

# Complicated command, so that LF and not CRLF and end of each line! => otherwhise it cannot be used for grep on Linux
f <- file("NestedTEs.txt", open="wb")
write.table(nestedTE, file = f, quote = F, eol = "\n", row.names = F, col.names = F)
close(f)

# Export simple cases
simplyNestedTE <- case1$V9$X1
simplyNestedTE <- gsub("ID=LTR_retrotransposon","LTR_RT", simplyNestedTE)
simplyNestedTE <- paste(simplyNestedTE, "_", sep = "")
f <- file("simplyNestedTEs.txt", open="wb")
write.table(simplyNestedTE, file = f, quote = F, eol = "\n", row.names = F, col.names = F)
close(f)

# Create gff file with split
# --------------------------

split <- read.table("PARSED_FL_LTR_RT.gff", sep = "\t")
split$V3 <- as.character(split$V3)
orig <- split$V3
split$V9 <- as.character(split$V9)
split[seq(4, nrow(split), 5),3] <- "split1"
split[seq(5, nrow(split), 5),3] <- "split2"


for (i in 1:nrow(split)){
  if (split$V3[i] != "split1" & split$V3[i] != "split2") next
  if (split$V3[i] == "split1"){
    split[i,5] <- split[i,4]-1
    split[i,4] <- split[i-3,4]
    split[i,9] <- paste(split[i-2,9], split[i,3], sep = ";")
  } else {
    split[i,4] <- split[i,5]+1
    split[i,5] <- split[i-4,5]
    split[i,9] <- split[i-2,9]
    split[i,9] <- paste(split[i-2,9], split[i,3], sep = ";")
  }
}

split$V3 <- orig

# Complicated command, so that LF and not CRLF and end of each line! => otherwhise it cannot be used for grep on Linux
f <- file("PARSED_FL_LTR_RT.gff", open="wb")
write.table(split, file = f, quote = F, eol = "\n", row.names = F, col.names = F, sep="\t")
close(f)

