args <- commandArgs(trailingOnly = TRUE)
if(length(args)<1) stop("didn't receive at least 1 argument")
gfffile <- args[1]
if(!file.exists(gfffile)) stop("1st argument isn't an existing file")


#setwd("P:/mahodel")
#gfffile <- "0.1_PARSED_FL_LTR_RT.pass.gff"

gff <- read.table(gfffile)

gff <- within(gff, V9 <- data.frame(do.call('rbind', strsplit(as.character(V9), ';', fixed=TRUE))))
gff[,9] <- gff$V9$X1 
gff[,9] <- gsub("ID=","", gff[,9])
gff[,9] <- gsub("Parent=", "", gff[,9])
gff[,9] <- as.factor(gff[,9])

longLTR <- vector()
for (i in 1:nlevels(gff[,9])){
  tmp <- gff[gff[,9] == levels(gff[,9])[i],]
  if (nrow(tmp) == 3){
    len <- tmp[1,5] - tmp[1,4]
    ltrlen <- tmp[2,5] - tmp[2,4]
    ratio <- ltrlen/len
  } else if (nrow(tmp) == 5) {
    split1len <- tmp[4,5] - tmp[4,4]
    split2len <- tmp[5,5] - tmp[5,4]
    len <- split1len + split2len
    ltrlen <- tmp[2,5] - tmp[2,4]
    ratio <- ltrlen/len
  } else {
    stop(paste(tmp[1,9], "does not contain 3 or 5 lines in gff file!", sep = " "))
  }
  if (ratio >= 0.4) {
    longLTR <- c(longLTR, as.character(tmp[1,9]))
  }
}

f <- file("LongLTR.lst", open="wb")
write.table(longLTR, file = f, quote = F, eol = "\n", row.names = F, col.names = F, sep="\t")
close(f)
