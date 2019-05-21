# ---------------------------------------------------------------------------------#
# Protein domains LTR-RT ----------------------------------------------------------#
# Get domains per LTR-RT according to ltrdigest gff output and write output table--#
# Marius Hodel --------------------------------------------------------------------#
# 2019-03-07 ----------------------------------------------------------------------#
#----------------------------------------------------------------------------------#

args <- commandArgs(trailingOnly = TRUE)
if(length(args)<1) stop("didn't receive 2 arguments")
gfffile <- args[1]
indexfile <- args[2]
if(!file.exists(gfffile)) stop("1st argument isn't an existing file")
if(!file.exists(indexfile)) stop("2nd argument isn't an existing file")

# Load/install required packages
.libPaths("/home/mahodel/bin/R")
if(! require("ape")) {
  install.packages("ape", repos="https://stat.ethz.ch/CRAN/", lib = "/home/mahodel/bin/R")
  require("ape")
}
if(! require("data.table")) {
  install.packages("data.table", repos="https://stat.ethz.ch/CRAN/", lib = "/home/mahodel/bin/R")
  require("data.table")
}
if(! require("gtools")) {
  install.packages("gtools", repos="https://stat.ethz.ch/CRAN/", lib = "/home/mahodel/bin/R")
  require("gtools")
}

# setwd("P:/mahodel")

dom <- read.gff(file = gfffile)
dom <- within(dom, attributes <- data.frame(do.call('rbind', strsplit(as.character(attributes), ';', fixed=F))))

dom$attributes$X1 <- gsub("ID=","", dom$attributes$X1)
dom$attributes$X1 <- gsub("Parent=", "", dom$attributes$X1)
dom$attributes$X1 <- as.factor(dom$attributes$X1)
dom$attributes$X1 <- factor(dom$attributes$X1, levels=mixedsort(levels(dom$attributes$X1)))

dom <- dom[!(dom$attributes$X3 %like% "CHR"),]
dom <- dom[!(dom$attributes$X3 %like% "GAGCOAT"),]
dom <- dom[!(dom$attributes$X3 %like% "ENV"),]

prot.domains <- as.data.frame(matrix(nrow = nlevels(dom$attributes$X1), ncol = 3))

for (i in 1:nlevels(dom$attributes$X1)){
  tmp <- dom[dom$attributes$X1 == levels(dom$attributes$X1)[i],]
  
  # If TE is host of inserted elment => Do not include protein domains of the nested elements
  if("split=split1" %in% tmp$attributes$X2) {
    split1.end <- tmp[tmp$attributes$X2 %like% "split1",]$end
    split2.start <- tmp[tmp$attributes$X2 %like% "split2",]$start
    tmp <- tmp[!(tmp$end < split2.start & tmp$end > split1.end),]
    tmp <- tmp[!(tmp$start < split2.start & tmp$start > split1.end),]
  }
  
  tmp <- tmp[tmp$attributes$X2 %like% "reading_frame",]
  X <- rle(gsub("_.*", "", tmp$attributes$X3))
  Y <- cumsum(c(1, X$lengths[-length(X$lengths)]))
  tmp <- tmp[Y,]
  tmp$attributes$X3 <- gsub("name=", "", tmp$attributes$X3)
  tmp$attributes$X3 <- gsub("_.*", "", tmp$attributes$X3)
  tmp$attributes$X3 <- gsub("Nase", "", tmp$attributes$X3)
  tmp <- do.call(data.frame, tmp)
  lidom <- tmp[,c(4,11)]
  lidom <- lidom[order(lidom[,1]),]
  # Get normal order if on - strand
  if (length(na.omit(tmp$strand[2])) > 0){
    if (na.omit(tmp$strand[2]) == "-") lidom <- lidom[order(-lidom[,1]),]
  }
  
  # Save protein domains and their order
  prot.domains[i,1] <- as.character(levels(dom$attributes$X1)[i])
  prot.domains[i,2] <- nrow(na.omit(lidom))
  prot.domains[i,3] <- paste(na.omit(lidom)[,2], collapse = "|")
  if(prot.domains[i,3] == "") prot.domains[i,3] <- NA
  
  if (i %% 100  == 0) message(paste(i, "Elements processed"))
}

# Get back original numbering of LTR-RT
index <- read.table(indexfile, sep = "\t")
prot.domains$ID <- index[,2]

# Add Gypsy and Copia according to order of domains
prot.domains$superfam <- NA

prot.domains$superfam <- ifelse(!is.na(prot.domains[,3]) 
                                & (prot.domains[,3] == "GAG|AP|RT|RH|INT" | prot.domains[,3] == "RT|RH|INT"
                                  | prot.domains[,3] == "RH|INT" | prot.domains[,3] == "RT|INT" | prot.domains[,3] == "GAG|RT|RH|INT"
                                  | prot.domains[,3] == "GAG|RH|INT" | prot.domains[,3] == "GAG|RT|INT" | prot.domains[,3] == "AP|RT|RH|INT"
                                  | prot.domains[,3] == "AP|RH|INT" | prot.domains[,3] == "AP|RT|INT" | prot.domains[,3] == "GAG|AP|RH|INT"
                                  | prot.domains[,3] == "GAG|AP|RT|INT"), "RLG", NA)

prot.domains$superfam <- ifelse(!is.na(prot.domains[,3]) 
                                & (prot.domains[,3] == "GAG|AP|INT|RT|RH" | prot.domains[,3] == "INT|RT|RH"| 
                                   prot.domains[,3] == "INT|RT"| prot.domains[,3] == "INT|RH"|
                                   prot.domains[,3] == "GAG|INT|RT|RH"| prot.domains[,3] == "GAG|INT|RT"| 
                                   prot.domains[,3] == "GAG|INT|RH"| prot.domains[,3] == "AT|INT|RT|RH"|
                                   prot.domains[,3] == "AT|INT|RT"| prot.domains[,3] == "AT|INT|RH"|
                                   prot.domains[,3] == "GAG|AT|INT|RT"| prot.domains[,3] == "GAG|AT|INT|RH"), "RLC", prot.domains$superfam)


f <- file("LTR_Protein_domains.txt", open="wb")
write.table(prot.domains, file = f, quote = F, eol = "\n", row.names = F, col.names = F, sep="\t")
close(f)

message("ALL DONE! LIFE SEEMS TO BE A BISCUIT FOR YOU!")



