# ---------------------------------------------------------------------------------#
# Parse RM-output -----------------------------------------------------------------#
# To find LTR-copies (FL, truncated, Solo) ----------------------------------------#
# Marius Hodel --------------------------------------------------------------------#
# 2019-05-06 ----------------------------------------------------------------------#
#----------------------------------------------------------------------------------#


# Command line arguments
args <- commandArgs(trailingOnly = TRUE)
if(length(args)<3) stop("didn't receive 3 arguments")
rm.out.file <- args[1]
parsed.gff.file <- args[2]
sizeltr <- args[3]
sizeltr <- as.numeric(sizeltr)
if(!file.exists(rm.out.file)) stop("1st argument isn't an existing file")

# setwd("P:/mahodel")

# Read Repeatmasker output
repma <- read.table(rm.out.file, skip = 3, fill = T, header = F, 
                    colClasses = c("integer", rep("numeric",3), "factor", rep("integer", 2), rep("character", 7), "integer"))
colnames(repma) <- c("score", "div", "del", "ins",  "sequence", "seq.begin", "seq.end", "seq.left", "strand", "TE", "class/family", "rep.begin", 
                     "rep.end", "rep.left", "ID")
tmp <- repma[repma$strand == "C",]$rep.begin
repma[repma$strand == "C",]$rep.begin <-  repma[repma$strand == "C",]$rep.end
repma[repma$strand == "C",]$rep.end <-  repma[repma$strand == "C",]$rep.left
repma[repma$strand == "C",]$rep.left <- tmp 
rm(tmp)
repma$rep.begin <- as.integer(repma$rep.begin)

repma$rep.end <- as.integer(repma$rep.end)
repma$strand <- as.factor(repma$strand)

# Hits have to be at least 80 bp long to get considered

repma <- repma[(repma$seq.end - repma$seq.begin) > 79,]

# Extract family
repma$fam <- gsub(pattern = "_int-int", replacement = "", x = repma$TE)
repma$fam <- gsub(pattern = "_LTR", replacement = "", x = repma$fam)


# Prepare data.frames & vectors
FL <- as.data.frame(matrix(nrow = 0, ncol = 9))
colnames(FL) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")

Trunc <- as.data.frame(matrix(nrow = 0, ncol = 9))
colnames(Trunc) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")

Solo <- as.data.frame(matrix(nrow = 0, ncol = 9))
colnames(Solo) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")

nrfl <- 1 # Number of FL-element
nrRfl <- 1 # Number of rows in gff of full length
nrtr <- 1 # Number of truncated elements
nrRtr <- 1 # Number of rows in gff of truncated elements
nrsol <- 1 # Number of Solo-element

for (i in 1:max(repma$ID)){
  TE <- repma[repma$ID == i, ]
  if (nrow(TE) == 0) next
  # Does the element contain an internal sequence?
  int <- TE[grepl("int-int", TE$TE),]
  if (nrow(int) > 0) {
    ltr <- TE[grepl("LTR", TE$TE),]
    if (nrow (ltr) > 0) {
      whichLTR <- grepl("LTR", TE$TE)
      # TE starts and ends with an LTR => FL element
      if (whichLTR[length(whichLTR)] & whichLTR[1]) {
        # Check if it is a complex insertion (containing additional LTR-sequences within the internal sequence)
        intstart <- min(TE[grepl("int", TE$TE),6])
        intend <- max(TE[grepl("int", TE$TE),7])
        ltr.within <- TE[grepl("LTR", TE$TE) & TE$seq.begin > intstart & TE$seq.end < intend,]
        if(nrow(ltr.within) > 1){ # complex
          feature <- "FL-LTR-RT-complex"
        } else {
          feature <- "FL-LTR-RT"
        }
        FL[nrRfl,] <- NA
        FL$seqname[nrRfl] <- as.character(TE$sequence[1])
        FL$source[nrRfl] <- "RepeatMasker_R"
        FL$feature[nrRfl] <- feature
        FL$start[nrRfl] <- min(TE$seq.begin)
        FL$end[nrRfl] <- max(TE$seq.end)
        FL$attribute[nrRfl] <- paste("ID=", TE$fam[1], ".", nrfl, sep = "")
        FL$score[nrRfl] <- sum(TE$score)
        if (TE$strand[1] == "+") {
          FL$strand[nrRfl] <- "+" 
        } else {
          FL$strand[nrRfl] <- "-"
        }
        # counters
        nrRfl <- nrRfl + 1
        
        # Add parts if the gaps between the hits are bigger than 79 bp (to allow nesting laTEr)
        nrpart <- 1
        partStart <- min(TE$seq.begin)
        for (k in 1:(nrow(TE)-1)) {
          inslen <- TE$seq.begin[k+1] - TE$seq.end[k]
          if (inslen > 79) {
            FL[nrRfl,] <- NA
            FL$seqname[nrRfl] <- as.character(TE$sequence[1])
            FL$source[nrRfl] <- "RepeatMasker_R"
            FL$feature[nrRfl] <- "part"
            FL$start[nrRfl] <- partStart
            FL$end[nrRfl] <- TE$seq.end[k]
            FL$attribute[nrRfl] <- paste("ID=", TE$fam[1], ".", nrfl, ":part:", nrpart, ";Parent=", TE$fam[1], ".", nrfl, sep = "")
            FL$score[nrRfl] <- sum(TE$score)
            # Start of next part
            partStart <- TE$seq.begin[k+1]
            # counters
            nrRfl <- nrRfl + 1
            nrpart <- nrpart + 1
          }
          # Last part
          if (k == (nrow(TE)-1)){
            FL[nrRfl,] <- NA
            FL$seqname[nrRfl] <- as.character(TE$sequence[1])
            FL$source[nrRfl] <- "RepeatMasker_R"
            FL$feature[nrRfl] <- "part"
            FL$start[nrRfl] <- partStart
            FL$end[nrRfl] <- TE$seq.end[k+1]
            FL$attribute[nrRfl] <- paste("ID=", TE$fam[1], ".", nrfl, ":part:", nrpart, ";Parent=", TE$fam[1], ".", nrfl, sep = "")
            FL$score[nrRfl] <- sum(TE$score)
            
            nrRfl <- nrRfl + 1
          }
        }
        nrfl <- nrfl + 1
      } else { # truncated TE (only one LTR and some internal sequence) => Add to gff file
        
        # Check if it is a complex insertion (containing mutliple LTRs within the internal region)
        intstart <- min(TE[grepl("int", TE$TE),6])
        intend <- max(TE[grepl("int", TE$TE),7])
        ltr.within <- TE[grepl("LTR", TE$TE) & TE$seq.begin > intstart & TE$seq.end < intend,]
        if (nrow(ltr.within) > 0) {
          feature <- "Trunc-LTR-RT-complex"
        } else {
          feature <- "Trunc-LTR-RT"
        }
        
        Trunc[nrRtr,] <- NA
        Trunc$seqname[nrRtr] <- as.character(TE$sequence[1])
        Trunc$source[nrRtr] <- "RepeatMasker_R"
        Trunc$feature[nrRtr] <- feature
        Trunc$start[nrRtr] <- min(TE$seq.begin)
        Trunc$end[nrRtr] <- max(TE$seq.end)
        Trunc$attribute[nrRtr] <- paste("ID=Trunc.", TE$fam[1], ".", nrtr, sep = "")
        Trunc$score[nrRtr] <- sum(TE$score)
        if (TE$strand[1] == "+") {
          Trunc$strand[nrRtr] <- "+" 
        } else {
          Trunc$strand[nrRtr] <- "-"
        }
        # counters
        nrRtr <- nrRtr + 1
        
        nrpart <- 1
        partStart <- min(TE$seq.begin)
        for (k in 1:(nrow(TE)-1)) {
          inslen <- TE$seq.begin[k+1] - TE$seq.end[k]
          if (inslen > 79) {
            Trunc[nrRtr,] <- NA
            Trunc$seqname[nrRtr] <- as.character(TE$sequence[1])
            Trunc$source[nrRtr] <- "RepeatMasker_R"
            Trunc$feature[nrRtr] <- "part"
            Trunc$start[nrRtr] <- partStart
            Trunc$end[nrRtr] <- TE$seq.end[k]
            Trunc$attribute[nrRtr] <- paste("ID=Trunc.", TE$fam[1], ".", nrtr, ":part:", nrpart, ";Parent=Trunc.", TE$fam[1], ".", nrtr, sep = "")
            Trunc$score[nrRtr] <- sum(TE$score)
            # Start of next part
            partStart <- TE$seq.begin[k+1]
            # counters
            nrRtr <- nrRtr + 1
            nrpart <- nrpart + 1
          }
          # Last part
          if (k == (nrow(TE)-1)){
            Trunc[nrRtr,] <- NA
            Trunc$seqname[nrRtr] <- as.character(TE$sequence[1])
            Trunc$source[nrRtr] <- "RepeatMasker_R"
            Trunc$feature[nrRtr] <- "part"
            Trunc$start[nrRtr] <- partStart
            Trunc$end[nrRtr] <- TE$seq.end[k+1]
            Trunc$attribute[nrRtr] <- paste("ID=Trunc.", TE$fam[1], ".", nrtr, ":part:", nrpart, ";Parent=Trunc.", TE$fam[1], ".", nrtr, sep = "")
            Trunc$score[nrRtr] <- sum(TE$score)
            
            nrRtr <- nrRtr + 1
          }
        }
        nrtr <- nrtr + 1
      }
    } else { # Only internal sequence => skip
      next
    }
  } else { # Add Solo LTR to gff file
    
    # Check length of Solo LTR-candidate so that small hits are not considered as Solo-LTR
    TE$len <- TE$seq.end - TE$seq.begin
    len <- sum(TE$len)
    if (len < (sizeltr*0.8)) next
    
    Solo[nrsol,] <- NA
    Solo$seqname[nrsol] <- as.character(TE$sequence[1])
    Solo$source[nrsol] <- "RepeatMasker_R"
    Solo$feature[nrsol] <- "Solo-LTR"
    Solo$start[nrsol] <- min(TE$seq.begin)
    Solo$end[nrsol] <- max(TE$seq.end)
    Solo$attribute[nrsol] <- paste("ID=Solo.", TE$fam[1], ".", nrsol, sep = "")
    Solo$score[nrsol] <- sum(TE$score)
    if (TE$strand[1] == "+") {
      Solo$strand[nrsol] <- "+" 
    } else {
      Solo$strand[nrsol] <- "-"
    }
    nrsol <- nrsol + 1
  }
}

# Check of truncated elements next to each other can be matched to a full length element with a big insertion (what Repeatmasker could have missed)
tru <- Trunc[Trunc$feature == "Trunc-LTR-RT" | Trunc$feature == "Trunc-LTR-RT-complex",]
if(nrow(tru) > 1) {
  for (i in 1:(nrow(tru) - 1)){
    # Same scaffold?
    if (tru$seqname[i]!= tru$seqname[i+1]) next 
    # Same strand?
    strand1 <- tru$strand[i]
    strand2 <- tru$strand[i+1]
    if (strand1 != strand2) next 
    # Right order
    if (!(grepl("LTR", repma[repma$seq.begin == tru$start[i] & repma$sequence == tru$seqname[i],]$TE))) next
    if (!(grepl("int", repma[repma$seq.begin == tru$start[i+1] & repma$sequence == tru$seqname[i],]$TE))) next
    # Distance
    if ((tru$start[i+1] - tru$end[i]) > 30000) next
    # Repeat
    tru1 <- repma[repma$seq.begin >= tru$start[i] & repma$seq.end <= tru$end[i] & repma$sequence == tru$seqname[i],]
    tru2 <- repma[repma$seq.begin >= tru$start[i+1] & repma$seq.end <= tru$end[i+1] & repma$sequence == tru$seqname[i],]
    if (strand1 == "-") {
      repdiff <- tail(tru1$rep.end,1) - tru2$rep.begin[1]
      if (repdiff < -150) next
    } else {
      repdiff <- tru2$rep.begin[1] - tail(tru1$rep.end,1)
      if (repdiff < -150) next
    }
    # Passed all tests => put together to a full length element
    # Check if complex
    if (tru$feature[i] == "Trunc-LTR-RT-complex" | tru$feature[i+1] == "Trunc-LTR-RT-complex"){
      feature <- "FL-LTR-RT-complex"
    } else {
      feature <- "FL-LTR-RT"
    }
    attr1 <- tru[i,9]
    attr2 <- tru[i+1,9]
    newFL <- Trunc[grepl(paste(attr1, ":part|", attr2, ":part", sep = ""), Trunc$attribute) & Trunc$feature == "part",]
    score <- sum(Trunc[Trunc$attribute == attr1 | Trunc$attribute == attr2,]$score)
    FL[nrRfl,] <- NA
    FL$seqname[nrRfl] <- as.character(newFL$seqname[1])
    FL$source[nrRfl] <- "RepeatMasker_R"
    FL$feature[nrRfl] <- feature
    FL$start[nrRfl] <- min(newFL$start)
    FL$end[nrRfl] <- max(newFL$end)
    FL$score[nrRfl] <- score
    FL$attribute[nrRfl] <- paste("ID=", TE$fam[1], ".", nrfl, sep = "")
    FL$strand[nrRfl] <- strand1
    nrRfl <- nrRfl + 1
    nrpart <- 1
    for (row in 1:nrow(newFL)){
      FL[nrRfl,] <- newFL[row,]
      FL[nrRfl,]$score <- score
      FL$attribute[nrRfl] <- paste("ID=", TE$fam[1], ".", nrfl, ":part:", nrpart, ";Parent=", TE$fam[1], ".", nrfl, sep = "")
      nrpart <- nrpart + 1
      nrRfl <- nrRfl + 1
    }
    nrfl <- nrfl + 1
    # Remove these truncated elements from the gff file
    ID1 <- tru$attribute[i]
    ID2 <- tru$attribute[i+1]
    Trunc <- Trunc[Trunc$attribute != ID2 & Trunc$attribute != ID1,]
    ID1 <- paste(ID1, ":part", sep = "")
    ID2 <- paste(ID2, ":part", sep = "")
    Trunc <- Trunc[!(grepl(paste(ID1, ID2, sep = "|"), Trunc$attribute)),]
  }
}

# Check if Solo LTR could be part of a truncated element => and therefore a FL-element
toRM <- vector()
if (nrow(Solo) > 0 & nrow(Trunc) > 0) {
  for (i in 1:nrow(Solo)){
    start <- Solo[i,]$start
    end <- Solo[i,]$end
    truLeft <- Trunc[(Trunc$feature == "Trunc-LTR-RT" | Trunc$feature == "Trunc-LTR-RT-complex") & Trunc$end > (start - 5000) & Trunc$end < start & Trunc$seqname == Solo$seqname[i],]
    truRight <- Trunc[(Trunc$feature == "Trunc-LTR-RT" | Trunc$feature == "Trunc-LTR-RT-complex") & Trunc$start < (end + 5000) & Trunc$start > end & Trunc$seqname == Solo$seqname[i],]
    if(nrow(truLeft) == 0 & nrow(truRight) == 0) next
    # Check left candidate
    if(nrow(truLeft) != 0) {
      truLeft <- truLeft[order(truLeft$end, decreasing = T),]
      if(grepl("int-int", repma[repma$seq.end == truLeft$end[1] & repma$sequence == truLeft$seqname[1],]$TE)){
        if(truLeft$strand[1] == Solo[i,]$strand){
          left <- T
        } else {
          left <- F
        }
      } else {
        left <- F
      }
    } else {
      left <- F
    }
    # Check right candidate
    if(nrow(truRight) != 0) {
      truRight <- truRight[order(truRight$start),]
      if(grepl("int-int", repma[repma$seq.begin == truRight$start[1] & repma$sequence == truRight$seqname[1],]$TE)) {
        if(truRight$strand[1] == Solo[i,]$strand){
          right <- T
        } else {
          right <- F
        }
      } else {
        right <- F
      } 
    } else {
      right <- F
    }
    # 2 possible candidates: Take the one closer
    if(right & left) {
      diffleft <- start - truLeft$end[1]
      diffright <- truRight$start[1] - end
      if (diffleft > diffright) left <- F else right <- F
    }
    # No candidates anymore
    if(!(right) & !(left)) next
    
    # Passed all tests =>  Put Solo-LTR and Truncated element together to a full length element
    if (left){
      attr1<- truLeft$attribute[1]
      newFL <- Trunc[Trunc$attribute == attr1,]
      # Check if complex
      if (newFL$feature[1] == "Trunc-LTR-RT-complex") {
        feature <- "FL-LTR-RT-complex"
      } else {
        feature <- "FL-LTR-RT"
      }
      attr1 <- paste(attr1, "part", sep =":")
      newFL <- rbind(newFL, Trunc[grepl(attr1,Trunc$attribute),], Solo[i,])
      # Check if LTR on both sides of new TE
      logFL <- grepl("LTR", repma[repma$seq.begin == min(newFL$start) & repma$sequence == newFL$seqname[1],]$TE)
    } else {
      attr1<- truRight$attribute[1]
      newFL <- Trunc[Trunc$attribute == attr1,]
      # Check if complex
      if (newFL$feature[1] == "Trunc-LTR-RT-complex") {
        feature <- "FL-LTR-RT-complex"
      } else {
        feature <- "FL-LTR-RT"
      }
      attr1 <- paste(attr1, "part", sep =":")
      newFL <- rbind(Solo[i,], newFL, Trunc[grepl(attr1,Trunc$attribute),])
      # Check if LTR on both sides of new TE
      logFL <- grepl("LTR", repma[repma$seq.end == max(newFL$end) & repma$sequence == newFL$seqname[1],]$TE)
    }
    score <- sum(newFL[newFL$feature != "part",]$score)
    # Add new FL-TE
    if (logFL) {
      FL[nrRfl,] <- NA
      FL$seqname[nrRfl] <- as.character(newFL$seqname[1])
      FL$source[nrRfl] <- "RepeatMasker_R"
      FL$feature[nrRfl] <- feature
      FL$start[nrRfl] <- min(newFL$start)
      FL$end[nrRfl] <- max(newFL$end)
      FL$score[nrRfl] <- score
      FL$attribute[nrRfl] <- paste("ID=", TE$fam[1], ".", nrfl, sep = "")
      FL$strand[nrRfl] <- Solo[i,]$strand
      nrRfl <- nrRfl + 1
      nrpart <- 1
      newFL <- newFL[newFL$feature != "Trunc-LTR-RT" & newFL$feature != "Trunc-LTR-RT-complex",]
      for (row in 1:nrow(newFL)){
        FL[nrRfl,] <- newFL[row,]
        FL$score[nrRfl] <- score
        FL$feature[nrRfl] <- "part"
        FL$attribute[nrRfl] <- paste("ID=", TE$fam[1], ".", nrfl, ":part:", nrpart, ";Parent=", TE$fam[1], ".", nrfl, sep = "")
        nrpart <- nrpart + 1
        nrRfl <- nrRfl + 1
      }
      nrfl <- nrfl + 1
      # Remove these truncated elements and Solo-LTRs from the gff file
      if (left) {
        ID <- Trunc[grepl(truLeft$attribute[1],Trunc$attribute) & (Trunc$feature == "Trunc-LTR-RT" | Trunc$feature == "Trunc-LTR-RT-complex"),]$attribute[1]
      } else {
        ID <- Trunc[grepl(truRight$attribute[1],Trunc$attribute) & (Trunc$feature == "Trunc-LTR-RT" | Trunc$feature == "Trunc-LTR-RT-complex"),]$attribute[1]
      }
      Trunc <- Trunc[!(grepl(paste(ID, ":", sep =""), Trunc$attribute)),]
      Trunc <- Trunc[Trunc$attribute != ID,]
    } else {
      # Add new part to Trunc-LTR-RT
      if (left) {
        ID <- Trunc[grepl(truLeft$attribute[1],Trunc$attribute) & (Trunc$feature == "Trunc-LTR-RT" | Trunc$feature == "Trunc-LTR-RT-complex"),]$attribute[1]
        Trunc[Trunc$attribute == ID,]$end <- max(newFL$end)
        Trunc[Trunc$attribute == ID,]$score <- score
        Trunc[(grepl(paste(ID, ":", sep =""), Trunc$attribute)),]$score <- score
        newPart <- tail(newFL,1)
        newPart$feature <- "part"
        newPart$score <- score
        newPart$attribute <- paste(head(newFL$attribute,1), ":part:", (nrow(newFL)-1), ";Parent=", gsub("ID=", "", head(newFL$attribute,1)), sep = "")
        Trunc <- rbind(Trunc, newPart)
      } else {
        ID <- Trunc[grepl(truRight$attribute[1],Trunc$attribute) & (Trunc$feature == "Trunc-LTR-RT" | Trunc$feature == "Trunc-LTR-RT-complex"),]$attribute[1]
        Trunc[Trunc$attribute == ID,]$start <- min(newFL$end)
        Trunc[Trunc$attribute == ID,]$score <- score
        Trunc[(grepl(paste(ID, ":", sep =""), Trunc$attribute)),]$score <- score
        newPart <- head(newFL,2)
        newPart$feature <- "part"
        newPart$score <- score
        newPart$attribute <- paste(newFL$attribute[2], ":part:", (nrow(newFL)-1), ";Parent=", gsub("ID=", "", head(newFL$attribute,1)), sep = "")
        Trunc <- rbind(Trunc, newPart)
      }
      
    }
    # Remember which Solo-LTRs to remove
    toRM <- c(toRM, i)
  } 
}


# Remove Solo-LTRs
if (length(toRM) > 0) Solo <- Solo[-toRM,]

all <- rbind(FL, Trunc, Solo)
all[is.na(all)] <- "."

f <- file(parsed.gff.file, open="wb")
write.table(all, file = f, quote = F, eol = "\n", row.names = F, col.names = F, sep="\t")
close(f)
