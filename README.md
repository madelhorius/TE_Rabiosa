# TE_Rabiosa

This repository includes all important scripts that were written for the Master`s thesis "Identification and characterization of long terminal repeat retrotransposons in Lolium multiflorum"

- `LTR-RT/validation/selfblast.sh` checks for repeats within LTR-RTs and identifies TEs with long LTRs.
- `LTR-RT/validation/domains.R` extracts protein domain regions from the ltrdigest gff-output.
- `LTR-RT/validation/nestedElements.R` finds nested TEs in a gff file, simple cases are exported, the host is split.
- `LTR-RT/validation/longLTRcheck.R` identifies LTR-RTs with long LTRs writes an output list.
- `LTR-RT/validation/geneOv.sh` finds different cases of overlaps between LTR-RTs and endogenous genes.
- `LTR-RT/family/extract.clust.vmatch.sh` extracts clusters from vmatch output.
- `LTR-RT/family/muscle.cons.seq.sh` aligns multiple sequences for each (sub)family and outputs their consensus sequences (requires `consensus.seq.R`)
- `LTR-RT/family/consensus.seq.R` creates a consensus sequence from a clustalw2 alignment file.
- `LTR-RT/family/extract.INT.dom.bedto.sh` extracts INT protein domains from ltrdigest output and creates a fasta file containing the protein sequence.
- `LTR-RT/family/extract.RT.dom.bedto.sh` extracts RT protein domains from ltrdigest output and creates a fasta file containing the protein sequence.
- `LTR-RT/family/RT_INT_phylo.sh` builds phylogenetic trees of full-length LTR-RTs using RT and INT protein domain sequences.
- `LTR-RT/insertions/RM.LTR.RT.insertion.finder.sh` finds additional insertions in the assembly and classifies them (requires `Parse_RepeatMaskerOutput.R`).
- `LTR-RT/insertions/Parse_RepeatMaskerOutput.R` classifies hits of RepeatMasker output.
- `LTR-RT/analysis/RT_INT_phylo_RLC.sh` and `../RT_INT_phylo_RLG.sh` build phylogenetic trees for Copia and Gypsy, respectively. 
- `LTR-RT/analysis/treeplotting.R` is used to plot phylogenetic trees.
- `LTR-RT/analysis/coverage.R` plots distribution of LTR-RT families across scaffolds and linkage groups.
- `LTR-RT/analysis/fullfamily.R` contains the code used for all the remaining plots and for statistical tests.
-  Required input files for R scripts in `LTR-RT/analysis` can be downloaded here: https://bit.ly/2O2Xnjf
 
