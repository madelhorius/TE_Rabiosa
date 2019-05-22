# TE_Rabiosa

This repository includes all important scripts that were written for the Master's thesis "dentification of full length transposable elements in Lolium multiflorum."

- 'LTR-RT/validation/selfblast.sh' checks for repeats within LTR-RTs and identifies TEs with long LTRs.
- 'LTR-RT/validation/domains.R' extracts protein domain regions from the ltrdigest gff-output.
- 'LTR-RT/validation/nestedElements.R' finds nested TEs in a gff file, simple cases are exported, whereby the host is split.
- 'LTR-RT/validation/longLTRcheck.R' identifies LTR-RTs with long LTRs writes an output list. 
