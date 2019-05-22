while read p
do
	grep -w $p all.CLSTR | cut -f 1 | grep -f - /data/mahodel/work/ltr/ltrretriever/PARSED_FL_LTR_RT.pass.code.v3.fa -A1 | sed '/--/d' > $p.fa
	grep -w $p all.CLSTR | cut -f 1 | sed 's/_$/;/g' | grep -f - /data/mahodel/work/ltr/ltrretriever/PARSED_FL_LTR_RT.pass.code.v2.gff > $p.gff
	grep "ID=" $p.gff | awk '{print $0 "\t" ($5 - $4)}' > tmp ; mv tmp $p.gff
	mean="$(grep -v N $p.gff | awk '{ total += $10 } END { print total/NR }')"
	awk -v mean="$mean" '$10 < (mean + mean/10) && $10 > (mean - mean/10)' $p.gff | cut -f 9 | awk -F ";" '{print $1}' | sed 's/ID=//g' | sed 's/$/_/g' | grep -f - $p.fa -A1 | sed '/--/d' > $p.len.fa
	/home/mahodel/miniconda2/pkgs/bbmap-38.22-h14c3975_1/bin/reformat.sh in=$p.len.fa out=tmp.fa samplereadstarget=100 ; mv tmp.fa $p.len.fa
	awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < $p.len.fa | tail -n+2 > tmp ; mv tmp $p.len.fa
	# Find sequences with a different orientation than the first sequence in the fasta file:
	head -2 $p.len.fa > head.fa
	blastn -query head.fa -subject $p.len.fa -outfmt 6 -max_hsps 1 | awk '$10 < $9' | cut -f 2 | grep -f - -A1 $p.len.fa | sed '/--/d' > to.rev.$p.len.fa
	blastn -query head.fa -subject $p.len.fa -outfmt 6 -max_hsps 1 | awk '$9 < $10' | cut -f 2 | grep -f - -A1 $p.len.fa | sed '/--/d' > norm.$p.len.fa
	# Get reverse complement
	revseq to.rev.$p.len.fa -reverse -complement -outseq rev.$p.len.fa
	sed -i 's/ Reversed://g' rev.$p.len.fa
	awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < rev.$p.len.fa | tail -n+2 > tmp ; mv tmp rev.$p.len.fa
	# Put reversed and normal sequences together
	cat norm.$p.len.fa > $p.len.fa
	cat rev.$p.len.fa >> $p.len.fa	
	# Align the sequences
	(muscle -in $p.len.fa -out $p.len.aln -clw &&
	# Consensus sequence of the alignment
	sed -i '/\*/d' $p.len.aln &&
	Rscript /home/mahodel/bin/R/consensus.seq.R $p.len.aln $p.cons &&
	sed -i 's/-//g' $p.cons &&
	mv $p.cons muscle/consensus.sequences &&
	mv $p.len.aln muscle/alignments) &
	# Restrict number of parallel processes
	background=( $(jobs -p) )
	if (( ${#background[@]} > 10 )); then
		wait -n
	fi
done < fam.lst
