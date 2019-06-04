# Create trees for for two different clustering results based on INT and RT protein sequence (only FL-TEs)
for f in 75.vmatchClusters.CLSTR  80.vmatchClusters.CLSTR
do
	# Connect ltrdigest output and families from clustering
	sed 's/_//g' $f | sed 's/LTRRT/LTR_RT/g' > tmp
	awk 'NR==FNR {code[$1]=$2; next} {for (i=1; i<=NF; i++) if ($i in code) $i=code[$i]; print}' tmp /data/mahodel/work/ltr/ltrretriever/ltrdigest/ltr.index.info > ltr.index.fam2.info
	awk '{print $1 "\t" $3 ";" $2}' ltr.index.fam2.info > tmp ; mv tmp ltr.index.fam2.info
	# Loop: 2 domains seperately
	for d in INT RT
	do
		sed 's/>//g' /data/mahodel/work/ltr/ltrretriever/ltrdigest/$d.domain.bedto.fa > $d.domain.bigfam.fa
		awk 'NR==FNR {code[$1]=$2; next} {for (i=1; i<=NF; i++) if ($i in code) $i=code[$i]; print}' ltr.index.fam2.info  $d.domain.bigfam.fa > tmp ; mv tmp $d.domain.bigfam.fa
		sed -i 's/LTR/>LTR/g' $d.domain.bigfam.fa
		# Only take the 30 biggest families
		grep ">" $d.domain.bigfam.fa | awk -F ";" '{print $7}' | sort | uniq -c | sort -nk1,1 | tail -30 | awk '{print $2}' > big.fam.FL.AUT.lst
		grep -w -f big.fam.FL.AUT.lst $d.domain.bigfam.fa -A1 | sed '/--/d' > tmp; mv tmp $d.domain.bigfam.fa
		# Remove nested elements
		grep ">" $d.domain.bigfam.fa | grep -v "TN\|tN\|Tn\|tn" > not.nested
		grep -f not.nested $d.domain.bigfam.fa -A1 | sed '/--/d' > tmp.fa ; mv tmp.fa $d.domain.bigfam.fa
		# Translate: DNA -> protein
		transeq $d.domain.bigfam.fa $d.domain.bigfam.pep
		sed -i 's/_1//g' $d.domain.bigfam.pep
	done
	sed -i 's/X//g' RT.domain.bigfam.pep INT.domain.bigfam.pep
	awk -F ";" '{print $1 ";" $7}' INT.domain.bigfam.pep | sed 's/;$//g' | sed 's/_1//g' > tmp ; mv tmp INT.domain.bigfam.pep
	awk -F ";" '{print $1 ";" $7}' RT.domain.bigfam.pep | sed 's/;$//g' | sed 's/_1//g' > tmp ; mv tmp RT.domain.bigfam.pep
	# Merge files
	awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < INT.domain.bigfam.pep | tail -n+2 > tmp ; mv tmp INT.domain.bigfam.pep
	awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < RT.domain.bigfam.pep | tail -n+2 > tmp ; mv tmp RT.domain.bigfam.pep
	f=$(sed 's/..\///g' <<< $f)
	f=$(sed 's/\.CLSTR//g' <<< $f)
	paste -d "_" RT.domain.bigfam.pep INT.domain.bigfam.pep | sed 's/_>LTR*.*$//g' | sed 's/_//g' | sed 's/LTRRT/LTR_RT/g' > RT.INT.domain.bigfam.pep
	# Take a subsample of the sequences to build the tree script from QIIME package
	subsample_fasta.py -i RT.INT.domain.bigfam.pep -o $f.subsample.RT.INT.domain.bigfam.pep -p 0.1
	# Tree building
	(clustalw2 -INFILE=$f.subsample.RT.INT.domain.bigfam.pep -ALIGN -OUTFILE=$f.subsample.RT.INT.domain.bigfam.aln; clustalw2 -INFILE=$f.subsample.RT.INT.domain.bigfam.aln -TREE) &
done
