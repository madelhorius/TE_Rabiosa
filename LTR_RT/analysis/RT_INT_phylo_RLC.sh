# Create a phylogenetic tree based on INT and RT protein sequence
# Loop: 2 domains seperately
for d in INT RT
do
	sed 's/>//g' /data/mahodel/work/ltr/ltrretriever/ltrdigest/$d.domain.bedto.fa > $d.domain.bigfam.fa
	awk 'NR==FNR {code[$1]=$2; next} {for (i=1; i<=NF; i++) if ($i in code) $i=code[$i]; print}' ltr.index.fam2.info  $d.domain.bigfam.fa > tmp ; mv tmp $d.domain.bigfam.fa
	sed -i '1~2s/^/>/' $d.domain.bigfam.fa
	# Only take the 30 biggest families
	grep ">" $d.domain.bigfam.fa | grep RLC | awk -F ";" '{print $7}' | sort | sed '/^$/d' | uniq -c | sort -nk1,1 | tail -30 | awk '{print $2}' > big.fam.FL.AUT.lst
	grep -w -Ff big.fam.FL.AUT.lst $d.domain.bigfam.fa -A1 | sed '/--/d' > tmp; mv tmp $d.domain.bigfam.fa
	# Remove nested elements
	grep ">" $d.domain.bigfam.fa | grep -v "TN\|tN\|Tn\|tn" > not.nested
	grep -Ff not.nested $d.domain.bigfam.fa -A1 | sed '/--/d' > tmp.fa ; mv tmp.fa $d.domain.bigfam.fa
	# Translate: DNA -> protein
	transeq $d.domain.bigfam.fa $d.domain.bigfam.pep
	sed -i 's/_1$//g' $d.domain.bigfam.pep
done
sed -i 's/X//g' RT.domain.bigfam.pep INT.domain.bigfam.pep
awk -F ";" '{print $1 ";" $7 ";" $6}' INT.domain.bigfam.pep | sed 's/;;$//g' > tmp ; mv tmp INT.domain.bigfam.pep
awk -F ";" '{print $1 ";" $7 ";" $6}' RT.domain.bigfam.pep | sed 's/;;$//g' > tmp ; mv tmp RT.domain.bigfam.pep
# Merge files
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < INT.domain.bigfam.pep | tail -n+2 > tmp ; mv tmp INT.domain.bigfam.pep
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < RT.domain.bigfam.pep | tail -n+2 > tmp ; mv tmp RT.domain.bigfam.pep
paste -d "_" RT.domain.bigfam.pep INT.domain.bigfam.pep | awk -F "_>" '{print $1}' > RT.INT.domain.bigfam.pep
lines="$(grep ">" RT.INT.domain.bigfam.pep | wc -l)"
prop="$(echo "800 / $lines" | bc -l)" 
# Take a subsample of the sequences to build the tree script from QIIME package
subsample_fasta.py -i RT.INT.domain.bigfam.pep -o RLC.subsample.RT.INT.domain.bigfam.pep -p $prop
cat RLG.compar.pep >> RLC.subsample.RT.INT.domain.bigfam.pep
# Tree building
muscle -in RLC.subsample.RT.INT.domain.bigfam.pep -out RLC.subsample.RT.INT.domain.bigfam.aln -clwstrict
clustalw2 -INFILE=RLC.subsample.RT.INT.domain.bigfam.aln -TREE
