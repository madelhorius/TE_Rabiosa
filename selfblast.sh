# Create list files
touch LTR_RT.discarded.lst LTR_RT.LTRlong.lst LTR_RT.pass.lst
# Loop over all cadidates
for i in {1..44292}
do
	name="$(grep ">" PARSED_FL_LTR_RT.fa | head "-$i" | tail -1 | sed 's/>//g')"
	grep -A1 $name PARSED_FL_LTR_RT.fa > ${name}.fa
	blastn -query ${name}.fa -subject ${name}.fa -evalue 1E-200 -outfmt 6 -out blastout.txt
	nrRep=$(wc -l blastout.txt | awk '{print $1}')
	# If more than 3 hits => repeats present.
	if (("$nrRep" >  3 ))
	then
		echo $name >> LTR_RT.discarded.lst 
	else
		if (("$nrRep" >  2 ))
		then
			tot=$(awk '{print $8}' blastout.txt | head -1)
			LTRstart=$(awk '{print $7}' blastout.txt | head -2 | tail -1)
			LTRend=$(awk '{print $8}' blastout.txt | head -2 | tail -1)
			LTRlen=$(($LTRend - $LTRstart))
			ratio=$(bc <<<"scale=2;$LTRlen/$tot*100")	
			if (( $(echo "$ratio > 40" |bc -l) ))
			then
				echo $name >> LTR_RT.LTRlong.lst
			else
			echo $name >> LTR_RT.pass.lst
			fi
		else
			echo $name >> LTR_RT.pass.lst
		fi

	fi
	rm ${name}.fa
	mod=$(($i%50))
	if (("$mod" == 0))
	then
		echo "Element" $i
	fi	
done
rm blastout.txt
