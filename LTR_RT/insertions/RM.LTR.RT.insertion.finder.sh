# Bash script to find LTR-copies in genome
# Using Repeatmasker, bedtools intersect and R

# declare an array variable
declare -a arr=("split_0" "split_1"  "split_2"  "split_3"  "split_4"  "split_5"  "split_6"  "split_7"  "split_8"  "split_9")
for i in "${arr[@]}"
do
	cd $i
	counter=1
	while read p
	do
    		# Get consensus sequences of families
    		grep -A1 $p /data/mahodel/work/ltr/ltrretriever/families/modnames.all.cons.ltr.int.fa > tmp.fa
    		sizeltr="$(grep _LTR tmp.fa -A1 | tail -1 | wc -c)"
    		# Run RepeatMasker with this one family as library
    		RepeatMasker -pa 16 -qq -norna -no_is -cutoff 250 -gccalc -engine ncbi -nolow -lib tmp.fa  $i_\Rabiosa_genome_1.0.fa 
    		# Remove output that is not needed
    		skip="$(wc -l split_9_Rabiosa_genome_1.0.fa.out | awk '{print $1}')"
    		if ((skip > 1))
    		then
    			rm $i\_Rabiosa_genome_1.0.fa.tbl $i\__Rabiosa_genome_1.0.fa.masked $i\_Rabiosa_genome_1.0.fa.cat.gz
    			p=$(sed 's/>//g' <<< $p)
    			mv $i\_Rabiosa_genome_1.0.fa.out $p.RM.out
    			sed -i 's/*//g' $p.RM.out
    			# Run Rscript to parse RM-Output
    			Rscript /home/mahodel/bin/R/Parse_RepeatMaskerOutput.R $p.RM.out $p.RM.parsed.gff $sizeltr

    			# Find overlaps and remove TE with lower score (Starting from round 2)
    			if ((counter > 1))
    			then
        			grep "part\|Solo" $p.RM.parsed.gff > parts.new.gff
        			grep "part\|Solo" LTR.copies.gff > parts.old.gff
        			bedtools intersect -a parts.old.gff -b parts.new.gff -wo -f 0.01 -r > intersect.txt
        			# If FL vs. Solo-LTR / Truncated TE keep FL
        			awk '$9 !~ /Trunc|Solo/  && $18 ~ /Trunc|Solo/' intersect.txt | awk '{print $18}' | awk -F ":" '{print $1}' > torm.new.lst
        			awk '$18 !~ /Trunc|Solo/  && $9 ~ /Trunc|Solo/' intersect.txt | awk '{print $9}' | awk -F ":" '{print $1}' > torm.old.lst
        			grep -Fwv -f torm.old.lst intersect.txt > tmp ; mv tmp intersect.txt
        			grep -Fwv -f torm.new.lst intersect.txt > tmp ; mv tmp intersect.txt
        			# Check score, keep the one with the higher score
        			awk '$6 > $15' intersect.txt | awk '{print $18}' | awk -F ":" '{print $1}' >> torm.new.lst
        			awk '$15 > $6' intersect.txt | awk '{print $9}' | awk -F ":" '{print $1}' >> torm.old.lst
        			grep -Fwv -f torm.old.lst LTR.copies.gff > tmp ; mv tmp LTR.copies.gff
        			grep -Fwv -f torm.new.lst $p.RM.parsed.gff >> LTR.copies.gff
    			else # First round(s): Nothing to be compared
        			cat $p.RM.parsed.gff > LTR.copies.gff
    			fi
    			mv $p.RM.parsed.gff $p.RM.out parsedRMoutput
    		else
    			mv $p.RM.out parsedRMoutput
			echo "No hits found for" $p
    		fi
    		# Counter
    		nrlines="$(wc -l LTR.copies.gff | awk '{print $1}')"
    		if ((nrlines > 0))
    		then
        		counter=$[$counter +1]
    		fi
    		echo $counter
	done < /data/mahodel/work/ltr/ltrretriever/families/modnames.all.cons.names.lst
	rm torm.old.lst torm.new.lst intersect.txt parts.old.gff parts.new.gff
	cd ..
done

