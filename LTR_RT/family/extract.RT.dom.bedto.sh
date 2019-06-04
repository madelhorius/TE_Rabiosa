# Extract RT domains from ltrdigest output and create a fasta file
# Only autonomous TEs are will be used
awk '$3 == "GAG|AP|RT|RH|INT" || $3 == "GAG|AP|INT|RT|RH" {print $1}' LTR_Protein_domains.txt | grep -w -f - ltrdigest.out.gff | grep "ID\|RT" > RT.FL.AUT.gff
# Count number of elements
n="$(grep "ID=" RT.FL.AUT.gff | wc -l)"
n="$(expr $n - 2)"
echo $n
cp RT.FL.AUT.gff tmp
for i in $(eval echo {0..$n} )
do
        ID1="$(grep "ID" tmp | head -1 | awk '{print $9}')"
        ID2="$(grep "ID" tmp | head -2 | tail -1 | awk '{print $9}')"
        sed -n "/$ID1/,/$ID2/p" tmp > tmpLTR
        grep "Parent" tmpLTR | awk '{print $0 "\t" ($5 - $4)}' | sort -nk10,10 | tail -1 > tmp.gff
        faheader="$(head -1 tmpLTR | awk '{print $9}' | awk -F ";" '{print $1}' | sed 's/ID=//g' | sed 's/^/>/g')"
        echo $faheader >> RT.domain.bedto.fa
        bedtools getfasta -fi /data/mahodel/data/Rabiosa_genome_1.0.fa -bed tmp.gff -s | tail -1 >> RT.domain.bedto.fa
        # Remove the first LTR of the gff file
        tail -n+2 tmp > tmp2
        sed -n "/$ID2/,\$p" tmp2 > tmp
        if (($i % 100 == 0))
        then
                echo $i "Elements"
        fi
done
# Last element
grep "Parent" tmp | awk '{print $0 "\t" ($5 - $4)}' | sort -nk10,10 | tail -1 > tmp.gff
faheader="$(head -1 tmp | awk '{print $9}' | awk -F ";" '{print $1}' | sed 's/ID=//g' | sed 's/^/>/g')"
echo $faheader >> RT.domain.bedto.fa
bedtools getfasta -fi /data/mahodel/data/Rabiosa_genome_1.0.fa -bed tmp.gff -s | tail -1 >> RT.domain.bedto.fa
rm tmp tmpLTR tmp2 tmp.gff
echo "DONE"
