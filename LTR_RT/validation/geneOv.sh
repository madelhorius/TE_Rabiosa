## Commands to find different cases of overlaps

# Case 1 & 2: Gene is completly inside an LTR-RT
bedtools intersect -a PARSED_FL_LTR_RT.pass.FL.gff -b genes.Rabiosa.1.0.gff -F 1.0 -wa -u > FL_LTR_RT.case1.case2.gff
# Case 1: gene is completely inside LTR-RT, gene does not contain any introns
bedtools intersect -a FL_LTR_RT.case1.case2.gff -b cds.Rabiosa.1.0.gff -c | awk 'BEGIN { OFS = "\t" } $10 == "1" {NF--; print}' > FL_LTR_RT.case1.gff
# Case 2: gene is completely inside LTR-RT. gene contains introns
bedtools intersect -a FL_LTR_RT.case1.case2.gff -b cds.Rabiosa.1.0.gff -c | awk 'BEGIN { OFS = "\t" } $10 != "1" {NF--; print}' > FL_LTR_RT.case2.gff
bedtools intersect -a FL_LTR_RT.case1.gff -b cds.Rabiosa.1.0.gff -wo > overlaps.case1.txt
bedtools intersect -a FL_LTR_RT.case2.gff -b cds.Rabiosa.1.0.gff -wo > overlaps.case2.txt

# Case 3,4,5 & 6: Gene is partially overlapping with LTR-RT
bedtools intersect -a PARSED_FL_LTR_RT.pass.FL.gff -b genes.Rabiosa.1.0.gff -wa -u | grep -v -f FL_LTR_RT.case1.case2.gff > FL_LTR_RT.case3.case4.case5.case6.gff
# Case 3: LTR-RT contains full coding sequence(s) of a gene, however it is just partially overlapping with the full gene 
bedtools intersect -a FL_LTR_RT.case3.case4.case5.case6.gff -b cds.Rabiosa.1.0.gff -F 1.0 -wa -u > FL_LTR_RT.case3.gff
# Case 4: LTR-RT is partially overlapping with a coding sequence of a gene
bedtools intersect -a FL_LTR_RT.case3.case4.case5.case6.gff -b cds.Rabiosa.1.0.gff -wa -u | grep -v -f FL_LTR_RT.case3.gff > FL_LTR_RT.case4.gff
grep -v -f FL_LTR_RT.case3.gff -f FL_LTR_RT.case4.gff FL_LTR_RT.case3.case4.case5.case6.gff > FL_LTR_RT.case5.case6.gff
# Case 6: LTR-RT just overlaps with the UTR of a gene
bedtools intersect -a FL_LTR_RT.case5.case6.gff -b utr.Rabiosa.1.0.gff -wa -u > FL_LTR_RT.case6.gff
# Case 5: LTR-RT is completly inside an intron of a gene
grep -v -f FL_LTR_RT.case6.gff FL_LTR_RT.case5.case6.gff > FL_LTR_RT.case5.gff
bedtools intersect -a FL_LTR_RT.case3.gff -b cds.Rabiosa.1.0.gff -wo > overlaps.case3.txt
bedtools intersect -a FL_LTR_RT.case4.gff -b cds.Rabiosa.1.0.gff -wo > overlaps.case4.txt
bedtools intersect -a FL_LTR_RT.case5.gff -b genes.Rabiosa.1.0.gff -wo > overlaps.case5.txt
bedtools intersect -a FL_LTR_RT.case6.gff -b utr.Rabiosa.1.0.gff -wo > overlaps.case6.txt

# Remove cases 2,3,4 because we cannot be sure if the LTR-RT or the gene is true or maybe their both false positive
mkdir gff
for file in gff/*; do awk '{print $9}' $file | awk -F ";" '{print $1}' | sed 's/ID=//g' | sed 's/$/;/' > ${file}.lst ; done
cp PARSED_FL_LTR_RT.pass.gff 2_PARSED_FL_LTR_RT.pass.gff
sed -i 's/$/;/' 2_PARSED_FL_LTR_RT.pass.gff
mv gff/* .
grep -v -f gene.ltr/FL_LTR_RT.case2.gff.lst -f gene.ltr/FL_LTR_RT.case3.gff.lst -f gene.ltr/FL_LTR_RT.case4.gff.lst 2_PARSED_FL_LTR_RT.pass.gff > tmp ; mv tmp 2_PARSED_FL_LTR_RT.pass.gff

