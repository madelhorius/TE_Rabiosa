# Extract Clusters from vmatch output
n="$(grep "\:" 70.vmatchClusters | wc -l)"
n="$(expr $n - 2)"
echo $n
for i in $(eval echo {0..$n} )
do
	ii="$(expr $i + 1)"
	sed -n "/^$i\:$/,/^$ii\:$/p" 70.vmatchClusters | head -n-1 | awk -F "sc" '{print $1}' |  tail -n+2 | sed 's/  //g' | awk -v i="$i" '{print $1 "\t" "F" i}' >> 70.vmatchClusters.CLSTR
done
# Last cluster
sed -n "/^$ii\:/,\$p" 70.vmatchClusters | tail -n+2 | awk -F "sc" '{print $1}' | sed 's/  //g' |  awk -v i="$ii" '{print $1 "\t" "F" i}' >> 70.vmatchClusters.CLSTR
