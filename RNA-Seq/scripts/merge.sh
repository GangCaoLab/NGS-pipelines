files=$@

echo -n -e 'gene\t'

ids=$(echo $files | sed 's/.count.txt//g')
echo "$(echo $ids | sed 's/ /\t/g')"

awk '
NF > 1 { a[$1] = a[$1]"\t"$2 }
END {for( i in a ) print i a[i]}' $files