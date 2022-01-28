# usage : bash ./remove_contig_files.sh "GENEPHENE2_*tsv" 10
# Do not forget to quote around the wildcar expansion

temp_file="temp"$(date +%F)".tmp"

for file in $1
do
    echo $file
    bash remove_contig.sh "$file" "$2" > $temp_file
    cat $temp_file > $file
done
rm $temp_file
