#usage : Add_D2V_header.sh "wildcard"


temp_file="temp"$(date +%F)".tmp"
temp_file2="temp"$(date +%F)"b.tmp"


for file in $1
do
    echo $file

    printf 'GenomeID' > $temp_file
    for i in $(seq 1 1 $(awk 'NR == 1 {print NF-1}' $file));
    do
        printf '\tD2V'$i >> $temp_file
    done
    printf '\n' >> $temp_file
    
    cat $temp_file $file > $temp_file2
    mv $temp_file2 $file
done

rm $temp_file
