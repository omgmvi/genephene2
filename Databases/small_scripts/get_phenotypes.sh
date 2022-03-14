# If you have a line of metabolic traits (columns of a table) and you want to get them as an array
cat Metabolic.tsv |head -n 1|cut -f 1-4 --complement|sed 's/[\t]/\","/g'|sed 's/^./"&/g'|sed 's/.$/&"/g'
