
Nmin="$2"
awk -v Nmin=$Nmin '{if (NF >= Nmin+1) print $_}' "$1"
