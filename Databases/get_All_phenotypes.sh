#tail -n +1 *phenotypes.txt > All_phenotypes.txt
for i in $(ls ./*/*phenotypes.txt) ;do echo \<===== $i ; cat $i|nl ; done > All_phenotypes.txt
