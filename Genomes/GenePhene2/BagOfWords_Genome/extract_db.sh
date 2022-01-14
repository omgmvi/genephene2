tar -zxvf GenePhene2_BoW.tar.gz

awk 'BEGIN {OFS = "\t"} NR == 1 {if ($1=="Genome_ID"){$1 = "GenomeID"}} 1' GenePhene2_BoW_COG.tsv > tmp
mv tmp GenePhene2_BoW_COG.tsv

awk 'BEGIN {OFS = "\t"} NR == 1 {if ($1=="Genome_ID"){$1 = "GenomeID"}} 1' GenePhene2_BoW_KEGG.tsv > tmp
mv tmp GenePhene2_BoW_KEGG.tsv

awk 'BEGIN {OFS = "\t"} NR == 1 {if ($1=="Genome_ID"){$1 = "GenomeID"}} 1' GenePhene2_BoW_pFam.tsv > tmp
mv tmp GenePhene2_BoW_pFam.tsv
