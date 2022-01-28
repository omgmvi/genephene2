tar -zxvf FAPROTAX_BoW.tar.gz

awk 'BEGIN {OFS = "\t"} NR == 1 {if ($1=="Genome_ID"){$1 = "GenomeID"}} 1' FAPROTAX_BoW_COG.tsv > tmp
mv tmp FAPROTAX_BoW_COG.tsv

awk 'BEGIN {OFS = "\t"} NR == 1 {if ($1=="Genome_ID"){$1 = "GenomeID"}} 1' FAPROTAX_BoW_KEGG.tsv > tmp
mv tmp FAPROTAX_BoW_KEGG.tsv

awk 'BEGIN {OFS = "\t"} NR == 1 {if ($1=="Genome_ID"){$1 = "GenomeID"}} 1' FAPROTAX_BoW_pFam.tsv > tmp
mv tmp FAPROTAX_BoW_pFam.tsv
