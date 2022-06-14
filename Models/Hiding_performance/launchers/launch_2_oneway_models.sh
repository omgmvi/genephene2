Rscript ./code/oneway_models.R DBs/Hiding_COG.reduced.tsv "catalase.positive" "COG\\d+"  "lm_2_COG_KatE.dat" &> lm_2_COG_KatE.log
Rscript ./code/oneway_models.R DBs/Hiding_KEGG.reduced.tsv "catalase.positive" "K\\d+"   "lm_2_KEGG_KatE.dat"&> lm_2_KEGG_KatE.log
Rscript ./code/oneway_models.R DBs/Hiding_pFam.reduced.tsv "catalase.positive" "PF.\\d+" "lm_2_pFam_KatE.dat"&> lm_2_pFam_KatE.log
Rscript ./code/oneway_models.R DBs/Hiding_COG.reduced.tsv "oxidase.positive" "COG\\d+"   "lm_2_COG_Ox.dat"   &> lm_2_COG_Ox.log
Rscript ./code/oneway_models.R DBs/Hiding_KEGG.reduced.tsv "oxidase.positive" "K\\d+"    "lm_2_KEGG_Ox.dat"  &> lm_2_KEGG_Ox.log
Rscript ./code/oneway_models.R DBs/Hiding_pFam.reduced.tsv "oxidase.positive" "PF\\d+"   "lm_2_pFam_Ox.dat"  &> lm_2_pFam_Ox.log
Rscript ./code/oneway_models.R DBs/Hiding_COG.reduced.tsv "acid.phosphatase" "COG\\d+"   "lm_2_COG_phos.dat" &> lm_2_COG_phos.log
Rscript ./code/oneway_models.R DBs/Hiding_KEGG.reduced.tsv "acid.phosphatase" "K\\d+"    "lm_2_KEGG_phos.dat"&> lm_2_KEGG_phos.log
Rscript ./code/oneway_models.R DBs/Hiding_pFam.reduced.tsv "acid.phosphatase" "PF\\d+"   "lm_2_pFam_phos.dat"&> lm_2_pFam_phos.log

