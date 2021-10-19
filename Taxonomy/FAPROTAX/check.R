

read.table("SpeciesTaxID.tsv",header = T,sep = "\t")->SpTxDB

read.table("../names.dmp",header = F,sep = "|",strip.white = T, comment.char = "",fill=T,quote = "")->TaxDB
TaxDB <- TaxDB[,1:2];
names(TaxDB) <- c("TaxID","Name")

subset(SpTxDB,Name !="Bacteria")->SpTxDB
merge(SpTxDB,TaxDB,by = "TaxID",all = F)->tocheck
