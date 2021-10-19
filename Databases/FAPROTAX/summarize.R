#Code similar to Hiding/summarize.R (Time to refactor? Yes)

read.table("FAPROTAX.tsv",header = T,sep = "\t") ->df

summary_df <- apply(df[,!names(df) %in% c("Genus","Species","Name","Strain")],2,
                    function(x){paste("Coverage: ",round(sum(!is.na(x),na.rm=T)/NROW(x)*100,digits = 2))})


write.table(summary_df,"summary_FAPROTAX.txt",sep = "\t",row.names = T)

balance_df<-apply(df[,!names(df) %in% c("Genus","Species","Name","Strain")],2,function(x){paste("Balance :",round(sum(x=="yes",na.rm=T)/sum(!is.na(x),na.rm=T)*100,digits =2),"%")})

write.table(balance_df,"balance_FAPROTAX.txt",sep = "\t",row.names = T)
