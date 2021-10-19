folder <- "~/Metadata/refseq/"
file <- "assembly_summary_refseq.txt"


read.table(file.path(folder,file),comment.char = "",header = F,nrow = 1,skip = 1,stringsAsFactor = F,sep = "\t")->header
header[1] <- "assembly_accession"
read.table(file.path(folder,file) ,skip = 2,sep = "\t",header = F,comment.char="",quote = "") ->Refseq
names(Refseq) <- header
Refseq[,c("taxid","species_taxid","organism_name","ftp_path","assembly_accession")] -> Refseq

read.table("~/Genomes/FAPROTAX/BagOfWords_Genome/KEGG_genome.tsv",sep = "\t")->Genome
Genome[1] ->Genome
names(Genome) <- "taxid"
Genome$GenomeID <- row.names(Genome)

read.table("~/Metadata/FAPROTAX/TaxID2Genome_metadata.tsv",sep = "\t",header = T)->Metadata


merge(Genome,Refseq,by.x = "GenomeID",by.y = "assembly_accession")->res1
merge(Genome,Metadata,by.x = "GenomeID",by.y = "assembly_accession")->res2

length(unique(Genome$GenomeID))#9.407
length(unique(Metadata$assembly_accession))#128.834

length(unique(res1$GenomeID)) # 8.542 Where the fuck are the rest? If these genomes are not in the Refseq?
length(unique(res2$GenomeID)) ##552


# What is driving me crazy is that it seems that the genome downloaded here from Chris previous version of the project does not looks like coming from the Refseq file
# in fact there are 
sum(Genome$GenomeID %in% Refseq$assembly_accession) # 8.542 genomes in the Refseq and 
sum(!Genome$GenomeID %in% Refseq$assembly_accession)# 865 not in there
# Out of these 865 genomes not found, there are 603 whose  taxid are neither in the Refseq taxid or species_taxid and 262 are there ...
sum(!unique(subset(Genome,!Genome$GenomeID %in% Refseq$assembly_accession)[["taxid"]]) %in% unique(c(Refseq$species_taxid,Refseq$taxid)))
# But from the original 9.407 genomes, 871 has a taxid that is not either in the Refseq species_taxid or taxid columns
#Not surprinsingly, 8536 are there
sum(unique(Genome[["taxid"]]) %in% unique(c(Refseq$species_taxid,Refseq$taxid)))

#This make me think that their original data set was broken???

# if 8536 genomes are in the refseq why my Metadata file does not match with them? - it may be my Metadata?

nrow(subset(res1,subset = taxid.x == taxid.y | taxid.x == species_taxid)) # 8.224 rows has a coincidence between taxid at Genome and the any of the Refseq taxid's - Curiously, 
nrow(subset(res1,subset = taxid.x == taxid.y))# 5.483 coincide in the taxid and 
nrow(subset(res1,subset = taxid.x == species_taxid)) # 8.219 do in the species_taxid (some would do in both
nrow(subset(res1,subset = taxid.x == taxid.y & taxid.x == species_taxid)) # 5478
# So, it happen that only 318 does not coincide in the taxids
# of wich 268 are not in the refseq id's at all
sum(!(subset(res1[1:4],subset = taxid.x != taxid.y & taxid.x != species_taxid))$taxid.x %in% c(Refseq$taxid,Refseq$species_taxid)) #268

# Ok, then what happen with my Metadata?
# Res2 contains those Genomes that are in both, my Metadata and the Genome
# And we can count how many of the taxids coincides
nrow(subset(res2,subset = taxid.x == taxid.y | taxid.x == species_taxid)) # 486 
nrow(subset(res2,subset = taxid.x != taxid.y & taxid.x != species_taxid))# 129
# But the question here is why there's only 615 genomeID coincidences? what happen with the others?

# Let's search those sequences that are not in the Metadata 
setdiff(unique(Genome$GenomeID),unique(Metadata$assembly_accession))->GenNotMet
#length(GenNotMet) 8855
# most of them (7990) has the Genome in the metadata and the taxids coincide
nrow(subset(res1,subset = GenomeID %in% GenNotMet))
# I'm not understanding what is going on, are those taxid the correct ones? are they in the FAPROTAX?
#I'll pick one 8519 GCF_900177705.1 1938898 1938898       1938898
subset(Refseq, taxid == '1938898')
#                 taxid     species_taxid   organism_name           assembly_accession
#         215055 1938898       1938898      Fibrobacter sp. UWB15    GCF_900177705.1
#AHAAAA, Fibrobacter sp. is a genera and I do not include them!

# So 1 reason to not appear, being a genus rather than a species.

# Let's try another one - 8526 GCF_900177815.1 1852522 1852522       1852522
subset(Refseq, taxid == '1852522')
#                taxid          species_taxid   organism_name               assembly_accession
#         163021 1852522       1852522          Paenibacillus aquistagni    GCF_012926515.1
#         215066 1852522       1852522          Paenibacillus aquistagni    GCF_900177815.1
#         215739 1852522       1852522          Paenibacillus aquistagni    GCF_900188525.1

subset(Metadata, taxid == '1852522')
#[1] assembly_accession taxid              species_taxid      organism_name
#<0 rows> (or 0-length row.names)

#So why this bacteria is not in the Metadata?
levels(Metadata$organism_name)[1250:1500][7:41]
#[1] "Paenibacillus amylolyticus"      "Paenibacillus anaericanus"
#[3] "Paenibacillus antarcticus"       "Paenibacillus apiarius"
#[5] "Paenibacillus barengoltzii"      "Paenibacillus brasilensis"
#[7] "Paenibacillus cellulosilyticus"  "Paenibacillus chibensis"
#[9] "Paenibacillus chitinolyticus"    "Paenibacillus cineris"
#[11] "Paenibacillus cookii"            "Paenibacillus curdlanolyticus"
#[13] "Paenibacillus darwinianus"       "Paenibacillus durus"
#[15] "Paenibacillus elgii"             "Paenibacillus favisporus"
#[17] "Paenibacillus glucanolyticus"    "Paenibacillus graminis"
#[19] "Paenibacillus jamilae"           "Paenibacillus kribbensis"
#[21] "Paenibacillus lautus"            "Paenibacillus macerans"
#[23] "Paenibacillus naphthalenovorans" "Paenibacillus odorifer"
#[25] "Paenibacillus pabuli"            "Paenibacillus peoriae"
#[27] "Paenibacillus phyllosphaerae"    "Paenibacillus polymyxa"
#[29] "Paenibacillus rhizosphaerae"     "Paenibacillus terrae"
#[31] "Paenibacillus thiaminolyticus"   "Paenibacillus turicensis"
#[33] "Paenibacillus validus"           "Paenibacillus wynnii"
#[35] "Paenibacillus xylanilyticus"
#Why it is not in the list of Paenibacillus? - I am checking manually just in case this name is mispelled - No, it is not there - explanation? The database is not for FAPROTAX

# Another example


subset(res1,subset = GenomeID %in% GenNotMet)[9,]
#                    GenomeID taxid.x taxid.y    species_taxid   organism_name
#         9 GCF_000006765.1     287  208964      287             Pseudomonas aeruginosa PAO1

subset(Refseq, taxid == '287')
# Lots of records I replicate the one that is like the res1
subset(Refseq, taxid == '287' & assembly_accession == 'GCF_000006765.1')
#[1] taxid              species_taxid      organism_name      assembly_accession
#<0 rows> (or 0-length row.names)
#It is not there?? 
subset(Refseq, species_taxid == '287' & assembly_accession == 'GCF_000006765.1')
#       taxid       species_taxid       organism_name                   assembly_accession
#    86 208964      287                 Pseudomonas aeruginosa PAO1    GCF_000006765.1

#In metadatadat there are many Pseudomonas aeruginosa but none with the GCF_000006765.1
subset(Metadata, (taxid == '287' | species_taxid == '287') )
#GCF_009830215.1   287           287 Pseudomonas aeruginosa

subset(Metadata, assembly_accession == 'GCF_000006765.1' )
#[1] assembly_accession taxid              species_taxid      organism_name
#<0 rows> (or 0-length row.names)

# So, there's an entry in Metadata for taxid = 287 and an entry in Refseq for 287 but the script did not include that genome?
setdiff(subset(Refseq,subset = taxid == 287 | species_taxid == 287)[4],subset(Metadata,subset = taxid == 287 | species_taxid == 287))
intersect(subset(Refseq,subset = taxid == 287 | species_taxid == 287)[4],subset(Metadata,subset = taxid == 287 | species_taxid == 287))
#data frame with 0 columns and 0 rows
setdiff(subset(Metadata,subset = taxid == 287 | species_taxid == 287),subset(Refseq,subset = taxid == 287 | species_taxid == 287)[4])
#There's something weird with my Metadata code - I am going to have a look
