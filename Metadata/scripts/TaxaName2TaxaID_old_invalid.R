# 16th July 2021 Odin Moron-Garcia Aberysytwyth
# Script to generate an intermediate file that links Bacterial names present in the phenotype databases with the assembly_ pressent in the genetic features built by Seb
# It accepts a single folder path and two file names
# [-f]:[not implemented] folder path to search for the files 
#  -g :[Not implemented] genome file - assumed to be a tab separated value with no header plus the columns being an ordered set such that {assembly_accession,species_taxid,ftp_path} See README_assembly_summary_refseq.txt
#  -n :[not implemented] names file  - assumed to be a tab separated value with no header plus the columns being an ordered set such that{Taxa_Name, species_taxid, ...} with ellipsis being a colection of strain names.
#  -v : verbose output

args <- commandArgs(trailingOnly = TRUE)
#parsing options
output <- stdout()
if(any(grepl(x =args,pattern = "-v"))){
	verbose <-  1
}else{
	verbose <-  0
}

if(any(grepl(x =args,pattern = "-o:(.*)"))){
	#output <- grep(x =args, pattern = "-o:(.*)")
	regexpr(text = args ,pattern = "(-o):(?<filename>.*\\b)",perl = T)->re
	 which(re==1)->pos
	substr(x =args,start=attr(re,"capture.start")[pos,][[2]],stop =attr(re,"capture.start")[pos,][[2]]+attr(re,"capture.length")[pos,][[2]])[pos]->output
}

folder <- "/home/ubuntu/Genomes/odin_data"

genome_file <- "genomes_metadata.tsv"
#genome_file <- "assembly_headerless.tsv"
#names_file <- "name_to_taxid.headless.tsv"
names_file <- "name_to_taxid.tsv"

genome_file <- file.path(folder,genome_file)
names_file <- file.path(folder,names_file)

read.table(genome_file,header = F,sep="\t",stringsAsFactors=F,comment.char = "") -> genome
names(genome) <- c("assembly_accession","species_taxid","ftp_path")
#names(genome) <- c("assembly_accession","taxid","species_taxid","TaxaName","Strain","ftp_path")

## if(verbose == 1){write("genome file read",file = stdout())}

read.table(names_file,header = F,sep="\t",stringsAsFactors=F) -> names
names(names)[1:2] <- c("TaxaName","species_taxid")

## if(verbose == 1){write("name file read",file = stdout())}

write.table(merge(genome[1:3],names[1:2],all = F),file = output,row.names =F,quote = F,sep = "\t")

if(verbose == 1){
	write("output: merged file",file = stdout())
	write("The next taxids were not found",file = stdout())
	names[names$species_taxid %in% setdiff(names$species_taxid,genome$species_taxid),]

}

