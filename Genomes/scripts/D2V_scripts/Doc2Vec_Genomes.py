
# Odin M. Moron-Garcia 
# Date of Creation 21st September, 2021
# Date of major modification 25th January, 2022

# For our purposes we need Doc2Vec since it is suitable for downstream machine learning applications
# And this script will refurbish the original Doc2Vec tutorial so the models are calculated for KEGG, COGs and pFAm
# datasets for the GenePhene2 bacterias data set - and soon for FAPROTAX, MDB and Hiding

#### Input files

# The genomes are generated by Sebastien Radigueau using some custom scripts that result in 
# BED like files that I have collected as "position" "Ortholog gene" files and now put together using

# for file in *KEGG*; do echo -n $file $'\t' >>../GenePhene2_Doc2Vec_KEGG.txt ;cat ${file}|cut -f 2|tail -n +2| tr '\n' $'\t' >> ../GenePhene2_Doc2Vec_KEGG.txt; echo '' >> ../GenePhene2_Doc2Vec_KEGG.txt ;done
#  for file in *COG*; do echo -n $file $'\t' >>../GenePhene2_Doc2Vec_COGS.txt ;cat ${file}|cut -f 2|tail -n +2| tr '\n' $'\t' >> ../GenePhene2_Doc2Vec_COGS.txt; echo '' >> ../GenePhene2_Doc2Vec_COGs.txt ;done
# for file in *pFam*; do echo -n $file $'\t' >>../GenePhene2_Doc2Vec_pFam.txt ;cat ${file}|cut -f 2|tail -n +2| tr '\n' $'\t' >> ../GenePhene2_Doc2Vec_pFam.txt; echo '' >> ../GenePhene2_Doc2Vec_pFam.txt ;done

# That contain a first column with the genome file and the rest of columns

# NOTE: Those scripts put made files like 
# GCF_000021425.1_COG.tsv   "COG0593"  "COG0592"  "COG1195"  "COG0187"  "COG0188" ...
# ...
# With a genome per row and regardless the contigs (I assumed they were long contigs and the effect of concatenating would be minimum)
# But then we notice that not all the contigs are big but there are many small - so this code modify previous one to filter by contig of certain size and treat every contig separately.
# The format is similar but now the same genome can be found in several rows in the file.

# To generate the input files from the BED-Like ones I have now some script named (find the name) that read the BED-like file and "jumps" to a new column every time a new contig appear.

#### Output files

	# Models as binary files like : GenePhene2_Doc2Vec_KEGG.dat
	# Exported also as text files for using in my previous R pipeline - colon separated valued 

#### SCRIPT ####


## Preamble : 


## Script copy-pasted from : https://radimrehurek.com/gensim/auto_examples/tutorials/run_doc2vec_lee.html
# With some modifications for tokenizing and tagging the genome

from os.path import join
import logging
import gensim
import time
import smart_open
import sys

logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO)

## Handling file and folder names - rustic version
Ndim =1000
Ortho = "KEGG"
Geno = "FAPROTAX"

if Geno == "GenePhene":
    Absolute_Folder = "/home/ubuntu/GenePhene2/Genomes/GenePhene2/Doc2Vec_Genome"
    if Ortho == "KEGG":
        File_Name = "GENEPHENE2_D2V_KEGG.tsv"
        Model_Filename = "GenePhene2_Doc2Vec_"+ str(Ndim) +"_KEGG.dat"
        Export_Filename = "GenePhene2_Doc2Vec_"+ str(Ndim) +"_KEGG.tsv"
    elif Ortho == "COGS":
        File_Name = "GENEPHENE2_D2V_COG.tsv"
        Model_Filename = "GenePhene2_Doc2Vec_"+ str(Ndim) +"_COGS.dat"
        Export_Filename = "GenePhene2_Doc2Vec_"+ str(Ndim) +"_COGS.tsv"

    elif Ortho == "pFam":
        File_Name = "GENEPHENE2_D2V_pFam.tsv"
        Model_Filename = "GenePhene2_Doc2Vec_"+ str(Ndim) +"_pFam.dat"
        Export_Filename = "GenePhene2_Doc2Vec_"+ str(Ndim) +"_pFam.tsv"
    else:
        print("No good option")
        exit()

if Geno == "FAPROTAX":  
    Absolute_Folder = "/home/ubuntu/GenePhene2/Genomes/FAPROTAX/Doc2Vec_Genome"
    if Ortho == "KEGG":
        File_Name = "FAPROTAX_SINGLE_D2V_KEGG.tsv"
        Model_Filename = "Faprotax_Doc2Vec_"+ str(Ndim) +"_KEGG.dat"
        Export_Filename = "Faprotax_Doc2Vec_"+ str(Ndim) +"_KEGG.tsv"
    elif Ortho == "COGS":
        File_Name = "FAPROTAX_SINGLE_D2V_COG.tsv"
        Model_Filename = "Faprotax_Doc2Vec_"+ str(Ndim) +"_COGS.dat"
        Export_Filename = "Faprotax_Doc2Vec_"+ str(Ndim) +"_COGS.tsv"

    elif Ortho == "pFam":
        File_Name = "FAPROTAX_SINGLE_D2V_pFam.tsv"
        Model_Filename = "Faprotax_Doc2Vec_"+ str(Ndim) +"_pFam.dat"
        Export_Filename = "Faprotax_Doc2Vec_"+ str(Ndim) +"_pFam.tsv"
    else:
        print("No good option")
        exit()
if Geno == "Hiding":  
    Absolute_Folder = "/home/ubuntu/GenePhene2/Genomes/Hiding/Doc2Vec_Genome"
    if Ortho == "KEGG":
        File_Name = "Hiding_D2V_KEGG.tsv"
        Model_Filename = "Hiding_Doc2Vec_"+ str(Ndim) +"_KEGG.dat"
        Export_Filename = "Hiding_Doc2Vec_"+ str(Ndim) +"_KEGG.tsv"
    elif Ortho == "COGS":
        File_Name = "Hiding_D2V_COG.tsv"
        Model_Filename = "Hiding_Doc2Vec_"+ str(Ndim) +"_COGS.dat"
        Export_Filename = "Hiding_Doc2Vec_"+ str(Ndim) +"_COGS.tsv"

    elif Ortho == "pFam":
        File_Name = "Hiding_D2V_pFam.tsv"
        Model_Filename = "Hiding_Doc2Vec_"+ str(Ndim) +"_pFam.dat"
        Export_Filename = "Hiding_Doc2Vec_"+ str(Ndim) +"_pFam.tsv"
    else:
        print("No good option")
        exit()


Genome_File = join(Absolute_Folder,File_Name)
Model_File = join(Absolute_Folder,Model_Filename)
Export_File = join(Absolute_Folder,Export_Filename)

# Check
print(Genome_File)
print(Model_File)
print(Export_File)


def read_corpus(fname, tokens_only=False):
    with smart_open.open(fname, encoding="iso-8859-1") as f:
        for i, line in enumerate(f):
            tokens = tokenize_genome(line)
            if tokens_only:
                yield tokens
            else:
                # For training data, add tags
                #yield gensim.models.doc2vec.TaggedDocument(tokens, [i])
                #modified to get the tag as the genome ID
                yield gensim.models.doc2vec.TaggedDocument(tokens[1:len(tokens)],[tokens[0]])
## Helper functions

def tokenize_genome(line):
		ll = line.split("\t") #Notice the change that the input files are tab files rather than colon separated
		return ll

##### SCRIPT ####


# Print the time to check how long does it takes

t = time.localtime()
current_time = time.strftime("%H:%M:%S", t)
print(current_time)


## Start processing the genomes

train_corpus = list(read_corpus(Genome_File))

model = gensim.models.doc2vec.Doc2Vec(vector_size=Ndim, min_count=2, epochs=40)

model.build_vocab(train_corpus)

t = time.localtime()
current_time = time.strftime("%H:%M:%S", t)
print(current_time)

# Save the model for later (time for Jupyter?)
model.save(Model_File)
print("The model has been saved at:",Model_File)

# ready for the export
with open(Export_File,'wt') as fd:
    for idx in range(0,len(model.docvecs)):
        fd.write(model.docvecs.index_to_key[idx]+ "\t" +  '\t'.join([ str(x) for x in model.docvecs[idx]]) + '\n')
    fd.close()
print("The model has been saved at:",Export_File)

print("Script Finished")
