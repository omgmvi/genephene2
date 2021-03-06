
# Odin M. Moron-Garcia 
# Date of Creation 21st September, 2021

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


logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO)


## Handling file and folder names - rustic version

Absolute_Folder = "/Users/moronga/Documents/Earlham/trainw2v/W2V_GenePhene2/"
Ortho = "COGS"


if Ortho == "KEGG":
    File_Name = "GenePhene2_Doc2Vec_KEGG.txt"
    Model_Filename = "GenePhene2_Doc2Vec_KEGG.dat"
    Export_Filename = "GenePhene2_Doc2Vec_KEGG.csv"
elif Ortho == "COGS":
    File_Name = "GenePhene2_Doc2Vec_COGS.txt"
    Model_Filename = "GenePhene2_Doc2Vec_COGS.dat"
    Export_Filename = "GenePhene2_Doc2Vec_COGS.csv"

elif Ortho == "pFam":
    File_Name = "GenePhene2_Doc2Vec_pFam.txt"
    Model_Filename = "GenePhene2_Doc2Vec_pFam.dat"
    Export_Filename = "GenePhene2_Doc2Vec_pFam.csv"
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


import smart_open

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

model = gensim.models.doc2vec.Doc2Vec(vector_size=50, min_count=2, epochs=40)

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
        fd.write(model.docvecs.index2entity[idx]+ ";" +  ';'.join([ str(x) for x in model.docvecs[idx]]) + '\n')
    fd.close()
print("The model has been saved at:",Export_File)

print("Script Finished")

