# Odin M. Moron-Garcia 
# Date of Creation 21st September, 2021


# Goal:
### Script and functions to pass BED-like files with Contig-position Ortholog columns onto a list of orthologs
# valid for Word2Vec

# Overview -  Word2Vec fitting algorithm expect a text written as sentences where very sentence is a list of strings - tokenized and the texts is then a list of list.
# Yet, according to documentation, it is recommended to use python generators to run over the sentences.
# We need to remaind the problem in our hands are around 3000 microorganims and 3000 genes and so far
# we consider a genome like a sentence (I am hoping this work with word2vec methodology) so 3000 thousand sentences 3000-5000 words long

# (Although, our first examples GenePhene2 only contain 45 bacterias and the same 3000 - 5000 genes)

# To convert all BED-like file to a single file with a genome per line and genes ordered by genome position
# is easier to do using bash

# When file are separated in folders acording to the ortholog type
# for file in *; do cat ${file}|cut -f 2|tail -n +2| tr '\n' ';' >> ../COGS.txt; echo $'' >> ../COGS.txt; done

# if all ortholog types are in the same folder

# for file in *KEGG*; do cat ${file}|cut -f 2|tail -n +2| tr '\n' ';' >> ../KEGG.txt; echo $'' >> ../KEGG.txt; done
# for file in *pFam*; do cat ${file}|cut -f 2|tail -n +2| tr '\n' ';' >> ../pFam.txt; echo $'' >> ../pFam.txt; done

# Following the example tutorial at Radim Rehurek's web for Gensim - making a Generator class is among the best option to read
# and make the coding easier.


######## Preamble #######
from gensim import utils
from gensim.models import Word2Vec
from os.path import join
import tempfile
import time


##### Input variables 


Absolute_Folder = "/Users/moronga/Documents/Earlham/trainw2v/W2V_GenePhene2/"

#File_Name = "GenePhene2_KEGG.txt"
#Model_Filename = "GenePhene2_Word2Vec_KEGG.dat"

#File_Name = "GenePhene2_COGS.txt"
#Model_Filename = "GenePhene2_Word2Vec_COGS.dat"

File_Name = "GenePhene2_pFam.txt"
Model_Filename = "GenePhene2_Word2Vec_pFam.dat"

Genome_File = join(Absolute_Folder,File_Name)
Model_File = join(Absolute_Folder,Model_Filename)

# Check
print(Genome_File)
print(Model_File)


###### Classes

# Copied from : https://radimrehurek.com/gensim/auto_examples/tutorials/run_word2vec.html#training-your-own-model

class Genetic_Corpus:
    """An iterator that yields sentences (lists of str)."""

    def __iter__(self):
        for line in open(Genome_File):
            # assume there's one document per line, tokens separated by whitespace
            yield tokenize_genome(line)

## Helper functions

def tokenize_genome(line):
		ll = line.split(";")
		return ll


##### SCRIPT ####


# Print the time to check how long does it takes

t = time.localtime()
current_time = time.strftime("%H:%M:%S", t)
print(current_time)


## Start processing the genes


sentences = Genetic_Corpus()
model = Word2Vec(sentences=sentences)


# Print again the time to checko how long did it take


t = time.localtime()
current_time = time.strftime("%H:%M:%S", t)
print(current_time)


# Save the model for later (time for Jupyter?)
model.save(Model_File)
print("The model has been saved at:",Model_File)
