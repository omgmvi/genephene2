# genephene2
TO DO: WRITE THE README

Introduction : GenePhene2 is meant to be the continuation of the project GenePhene from Prof. Christopher Quince and its collaborators.

The objective is to elaborate a machine learning model able to predict bacteria and archaea metabolic phenotypes from their genomic functional annotation.
More specifically, many bacterias and archaeas cannot be cultivate using standard laboratory technique and therefore knowledge of the phenotypic characteristics are mostly unknown.
On the other hand modern sequencing technologies and the bioinformatic tools associated has facilitated the recovery and reconstruction of 'metagenomic assembled genomes". 
That would mean that we are able to observe the genomes of unknown microorganims but our knowledge of their activities is only provided by comparison with other siblings taxa through RNA16S or similar barcoding.

The approach to be taken here is to link the functionally annotated genomes, through ORF finding and Ortholog detection algorithms, to more or less specific phenotypes such like 'aerobe' or 'catalase activity'. This linkage is being built on machine learning algorithms for supervised classification where the phenotypes play the role of predicted variables and the Ortholog gene name the role of predictors. To advance the current tools we plant to transform the annotated gene names in a continuous numberical variable using NeuroLinguistic programming tools like Word2Vec and Doc2Vec instead of the so called bag of words approximation.
