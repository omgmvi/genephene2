### Odin Moron-Garcia - 19th May 2022 ### Modification from lm_oneway_model on 23rd May 2022

### Objective #####
# The mini-goal, just perform a logistic regresion univariate on each feature (gene) separately


########################################
# FUNCTIONS
logging <- function(...)cat(sprintf(paste0(...)),sep = '\n',file = stderr())

df_replace_NA <- function(df,column,replace)
{as.character(df[[column]])->tmp; tmp[is.na(tmp)]<- replace;df[[column]] <-tmp;df[[column]]<-as.factor(df[[column]]);df}

mI_categorical <-function(x,y){
    tx <- prop.table(table(x))
    ty <- prop.table(table(y))
    txy <- prop.table(table(x,y))

    X <- as.character(unique(x))
    Y <- as.character(unique(y))

    listXY <- unlist(apply(expand.grid(X,Y),1,list),recursive = F)

    lapply(listXY,FUN = function(p,px,py,pxy){
        if(pxy[p[1],p[2]] !=0){
            pxy[p[1],p[2]]*log(pxy[p[1],p[2]]/(px[ p[1] ]*py[ p[2] ]),base =2)
        }else{0}
    },px = tx,py = ty,pxy = txy)->elem
    sum(unlist(elem))
}
######################################
## Options - To be parsed from CommandArgs()

commandArgs(trailingOnly = T)->opt
if(length(opt) !=4){
    print("usage: Rscript oneway_models.R : Rscript oneway_models.R ./DBs/Hiding_COG_reduce.tsv \'catalase.positive\' \'COG\\d+\' output_COG_catalase.dat")
    stop()
}
DB <- opt[1]
Phenotype<- opt[2]
Genotype_regexp <- opt[3]
output_models <- opt[4]

logging("Database:",DB)
logging("Phenotype:",Phenotype)
logging("Regexp for find gene columns:",Genotype_regexp)
logging("Output file:",output_models)

########################################
## checks ###
stopifnot(file.exists(DB))
###########

##### Reading data

read.table(file = DB,header = T,sep="\t") ->db
# collecting the genome names

grep(pattern = Genotype_regexp,x = names(db),value =T)->Genome

#### checks ####
# db is not empty

# The phenotype is on the columns
stopifnot(Phenotype %in% names(db))
## The Genotype names is not empty
stopifnot(length(Genome)>0)
###############
#Some preprocessing
db <- df_replace_NA(db,Phenotype,"no")
#In principle, logistic regression works fine with categories in factors.
#db[[Phenotype]]<-as.numeric(db[[Phenotype]])-1

cat(c("Gene","Slope","\"Std.Error\"","p.val","R2","MI"),file = output_models,fill=T,sep = "\t")

for(Gene in Genome){
    if(length(unique(db[[Gene]]))<2){next}
    cur_db<-db[,c(Phenotype,Gene)]
    mI_categorical(cur_db[[Phenotype]],cur_db[[Gene]])->MI
    stop()
    glm(cur_db,family=binomial(link="logit"))->model
    summary(model)->sm
    output <-  as.vector(c(Gene = Gene,sm$coefficients[Gene,c("Estimate","Std. Error","Pr(>|t|)")],r2 = sm$r.squared,MI = MI))
    cat(output,"\n",file =output_models,sep = c("\t"),append = T,fill=F)
}
