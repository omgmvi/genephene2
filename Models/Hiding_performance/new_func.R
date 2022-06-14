## Function to produce a confusion matrix and the metrics

## Reinventing the wheel

## This function will tackle a single problem, given the outcome of a logistic model or similar
## where the output is a vector of probabilities? or at least a number between 0 and 1
## enter a cutoff /threshold value and calculate the classes.

# I am not sure if this is exactly what I want: I thought that those algorithm fits the parameters so you do not need to select a threshold, or in other way, it is fixed that the threshold would be 0.5 so the fit is according to choosing this value later (this seems not happening in glm)

## Full goal:

## Evaluate the accuracy, precision, recall, TPR, TFR, specificity, sensitivity, F1, ... given a predicted and a reference
## Find the best cutoff value for a given parameter, e.g. precision.
## Calculate the ROC AUC and the PR AUC as summarizing


## Therefore, to start, 

# Bibliography -


#### How to perform a logistic regression using glm: https://www.r-bloggers.com/2015/09/how-to-perform-a-logistic-regression-in-r/
###  basically : glm(formula,data, family = binomial(link("logit")):

#   glm, returns fitted.values with the results of the sigmoid function 1/(1+exp(-(a*x+b)))
# since it is a binomial - I only need to select a class for reference

#### Use of 0.5 as threshold : https://www.yourdatateacher.com/2021/06/14/are-you-still-using-0-5-as-a-threshold/
# Good explanation on how to choose the threshold based on the properties of the ROC AUC

#### Tuning the 0.5 threshold for ROC and PR curves : https://machinelearningmastery.com/threshold-moving-for-imbalanced-classification/
# It makes the point explicitly for imbalanced classification (yeah, our problem!). 
# Describes a procedure to choose the threshold (decision threshold) as a sort of hyper parameters to later use in the validation sets 

#### ROC CURVES and PR curves + metrics
#### https://machinelearningmastery.com/roc-curves-and-precision-recall-curves-for-imbalanced-classification/
#### https://machinelearningmastery.com/roc-curves-and-precision-recall-curves-for-classification-in-python/
#### https://scikit-learn.org/stable/auto_examples/model_selection/plot_precision_recall.html
#### https://machinelearningmastery.com/threshold-moving-for-imbalanced-classification/
####
####
# FUNCTIONS #

# Convert a factor in its numeric representation starting from 0
#This function could be mess with a function that has numbers as levels and one would like to make a number vector - that would be as.numeric(as.character(fac))
factor2num <- function(fac){
    if(is.factor(fac)){as.numeric(fac)-1}else{stop("input is not a factor")}
}

#The reverse, pass a vector of nums to a factor given the levels of another factor - I feel this is dangerouse when the numbers do not coincide
num2factor <- function(num,fac){
    if(is.factor(fac)){factor(num,labels = levels(fac))}else{stop("input is not a factor")}
}


# Given a vector of 'probabilities' or [0,1] and a cutoff/ decision threshold value, produce a vector of predicted values
pass_threshold <- function(probs,thresh,fac=NULL){
    tmp <- as.numeric(probs>=thresh)
    if(fac != NULL){
        tmp <- num2factor(tmp,fac)
    }
    tmp
}

# compare two vectors, namely predicted, real and a reference class to detect TP, FP, TN, FN
confusion_matrix(pred,supervised,refclass){
    conf_matrix <- list(TP = 0,TN = 0,FP = 0,FN = 0)
    if(!any(pred %in% refclass) && any(supervised %in% refclass)){stop("Reference class not found in predicted values")}
    if(nrows(pred)==nrows(supervised)){
        is_ref <- supervised == refclass
        is_eq  <- predicted == supervised
        # 4cases
        conf_matrix$TP <- is_ref & is_eq
        conf_matrix$TN <- !is_ref & is_eq
        conf_matrix$FP <- !(is_ref |  is_eq)
        conf_matrix$FN <- is_ref & !is_eq
    }
}
