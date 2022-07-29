function (response, predictor, controls, cases, density.controls, 
    density.cases, levels = base::levels(as.factor(response)), 
    percent = FALSE, na.rm = TRUE, direction = c("auto", "<", 
        ">"), algorithm = 6, quiet = FALSE, smooth = FALSE, auc = TRUE, 
    ci = FALSE, plot = FALSE, smooth.method = "binormal", smooth.n = 512, 
    ci.method = NULL, density = NULL, ...) 
{
    direction <- match.arg(direction)
    if (!missing(response) && !is.null(response) && !missing(predictor) && 
        !is.null(predictor)) {
        if ((!missing(cases) && !is.null(cases)) || (!missing(controls) && 
            !is.null(controls))) {
            stop("'cases/controls' argument incompatible with 'response/predictor'.")
        }
        if ((!missing(density.cases) && !is.null(density.cases)) || 
            (!missing(density.controls) && !is.null(density.controls))) {
            stop("'density.*' arguments incompatible with 'response/predictor'.")
        }
        original.predictor <- predictor
        original.response <- response
        if (missing(levels)) {
            if (length(levels) > 2) {
                warning("'response' has more than two levels. Consider setting 'levels' explicitly or using 'multiclass.roc' instead")
                levels <- levels[1:2]
            }
            else if (length(levels) < 2) {
                stop("'response' must have two levels")
            }
            ifelse(quiet, invisible, message)(sprintf("Setting levels: control = %s, case = %s", 
                levels[1], levels[2]))
        }
        else if (length(levels) != 2) {
            stop("'levels' argument must have length 2")
        }
        if (!is.numeric(predictor)) {
            if (is.ordered(predictor)) {
                predictor <- tryCatch({
                  as.numeric(as.character(predictor))
                }, warning = function(warn) {
                  warning("Ordered predictor converted to numeric vector. Threshold values will not correspond to values in predictor.")
                  return(as.numeric(predictor))
                })
            }
            else {
                stop("Predictor must be numeric or ordered.")
            }
        }
        if (is.matrix(predictor)) {
            warning("Deprecated use a matrix as predictor. Unexpected results may be produced, please pass a numeric vector.")
        }
        if (is.matrix(response)) {
            warning("Deprecated use a matrix as response. Unexpected results may be produced, please pass a vector or factor.")
        }
        if (length(predictor) != length(response)) {
            stop("Response and predictor must be vectors of the same length.")
        }
        if (na.rm) {
            nas <- is.na(response) | is.na(predictor)
            if (any(nas)) {
                na.action <- grep(TRUE, nas)
                class(na.action) <- "omit"
                response <- response[!nas]
                attr(response, "na.action") <- na.action
                predictor <- predictor[!nas]
                attr(predictor, "na.action") <- na.action
            }
        }
        else if (any(is.na(c(predictor[response == levels[1]], 
            predictor[response == levels[2]], response)))) 
            return(NA)
        splitted <- split(predictor, response)
        controls <- splitted[[as.character(levels[1])]]
        if (length(controls) == 0) 
            stop("No control observation.")
        cases <- splitted[[as.character(levels[2])]]
        if (length(cases) == 0) 
            stop("No case observation.")
        patients.in.levels <- response %in% levels
        if (!all(patients.in.levels)) {
            response <- response[patients.in.levels]
            predictor <- predictor[patients.in.levels]
        }
        if (any(which <- is.infinite(predictor))) {
            warning("Infinite values(s) in predictor, cannot build a valid ROC curve. NaN returned instead.")
            return(NaN)
        }
    }
    else if (!missing(cases) && !is.null(cases) && !missing(controls) && 
        !is.null(controls)) {
        if ((!missing(density.cases) && !is.null(density.cases)) || 
            (!missing(density.controls) && !is.null(density.controls))) {
            stop("'density.*' arguments incompatible with 'response/predictor'.")
        }
        if (na.rm) {
            if (any(is.na(controls))) 
                controls <- na.omit(controls)
            if (any(is.na(cases))) 
                cases <- na.omit(cases)
        }
        else if (any(is.na(c(controls, cases)))) 
            return(NA)
        if (length(controls) == 0) 
            stop("No control observation.")
        if (length(cases) == 0) 
            stop("No case observation.")
        if (is.ordered(cases)) {
            if (is.ordered(controls)) {
                if (identical(attr(cases, "levels"), attr(controls, 
                  "levels"))) {
                  original.predictor <- ordered(c(as.character(cases), 
                    as.character(controls)), levels = attr(controls, 
                    "levels"))
                  predictor <- as.numeric(original.predictor)
                  controls <- as.numeric(controls)
                  cases <- as.numeric(cases)
                }
                else {
                  stop("Levels of cases and controls differ.")
                }
            }
            else {
                stop("Cases are of ordered type but controls are not.")
            }
        }
        else if (is.numeric(cases)) {
            if (is.numeric(controls)) {
                predictor <- c(controls, cases)
                original.predictor <- predictor
            }
            else {
                stop("Cases are of numeric type but controls are not.")
            }
        }
        else {
            stop("Cases and controls must be numeric or ordered.")
        }
        if (any(which <- is.infinite(predictor))) {
            warning("Infinite values(s) in predictor, cannot build a valid ROC curve. NaN returned instead.")
            return(NaN)
        }
        response <- c(rep(0, length(controls)), rep(1, length(cases)))
        original.response <- response
        levels <- c(0, 1)
    }
    else if (!missing(density.cases) && !is.null(density.cases) && 
        !missing(density.controls) && !is.null(density.controls)) {
        if (!is.numeric(density.cases) || !is.numeric(density.controls)) 
            stop("'density.cases' and 'density.controls' must be numeric values of density (over the y axis).")
        if (direction == "auto") 
            dir <- "<"
        else dir <- direction
        smooth.roc <- smooth.roc.density(density.controls = density.controls, 
            density.cases = density.cases, percent = percent, 
            direction = dir)
        class(smooth.roc) <- "smooth.roc"
        smooth.roc <- sort(smooth.roc)
        smooth.roc$specificities <- c(0, as.vector(smooth.roc$specificities), 
            ifelse(percent, 100, 1))
        smooth.roc$sensitivities <- c(ifelse(percent, 100, 1), 
            as.vector(smooth.roc$sensitivities), 0)
        smooth.roc$percent <- percent
        smooth.roc$direction <- direction
        smooth.roc$call <- match.call()
        if (auc) {
            smooth.roc$auc <- auc(smooth.roc, ...)
            if (direction == "auto" && smooth.roc$auc < roc.utils.min.partial.auc.auc(smooth.roc$auc)) {
                smooth.roc <- roc.default(density.controls = density.controls, 
                  density.cases = density.cases, levels = levels, 
                  percent = percent, direction = ">", auc = auc, 
                  ci = ci, plot = plot, ...)
                smooth.roc$call <- match.call()
                return(smooth.roc)
            }
        }
        if (ci) 
            warning("CI can not be computed with densities.")
        if (plot) 
            plot.roc(smooth.roc, ...)
        return(smooth.roc)
    }
    else {
        stop("No valid data provided.")
    }
    if (direction == "auto" && median(controls) <= median(cases)) {
        direction <- "<"
        ifelse(quiet, invisible, message)("Setting direction: controls < cases")
    }
    else if (direction == "auto" && median(controls) > median(cases)) {
        direction <- ">"
        ifelse(quiet, invisible, message)("Setting direction: controls > cases")
    }
    if (smooth) {
        if (missing(density.controls)) 
            density.controls <- density
        if (missing(density.cases)) 
            density.cases <- density
    }
    if (isTRUE(algorithm == 6)) {
        if (is.numeric(predictor)) {
            algorithm <- 2
        }
        else {
            algorithm <- 3
        }
    }
    else if (isTRUE(algorithm == 0)) {
        load.suggested.package("microbenchmark")
        cat("Starting benchmark of algorithms 2 and 3, 10 iterations...\n")
        thresholds <- roc.utils.thresholds(c(controls, cases), 
            direction)
        benchmark <- microbenchmark::microbenchmark(`2` = roc.utils.perfs.all.fast(thresholds = thresholds, 
            controls = controls, cases = cases, direction = direction), 
            `3` = rocUtilsPerfsAllC(thresholds = thresholds, 
                controls = controls, cases = cases, direction = direction), 
            times = 10)
        print(summary(benchmark))
        if (any(is.na(benchmark))) {
            warning("Microbenchmark returned NA. Using default algorithm 1.")
            algorithm <- 2
        }
        algorithm <- as.integer(names(which.min(tapply(benchmark$time, 
            benchmark$expr, sum))))
        cat(sprintf("Selecting algorithm %s.\n", algorithm))
    }
    else if (isTRUE(algorithm == 5)) {
        thresholds <- length(roc.utils.thresholds(c(controls, 
            cases), direction))
        if (thresholds > 55) {
            algorithm <- 2
        }
        else {
            algorithm <- 3
        }
    }
    if (isTRUE(algorithm == 2)) {
        fun.sesp <- roc.utils.perfs.all.fast
    }
    else if (isTRUE(algorithm == 3)) {
        fun.sesp <- rocUtilsPerfsAllC
    }
    else if (isTRUE(algorithm == 1)) {
        fun.sesp <- roc.utils.perfs.all.safe
    }
    else if (isTRUE(algorithm == 4)) {
        fun.sesp <- roc.utils.perfs.all.test
    }
    else {
        stop("Unknown algorithm (must be 0, 1, 2, 3, 4 or 5).")
    }
    roc <- roc.cc.nochecks(controls, cases, percent = percent, 
        direction = direction, fun.sesp = fun.sesp, smooth = smooth, 
        density.cases = density.cases, density.controls = density.controls, 
        smooth.method = smooth.method, smooth.n = smooth.n, auc, 
        ...)
    roc$call <- match.call()
    if (smooth) {
        attr(roc, "roc")$call <- roc$call
        attr(roc, "roc")$original.predictor <- original.predictor
        attr(roc, "roc")$original.response <- original.response
        attr(roc, "roc")$predictor <- predictor
        attr(roc, "roc")$response <- response
        attr(roc, "roc")$levels <- levels
    }
    roc$original.predictor <- original.predictor
    roc$original.response <- original.response
    roc$predictor <- predictor
    roc$response <- response
    roc$levels <- levels
    if (auc) {
        attr(roc$auc, "roc") <- roc
    }
    if (ci) 
        roc$ci <- ci(roc, method = ci.method, ...)
    if (plot) 
        plot.roc(roc, ...)
    return(roc)
}
<bytecode: 0x55d60ccd2ec8>
<environment: namespace:pROC>
