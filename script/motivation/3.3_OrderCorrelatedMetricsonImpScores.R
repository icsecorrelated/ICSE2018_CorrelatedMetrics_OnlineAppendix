


# Motivation analysis
# Section 3.3 The impact of the order of correlated metrics on importance scores

# Import packages
library(DefectData)
library(rms)
library(Hmisc)
library(effsize)
library(car)
library(randomForest)
library(ggplot2)
library(caret)
library(pROC)
library(mccr)

# List of functions

# VarClus (Variable Clustering - see section 2.1)
applyVarClus <- function(dataset, sw.metrics, VarClus.threshold) {
    print('Variable Clustering : START')
    output <- {
        
    }
    # Apply variable clustering
    vc <-
        varclus(~ .,
                similarity = 'spearman',
                data = dataset[, sw.metrics],
                trans = "abs")
    
    # Apply cutoff threshold
    var.clustered <-
        cutree(vc$hclust, h = (1 - VarClus.threshold))
    
    # Get surviving metrics
    non.correlated.metrics <-
        names(var.clustered)[var.clustered %in% names(table(var.clustered)[table(var.clustered) == 1])]
    print(paste0(
        length(non.correlated.metrics),
        ' non-correlated metrics : ',
        paste0(non.correlated.metrics, collapse = ', ')
    ))
    
    # Get correlated clusters index
    correlated.clusters.index <-
        names(table(var.clustered)[table(var.clustered) > 1])
    
    print(paste0((length(sw.metrics) - length(non.correlated.metrics)),
                 ' correlated metrics from ',
                 length(correlated.clusters.index),
                 ' clusters'
    ))
    
    # For each cluster of correlated metrics, print out for manual selection
    cluster.count <- 1
    for (cluster.index in correlated.clusters.index) {
        print(paste0(
            'Cluster ',
            cluster.count,
            ' : ',
            paste0(names(var.clustered)[var.clustered == cluster.index], collapse = ', ')
        ))
        output <- rbind(output,
                        data.frame(
                            CID = cluster.count,
                            CorrelatedMetrics = paste0(names(var.clustered)[var.clustered == cluster.index], collapse = ',')
                        ))
        cluster.count <- cluster.count + 1
    }
    print('Variable Clustering : END')
    return (output)
}

# Cliff's Delta Estimates
getCliffDelta <- function(dataset, sw.metrics, defect) {
    cliff.delta.list <- {
        
    }
    for (index in 1:length(sw.metrics)) {
        cd <-
            cliff.delta(dataset[dataset[, defect] == TRUE, sw.metrics[index]], dataset[dataset[, defect] == FALSE, sw.metrics[index]])
        cliff.delta.list <- rbind(
            cliff.delta.list,
            data.frame(
                name = sw.metrics[index],
                estimate = cd$estimate,
                magnitude = cd$magnitude
            )
        )
    }
    
    cliff.delta.list <-
        cliff.delta.list[order(-abs(cliff.delta.list$estimate)),]
    
    return(cliff.delta.list)
    
}

getPerformanceLR <-
    function(lr.model,
             testing.data,
             sw.metrics,
             defect) {
        testing.data[, defect] <- as.factor(testing.data[, defect])
        
        lr.testing.probs <-
            predict(lr.model, newdata = testing.data[, c(sw.metrics, defect)], type = "response")
        lr.testing.probs.factor <-
            factor(ifelse(lr.testing.probs >= 0.5, "TRUE", "FALSE"))
        
        prob01 <- lr.testing.probs.factor
        levels(prob01) <- c(1, 0)
        test01 <- testing.data[, defect]
        levels(test01) <- c(1, 0)
        lr.mcc <- mccr(test01, prob01)
        lr.precision <-
            posPredValue(lr.testing.probs.factor, testing.data[, defect])
        lr.recall <-
            sensitivity(lr.testing.probs.factor, testing.data[, defect])
        lr.f1 <-
            (2 * lr.precision * lr.recall) / (lr.precision + lr.recall)
        lr.roc <- roc(testing.data[, defect], lr.testing.probs)
        lr.auc <- as.numeric(lr.roc$auc)
        lr.misclassification <-
            mean(lr.testing.probs.factor != testing.data[, defect])
        
        return(
            data.frame(
                precision = lr.precision,
                recall = lr.recall,
                f1 = lr.f1,
                auc = lr.auc,
                misclassification = lr.misclassification,
                mcc = lr.mcc
            )
        )
    }

getPerformanceRF <-
    function(rf.model,
             testing.data,
             sw.metrics,
             defect) {
        testing.data[, defect] <- as.factor(testing.data[, defect])
        
        rf.testing.probs <-
            predict(rf.model, newdata = testing.data[, c(sw.metrics, defect)], type = 'prob')[, 1]
        rf.testing.probs.factor <-
            factor(ifelse(rf.testing.probs >= 0.5, "TRUE", "FALSE"))
        
        prob01 <- rf.testing.probs.factor
        levels(prob01) <- c(1, 0)
        test01 <- testing.data[, defect]
        levels(test01) <- c(1, 0)
        rf.mcc <- mccr(test01, prob01)
        rf.precision <-
            posPredValue(rf.testing.probs.factor, testing.data[, defect])
        rf.recall <-
            sensitivity(rf.testing.probs.factor, testing.data[, defect])
        rf.f1 <-
            (2 * rf.precision * rf.recall) / (rf.precision + rf.recall)
        rf.roc <- roc(testing.data[, defect], rf.testing.probs)
        rf.auc <- as.numeric(rf.roc$auc)
        rf.misclassification <-
            mean(rf.testing.probs.factor != testing.data[, defect])
        
        return(
            data.frame(
                precision = rf.precision,
                recall = rf.recall,
                f1 = rf.f1,
                auc = rf.auc,
                misclassification = rf.misclassification,
                mcc = rf.mcc
            )
        )
        
    }


# Main

# Init thresthold values
VarClus.threshold <- 0.7

# Load eclipse-2.0 dataset
Data <- loadData('eclipse-2.0')
# Retrieve data, name of defect-proneness and software metrics
dataset <- Data$data
defect <- Data$dep
sw.metrics <- Data$indep

# Get cliff's delta estimate for each metric, i.e., determine the relationship of each metric with defect-proneness
cliff.delta.list <- getCliffDelta(dataset, sw.metrics, defect)
# We find that TLOC is the metric that shares the strongest relationship with defect-proneness

# Find metrics (correlated.metrics) that are correlated with the metric that shares the strongest relationship with defect-proneness (TLOC)
metrics.correlation <-
    applyVarClus(dataset, sw.metrics, VarClus.threshold)
correlated.metrics <-
    as.character(metrics.correlation$CorrelatedMetrics[grepl('TLOC', metrics.correlation$CorrelatedMetrics)])
correlated.metrics <- strsplit(correlated.metrics, '[,]')[[1]]

# We randomly select 5 metrics, i.e., TLOC, MLOC_sum, FOUT_sum, VG_sum, and NBD_sum, in order to simplify the demonstration
selected.metrics <- correlated.metrics[c(10, 6, 3, 13, 9)]

# From selected metrics, we generate 2 orders of metrics.
orders <- list(# - Order1: TLOC appears at the first position
    order.1 = selected.metrics,
    # - Order2: TLOC appears at the last position
    order.2 = selected.metrics[c(2:5, 1)])

# Generate a bootstrap sample
set.seed(1)
train.data.index <- sample(nrow(dataset), replace = TRUE)

outputList <- list()
for (i in 1:length(orders)) {
    print(paste0('Model formula ', i))
    print(paste0(defect, ' ~ ', paste0(orders[[i]], collapse = ' + ')))
    
    
    # Build a logistic regression model
    options(contrasts = c("contr.sum", "contr.poly"))
    lr.model <-
        glm(as.formula(paste(
            defect, '~', paste0(orders[[i]], collapse = '+')
        )),
        data = dataset[train.data.index,],
        family = binomial())
    
    # Build a random forest model (Please note that randomForest might be difference for each execution)
    rf.model <- randomForest(
        x = dataset[train.data.index, orders[[i]]],
        # In order to let a random forest model to perform classification instead of regression, we have to convert data type to factor
        y = as.factor(dataset[train.data.index, defect]),
        ntree = 100,
        importance = TRUE,
        keep.forest = TRUE
    )
    
    # ANOVA Type-I
    anova.typeI <-
        anova(lr.model, test = c('Chisq'))[2][-1, ] # Get each metric sum square
    anova.typeI <-
        anova.typeI / sum(anova.typeI) * 100 # Normalize to percentage
    names(anova.typeI) <- orders[[i]] # Label with metric names
    anova.typeI <- anova.typeI[orders[[1]]]  # Sort to Model 1 order
    
    # Gini Importance
    gini.importance <-
        importance(rf.model, type = 2)[, 1] # Get each metric sum square
    gini.importance <-
        gini.importance / sum(gini.importance) * 100 # Normalize to percentage
    gini.importance <-
        gini.importance[orders[[1]]] # Sort to Model 1 order
    
    outputList[[i]] <- list(
        lr.interpretation = anova.typeI,
        lr.performance = getPerformanceLR(lr.model, dataset[-train.data.index, ], orders[[i]], defect),
        rf.interpretation = gini.importance,
        rf.performance = getPerformanceRF(rf.model, dataset[-train.data.index, ], orders[[i]], defect)
    )
    
}
