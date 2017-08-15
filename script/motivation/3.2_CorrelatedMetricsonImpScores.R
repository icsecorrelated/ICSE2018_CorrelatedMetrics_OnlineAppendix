



# Motivation analysis
# Section 3.2 The impact of the number of correlated metrics on importance scores

# Import packages
library(DefectData)
library(rms)
library(Hmisc)
library(effsize)
library(car)
library(randomForest)
library(ggplot2)

# List of functions

# Remove correlated metrics
removeCorrelatedMetrics <-
    function(dataset,
             sw.metrics,
             defect,
             VarClus.threshold,
             VIF.threshold) {
        # Remove constant metrics, i.e., metrics that are a constant value (max value == min value)
        mitigated.sw.metrics <-
            sw.metrics[!apply(dataset[, sw.metrics], 2, function(metric)
                max(metric) == min(metric))]
        
        # Apply variable clustering with the threshold Spearman correlation of 0.7
        mitigated.sw.metrics <-
            applyVarClus(dataset, sw.metrics, VarClus.threshold)
        # VarClus 1st iteration
        # We find that there are 7 clusters of correlated metrics and 3 non-correlated metrics.
        # We perform a manual selection to select a representative metric for each cluster
        # Varclus 1st iteration
        # [1] "Cluster 1 : FOUT_avg, FOUT_max, FOUT_sum, MLOC_avg, MLOC_max, MLOC_sum, NBD_avg, NBD_max, NBD_sum, TLOC, VG_avg, VG_max, VG_sum"
        # We select TLOC from Cluster 1
        # [1] "Cluster 2 : NOF_avg, NOF_max, NOF_sum"
        # We select NOF_sum from Cluster 2
        # [1] "Cluster 3 : NOI, NOT"
        # We select NOT from Cluster 3
        # [1] "Cluster 4 : NOM_avg, NOM_max, NOM_sum"
        # We select NOM_sum from Cluster 4
        # [1] "Cluster 5 : NSF_avg, NSF_max, NSF_sum"
        # We select NSF_sum from Cluster 5
        # [1] "Cluster 6 : NSM_avg, NSM_max, NSM_sum"
        # We select NOF_sum from Cluster 6
        # [1] "Cluster 7 : PAR_max, PAR_sum"
        # We select PAR_sum from Cluster 7
        # Together with 3 non-correlated metrics (i.e., pre, ACD, PAR_avg), we have 10 metrics as follows:
        # pre ACD NOF_sum NOM_sum NOT NSF_sum NSM_max PAR_avg PAR_sum TLOC
        mitigated.sw.metrics <-
            c(
                'pre',
                'ACD',
                'NOF_sum',
                'NOM_sum',
                'NOT',
                'NSF_sum',
                'NSM_max',
                'PAR_avg',
                'PAR_sum',
                'TLOC'
            )
        
        # To ensure that all of the correlated metrics are removed, we reapply variable clustering
        mitigated.sw.metrics <-
            applyVarClus(dataset, mitigated.sw.metrics, VarClus.threshold)
        # VarClus 2nd iteration
        # We find that there are 1 clusters of correlated metrics and 7 non-correlated metrics.
        # We perform a manual selection to select a representative metric for each cluster
        # [1] "Cluster 1 : NOM_sum, PAR_sum, TLOC"
        # We select TLOC from Cluster 1
        # Together with 7 non-correlated metrics (i.e., pre, ACD, NOF_sum, NOT, NSF_sum, NSM_max, PAR_avg), we have 8 metrics as follows:
        # pre ACD NOF_sum NOT NSF_sum NSM_max PAR_avg TLOC
        mitigated.sw.metrics <-
            c('pre',
              'ACD',
              'NOF_sum',
              'NOT',
              'NSF_sum',
              'NSM_max',
              'PAR_avg',
              'TLOC')
        
        # Apply variance inflation factor analysis with the threshold value of 5
        mitigated.sw.metrics <-
            applyVIF(dataset, mitigated.sw.metrics, defect, VIF.threshold)
        
        return (mitigated.sw.metrics)
        
    }

# VarClus (Variable Clustering - see section 2.1)
applyVarClus <- function(dataset, sw.metrics, VarClus.threshold) {
    print('Variable Clustering : START')
    output <- {
        
    }
    # Apply variable clustering
    vc <-
        varclus( ~ .,
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

# VIF (Variance inflation factor - see Section 2.1)
applyVIF <- function(dataset,
                     sw.metrics,
                     defect,
                     VIF.threshold) {
    print('Variance Inflation Factor : START')
    glm.model <-
        glm(as.formula(paste(
            defect, '~', paste0(sw.metrics, collapse = '+')
        )),
        data = dataset,
        family = binomial())
    VIF <- rms::vif(glm.model)
    inter.correlated.metrics <- names(VIF[VIF >= VIF.threshold])
    mitigated.metrics <-
        sw.metrics[!(sw.metrics %in% inter.correlated.metrics)]
    
    print(paste0('nMetrics : ', length(sw.metrics)))
    print(paste0(
        'Mitigated metrics (',
        length(mitigated.metrics),
        '): ',
        paste0(mitigated.metrics, collapse = ', ')
    ))
    print(paste0(
        'Inter-correlated metrics (',
        length(inter.correlated.metrics),
        '): ',
        paste0(inter.correlated.metrics, collapse = ', ')
    ))
    
    # return original software metrics if all of the metrics are removed by VIF
    if (length(mitigated.metrics) == 0) {
        return(sw.metrics[sw.metrics %in% sw.metrics])
    }
    
    # reapply VIF to ensure that there is no presence of inter-correlated metrics
    if (length(inter.correlated.metrics) != 0) {
        mitigated.metrics <-
            VIF(dataset, mitigated.metrics, defect, VIF.threshold)
    }
    
    print('Variance Inflation Factor : END')
    return(sw.metrics[sw.metrics %in% mitigated.metrics])
    
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
        cliff.delta.list[order(-abs(cliff.delta.list$estimate)), ]
    
    return(cliff.delta.list)
    
}


# Main

# Init thresthold values
VarClus.threshold <- 0.7
VIF.threshold <- 5

# Load eclipse-2.0 dataset
Data <- loadData('eclipse-2.0')
# Retrieve data, name of defect-proneness and software metrics
dataset <- Data$data
defect <- Data$dep
sw.metrics <- Data$indep

# Remove correlated metrics
mitigated.sw.metrics <-
    removeCorrelatedMetrics(dataset, sw.metrics, defect, VarClus.threshold, VIF.threshold)

# Get cliff's delta estimate for each metric, i.e., determine the relationship of each metric with defect-proneness
cliff.delta.list <- getCliffDelta(dataset, sw.metrics, defect)
# We find that TLOC is the metric that shares the strongest relationship with defect-proneness

# Find metrics (correlated.metrics) that are correlated with the metric that shares the strongest relationship with defect-proneness (TLOC)
metrics.correlation <-
    applyVarClus(dataset, sw.metrics, VarClus.threshold)
correlated.metrics <-
    as.character(metrics.correlation$CorrelatedMetrics[grepl('TLOC', metrics.correlation$CorrelatedMetrics)])
correlated.metrics <- strsplit(correlated.metrics, '[,]')[[1]]
correlated.metrics <-
    as.character(cliff.delta.list$name[cliff.delta.list$name %in% correlated.metrics])[-1]

# Generate orders of metrics (see Section 3.2)
mitigated.sw.metrics <-
    mitigated.sw.metrics[mitigated.sw.metrics != 'TLOC']
orders <- list()

# First order: TLOC + Control metrics (Baseline model)
orders[[1]] <- c('TLOC', mitigated.sw.metrics)
for (i in 1:length(correlated.metrics)) {
    # Following orders: Metrics that are correlated with TLOC + TLOC + Control metrics (Mitigated metrics except TLOC)
    orders[[i + 1]] <-
        c(correlated.metrics[1:i], 'TLOC', mitigated.sw.metrics)
}

# Examine the importance scores of TLOC for each order using:
# - ANOVA Type-I
# - ANOVA Type-II/III (LR option)
# - Gini Importance
# - Permutation Importance
TLOC.importance.scores <-
    list(
        ANOVA.TypeI = {
            
        },
        ANOVA.TypeII.III.LR = {
            
        },
        Gini.Importance = {
            
        },
        Permutation.Importance = {
            
        }
    )


# Generate a bootstrap sample
set.seed(1)
train.data.index <- sample(nrow(dataset), replace = TRUE)
for (i in 1:length(orders)) {
    print(paste0('Model formula ', i))
    print(paste0(defect, ' ~ ', paste0(orders[[i]], collapse = ' + ')))
    
    
    # Build a logistic regression model
    options(contrasts = c("contr.sum", "contr.poly"))
    lr.model <-
        glm(as.formula(paste(
            defect, '~', paste0(orders[[i]], collapse = '+')
        )),
        data = dataset[train.data.index, ],
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
        anova(lr.model, test = c('Chisq'))[2][-1,] # Get each metric sum square
    anova.typeI <-
        anova.typeI / sum(anova.typeI) * 100 # Normalize to percentage
    TLOC.importance.scores$ANOVA.TypeI <-
        rbind(TLOC.importance.scores$ANOVA.TypeI, anova.typeI[match('TLOC', orders[[i]])])
    
    # ANOVA Type-II/III (LR option)
    anova.typeII.III.LR <-
        Anova(lr.model, test.statistic = 'LR')[, 1] # Get each metric sum square
    anova.typeII.III.LR <-
        anova.typeII.III.LR / sum(anova.typeII.III.LR) * 100 # Normalize to percentage
    TLOC.importance.scores$ANOVA.TypeII.III.LR <-
        rbind(TLOC.importance.scores$ANOVA.TypeII.III.LR,
              anova.typeII.III.LR[match('TLOC', orders[[i]])])
    
    # Gini Importance
    gini.importance <-
        importance(rf.model, type = 2)[, 1] # Get each metric sum square
    gini.importance <-
        gini.importance / sum(gini.importance) * 100 # Normalize to percentage
    TLOC.importance.scores$Gini.Importance <-
        rbind(TLOC.importance.scores$Gini.Importance,
              gini.importance[match('TLOC', orders[[i]])])
    
    # Permutation Importance
    permutation.importance <-
        importance(rf.model, type = 1, scale = F)[, 1] # Get each metric sum square
    permutation.importance <-
        permutation.importance / sum(permutation.importance) * 100 # Normalize to percentage
    TLOC.importance.scores$Permutation.Importance <-
        rbind(TLOC.importance.scores$Permutation.Importance,
              permutation.importance[match('TLOC', orders[[i]])])
    
}

# Normalize TLOC scores for each model with the score of the baseline model
TLOC.importance.scores <-
    lapply(TLOC.importance.scores, function(x)
        x / x[1] * 100)

# Preprocess to plot a graph
TLOC.plot.data <- melt(TLOC.importance.scores)
TLOC.plot.data$Var1 <- TLOC.plot.data$Var1 - 1
TLOC.plot.data$L1 <- factor(TLOC.plot.data$L1)
levels(TLOC.plot.data$L1) <-
    c('ANOVA Type I',
      'ANOVA Type II/III (LR)',
      'Gini Importance',
      'Permutation Importance')

# Plot a graph
ggplot(data = TLOC.plot.data,
       aes(x = Var1, y = value, fill = L1)) +
    geom_bar(colour = "black",
             position = "dodge",
             stat = "identity") +
    labs(x = "Number of correlated metrics in a model", y = "Relative difference of importance score (%)") +
    scale_fill_manual(name = "Model Interpretation\nTechnique",
                      values = c("#d73027",
                                 "#fee090",
                                 "#91bfdb",
                                 "#4575b4")) +
    coord_cartesian(ylim = c(0, 100)) +
    facet_wrap(~ L1, ncol = 2) +
    theme(legend.position = "none")

# Export a graph
ggsave(filename = '3.2_plot.pdf',
       width = 6.5,
       height = 4)
