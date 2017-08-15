
# Install required packages

package.list <- c(
    'rms',
    'Hmisc',
    'effsize',
    'car',
    'randomForest',
    'ggplot2',
    'caret',
    'pROC',
    'mccr',
    'devtools'
)

installed.package.list <- installed.packages()[, "Package"]
new.packages <-
    package.list[!(package.list %in% installed.package.list)]
if (length(new.packages))
    install.packages(new.packages, repos = "http://cran.rstudio.com/",  dependencies =
                         TRUE)
if (!('DefectData' %in% installed.package.list)) {
    devtools::install_github("klainfo/DefectData")
}

print('All of the required packages have been installed.')