# R and Python package for implementing ensemble methods of rank-based trees



## Getting Started

### Dependencies

#### For R package:

* Prerequisites: caret, data.table, data.tree, DiagrammeR, dplyr, gbm, ggplot2,
                 parallel, randomForestSRC
* R (>= 4.0.0)

### Installing

#### For R package:

* Download the R package from /download/EnsembleRankTrees_0.1.0.tar.gz
* install.packages("EnsembleRankTrees")


### Examples:

* Build a Random Rank Forest:
```
library(EnsembleRankTrees)
# load the example data:
data(liver_data)
model <- rrf(liver_data$train_x_liver, liver_data$train_y_liver)
```
* Build a Boosting with LogitBoost Cost model:
```
library(EnsembleRankTrees)
# load the example data:
data(liver_data)
model <- blc(liver_data$train_x_liver, liver_data$train_y_liver)
```


## Authors

Ruijie Yin (ruijieyin428@gmail.com)



## Version History

* 0.1.1
    * Added decision rule-based methods as a reference to explain the results

* 0.1.0
    * Initial Release

## License

This project is licensed under the GPL 3.0 License - see the LICENSE.md file for details



