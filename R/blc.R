require(randomForestSRC)
require(gbm)
require(caret)



# build a Boosting with LogitBoost Cost model

blc <- function (train_x, train_y = NULL, n.tree = 500, n.minobsinnode = 10) {
  options(warn = -1)

  # function: check if Check if a dataset has been row standardized (sample-wise)
  is_row_standardized <- function(data) {
    row_means <- rowMeans(data)
    row_vars <- apply(data, 1, var)
    all(abs(row_means) < 1e-10) && all(abs(row_vars - 1) < 1e-10)
  }


  # function: convert.genepairs
  # convert the raw expression values into gene pairs
  convert.genepairs <- function(input_dat = train_x) {
    # note: columns should be genes, rows are samples
    # the last column are the class labels
    # get gene pair data:
    input_expr = input_dat
    # get the number of columns
    n_col = dim(input_expr)[2]
    converted_input_expr = NULL
    for (i in 1:(n_col - 1)) {
      for (j in (i + 1):n_col) {
        converted_input_expr = cbind(converted_input_expr, as.numeric(input_expr[, i] < input_expr[, j]))
      }
    }
    data.list <- list("data.converted" = converted_input_expr)
    return(data.list)
  }

  if (!is.data.frame(train_x) && !is.matrix(train_x)) {
    # sanity check 1: if the input training data is a data frame or a matrix:
    stop("Please convert the input dataset to a data.frame or matrix object")
  } else if (is_row_standardized(train_x)) {
    # sanity check 2: if the input training data has been row standardized:
    stop("The dataset has been row standardized. Please do not row standardize the data.")
  } else {
    # remove symbols and special characters in column names, and replace them with "_"
    gsub("[[:punct:]]", "_", colnames(train_x))
    # convert input_dat to gene pairs:
    dat = convert.genepairs(train_x)$data.converted
    # merge expr and labels:
    data_train = cbind.data.frame(dat, as.factor(train_y))
    colnames(data_train)[dim(data_train)[2]] <- "subtype"
    blc.model = gbm(subtype ~.,
                    data = data_train,
                    distribution = "multinomial",
                    shrinkage = .01,
                    n.minobsinnode = n.minobsinnode,
                    n.trees = n.tree)
  }
  model.list <- list("model" = blc.model)
  return(blc.model)
}















